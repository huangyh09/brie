"""computes correlation between pseudotime and BRIE WX (prediction)

This script computes pseudotime for each cell in given cells'directory, then
runs BRIE analysis on each cell and compute the correlation coefficient
(Pearson's r) between BRIE predictor WX (with W the average of learned weights
and X corresponding cell's features) and pseudotime for each gene. It also
performs a pseudotime brie analysis if --brie-pseudotime-path is provided with
the right path.
This script only support the two isoforms per gene scenario.
To see how to use this program, run `python3 pseudotime_gene_correlation.py -h`
BRIE paper can be consulted here:
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1248-5
"""

import os # to navigate and create directories
import argparse # to parse script arguments
import numpy as np # required for pseudotime analysis
import pandas as pd # required for pseudotime analysis
import scanpy.api as sc # to perform pseudotime analysis
import scipy.stats as st # to compute Pearson's r correlation coefficient
import subprocess as sub # to run external programs
import csv # to write and read csv files
import pseudotime_auxiliary # pseudotime script to run main pseudotime brie
from utils.gtf_utils import loadgene # get list of Gene objects from gff/gtf
from math import exp # for logistic

def logistic(x):
    return exp(x)/(1+exp(x))

def parse_arguments():
    """ parse arguments of this script

    Returns
    -------
    args: Namespace
        Namespace of arguments (object with key-value pairs).
    """
    
    # default will be displayed in help thanks to formater_class:
    parser = argparse.ArgumentParser(description='computes correlation between'
                                     + ' pseudotime and BRIE WX (prediction)',
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # mandatory arguments
    parser.add_argument('output_dir',
                        help='output directory.')
    
    parser.add_argument("sam_dir",
                        help="Directory where are stored each sorted and "
                        "indexed bam/sam files. Files are assumed to be unique,"
                        ", sorted and with name exact name "
                        "'$(cell_id).sorted.[b|s]am'")

    parser.add_argument('annotation_file',
                        help='annotation gtf or gff3 file (where genes and '
                        'transcripts to consider are stored).')

    parser.add_argument('factor_file',
                        help='Features in csv.gz file to predict isoform '
                        'expression.')

        
    # optionnal arguments
    parser.add_argument('--brie-arguments', default="",
                        help='Other BRIE arguments given inside a string '
                        '(written as you would do calling BRIE script). Do not '
                        'forget to surround string arguments with quotes.'
                        '[beware: possible issue with long arguments, use short'
                        ' ones to be sure, or put spaces surrounding "=" sign]')

    parser.add_argument('--overwrite-brie', action='store_true',
                        help='If this argument is present, already computed '
                        'brie results in output_dir will be erased, else brie '
                        'computation will not be performed.')

    parser.add_argument('--counts-dir',
                        help="Directory where are stored precomputed counts' "
                        "files. Counts files are assumed to be unique for each "
                        "cell and with exact name '$(cell_id).gene_counts'. If"
                        " not provided, a directory will be created in output "
                        "directory and matrices computed and stored there. If "
                        "this directory is empty or do not exist, it is "
                        "interpreted as target directory where to store counts"
                        "to compute. "
                        "Ignored if --matrices-of-counts is provided.")

    parser.add_argument('--matrix-of-counts',
                        help="matrix of counts (gene expression level) csv "
                        "file.")

    #parser.add_argument('--brie-pseudotime-path',
    #                    help="path to brie_pseudotime.py script.")

    return parser.parse_args()

def counts_from_sam(sam_file, output_dir, annotation_file):
    """create file of counts (gene expression level) with htseq-count.

    Sam or bam file is assumed to be sorted and with exact name
    '$(cell_id).sorted.[b|s]am'. Output file will have name
    '$(cell_id).gene_counts'.

    Parameters
    ----------
    sam_file: string
        sam or bam file.
     output_dir: string
        output directory where are stored gene expression counts as
        '$(cell_id).gene_counts'.
    annotation_file: string
        genome annotation file in ggf or gff3 format.

    Returns
    -------
    int : 0
    """
    extension = '' # nature of sam_file: sam or bam

    if sam_file[-11:] == ".sorted.bam":
        extension = 'bam'
    elif sam_file[-11:] == ".sorted.sam":
        extension = 'sam'

    if extension: # if sam_file have a suitable extension
        # compute gene expression counts:
        sub.run(['htseq-count', '-f', extension, '-t', 'gene', '-q', sam_file,
                 annotation_file, '>', os.path.join(output_dir, sam_file[:-11])
                 + '.gene_counts'])
    else:
        raise ExtensionError("sam/bam file %s have not the correct extension. "
                             "Needed extension is '.sorted.[b|s]am'"%(sam_file))
    
    return 0

def create_matrix_of_counts(output_file, count_dir):
    """create matrix of counts from counts files (in count_dir) in output file.

    Output matrix of counts will be stored in output_file in csv format.
    
    Parameters
    ----------
    output_file: string
        output directory.
    count_dir: string
        directory where are stored gene expression counts as
        '$(cell_id).gene_counts'.
        
    Returns
    -------
    int : 0
    """

    # build matrix of counts file in csv format:
    with open(output_file, 'w') as matrix:
        writer = csv.writer(matrix, delimiter=',', lineterminator='\n')
        # create csv header:
        header = ['cells']
        head = os.listdir(count_dir)[0]
        if head[-12:] != ".gene_counts":
            raise ExtensionError("file %s doesn't have .gene_counts extension"
                                 %(head))
        with open(os.path.join(count_dir, head), 'r') as f:
            reader = csv.reader(f, delimiter='\t', lineterminator='\n')
            for row in reader:
                if row[0][:2] != "__": # exclude meta information
                    header.append(row[0]) # add gene id to header
        writer.writerow(header) # write the header
        for fichier in os.listdir(count_dir):
            if fichier[-12:] == ".gene_counts":
                with open(os.path.join(count_dir, fichier), 'r') as f:
                    reader = csv.reader(f, delimiter='\t', lineterminator='\n')
                    cell_row = [fichier[:-12]] # add cell id
                    for row in reader:
                        if row[0][:2] != "__": # exclude meta information
                            cell_row.append(row[1])
                    writer.writerow(cell_row) # write the row of current cell
    
    return 0

def extract_brie_psi_matrix(dict_of_cells, matrix_file):
    """extract psi from simple brie results and write results in a csv file.

    Psi is extracted from fractions.tsv.

    Parameters
    ---------
    dict_of_cells: dict
    dictionnary with cells id as keys and according brie output directories as
    values

    Returns
    -------
    numpy.array: matrix of predicted psi for each cell and genes
    """
    #build matrix of psis file in csv format:
    with open(matrix_file, 'w') as matrix:
        #writer = csv.DictWriter(matrix, delimiter=',', lineterminator='\n', fieldnames=fieldnames)
        #writer = csv.writer(matrix, delimiter=',', lineterminator='\n')
        
        header = ['cell'] # header for csv file
        for cell in dict_of_cells: # for each cell
            fractions_file = os.path.join(dict_of_cells[cell], 'fractions.tsv')
            
            if header == ['cell']: # if header had not be written yet
                with open(fractions_file, 'r') as f:
                    reader = csv.reader(f, delimiter='\t', lineterminator='\n')
                    ## create csv header
                    # every one over two transcripts:
                    for row in reader: # for each transcript
                        # if that transcript is exon-included:
                        if row[0][-3:] == ".in":
                            header += [row[0]] # add transcript id
                    #header += 'pseudotime'
                    writer = csv.DictWriter(matrix, delimiter=',',
                                        lineterminator='\n', fieldnames=header)
                    writer.writeheader()
                    #writer.writerow(header) # write header

            # extract psi from current cell
            with open(fractions_file, 'r') as f:
                reader = csv.reader(f, delimiter='\t', lineterminator='\n')
                d = { 'cell': cell } # dict that describes current cell row
                # every one over two transcripts:
                for row in reader:
                    if row[0][-3:] == ".in": # exon inclusion transcript
                        d[row[0]] = row[5]#d[transcript_id] = corresponding psi
                writer.writerow(d) # write the row of current cell
                #writer.writerow({'cells': cell, gene_row[0]: gene_row[5]})
                            
    return

def compute_correlation(matrix_file, pseudotimes):
    """compute for each gene Person's r between pseudotime and psi for each cell

    Returns results in a dict. WX = logit(psi).

    Parameters
    ----------
    matrix_file (string): name of the matrix of brie WX for each cell and gene.

    pseudotimes (panda.DataFrame): array where first column store cells'ids and
        second column stores corresponding pseudotimes (floats between 0 and 1).

    Returns
    -------
    dict : keys are transcript names and values are computed Pearson's r.
    """
    with open(matrix_file, 'r') as f:
        reader = csv.reader(f, delimiter=',', lineterminator='\n')
        header = reader.__next__() # get first row, ie header
        genes = header[1:] # get rid of 'cell' column attribute

    gene_corr_dict = {} # dict that store results
    for gene in genes: # for each gene
        WX = [] # predicted psi values for gene according to brie
        pseudotime = [] # pseudotimes given in the right order
        with open(matrix_file, 'r') as f:
            reader = csv.DictReader(f, delimiter=',', lineterminator='\n')
            for row in reader:
                WX.append(logistic(float(row[gene]))) # because row[gene] == psi
                pseudotime.append(pseudotimes[row['cell']])
                # p['dpt_pseudotime']['ERR1147410']
        correlation = st.pearsonr(pseudotime, WX)[0]
        gene_corr_dict[gene] = correlation
        
    return gene_corr_dict

def store_pseudotime(storage_file, pseudotimes):
    """store pseudotimes for each cell in a tsv file

    structure of file:
    cell_id    pseudotime_value

    Parameters
    ----------
    storage_file (string): file where to store pseudotime values for each cell.

    pseudotimes (panda.DataFrame): array where first column store cells'ids and
        second column stores corresponding pseudotimes (floats between 0 and 1).

    Returns
    -------
    """
    with open(storage_file, 'w') as f:
        header = ['cell', 'pseudotime']
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=header)
        writer.writeheader()
        # for each cell and corresponding pseudotime t
        for cell, t in pseudotimes.items():
            writer.writerow({"cell": cell, "pseudotime": t})
            # print("cell: ", cell)
            # print("t: ", t)
        
    return

def pseudotime(matrix_of_counts):
    """Compute pseudotimes of each cell of file of name matrix_of_counts

    Pseudotime analysis is computed through a diffusion map.

    Args:
    matrix_of_counts (string): name of the matrix of counts file

    Returns:
    pseudotimes (panda.DataFrame): array where first column store cells'ids and
        second column stores corresponding pseudotimes (floats between 0 and 1).
    """
    adata = sc.read(matrix_of_counts) # read matrix of counts
    adata.uns['iroot'] = 1 # choose an arbitrary root cell
    # compute neighbors to prepare for diffusion map:
    sc.pp.neighbors(adata, n_neighbors=5, method='gauss', knn=False, use_rep='X')
    # compute branching and diffusion pseudotime:
    sc.tl.dpt(adata, n_branchings=1)
    # extract pseudotimes (two first columns):
    pseudotimes = adata.obs['dpt_pseudotime']
    
    return pseudotimes


def main():
    """compute correlation between pseudotime and BRIE psi prediction

    Pearson's r correlation coefficient between WX (where W is the matrix of
    weights computed from BRIE and X given feature matrix) and pseudotime is
    computed for each gene.
    """

    args = parse_arguments()

    # extract genes and isoforms transcrits from annotation file
    anno_file = args.annotation_file
    genes = loadgene(anno_file) # load gene data from annotation file
    # gene_ids[i] store the gene id of transcript given by tran_ids[i]
    # hence, in 2 isoforms per gene scenario, genes[i].geneID==gene_ids[2*i]

    gene_ids, tran_ids = [], []
    for g in genes: # for each gene among the gene_nb first ones
        for t in g.trans:
            tran_ids.append(t.tranID)
            gene_ids.append(g.geneID)
    
    out_dir = args.output_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ## compute counts and/or if needed
    if args.matrix_of_counts != None: # if matrix of counts is provided
        matrix_of_counts = args.matrix_of_counts  
    else:
        if ((args.counts_dir == None) # if no counts' directory provided,
            or (not os.path.isdir(args.counts_dir)) # or it doesn't exist yet,
            or os.listdir(args.counts_dir)): # or it is empty:
            
            if (args.counts_dir == None): # if no counts' directory provided:
                counts_dir = os.path.join(args.output_dir, "counts")
            else:
                counts_dir = args.counts_dir
                
            if not os.path.isdir(counts_dir): # if directory do not exist
                os.makedirs(counts_dir) # create it
                
            # compute counts
            for sam_file in os.listdir(args.sam_dir): # for each file in sam_dir
                if sam_file[-11:] == ".sorted.bam": # if it is a sorted sam file
                    sam_file = os.path.join(args.sam_dir, sam_file)
                    # compute counts from this file
                    counts_from_sam(sam_file , counts_dir, args.annotation_file)

        # compute matrix of counts:
        matrix_of_counts = os.path.join(args.output_dir, "matrix_of_counts.csv")
        create_matrix_of_counts(matrix_of_counts, counts_dir)
        
    ## compute pseudotimes
    pseudotimes = pseudotime(matrix_of_counts)

    ## compute BRIE and get psi from fraction.tsv
    output_dir = os.path.join(args.output_dir, "brie_outputs")
    if not os.path.exists(output_dir): # if brie output_dir does not exist yet
        os.makedirs(output_dir) # create directory for brie outputs
    output_dict = {} # dict { cell ids: brie output directory }
    for sam_file in os.listdir(args.sam_dir): # for each file in sam_dir
        if (sam_file[-11:] == ".sorted.bam" # if it is a sorted bam file
            or sam_file[-11:] == ".sorted.sam"): # if it is a sorted sam file
            sam_file_location = os.path.join(args.sam_dir, sam_file)
            # create output directory for each single brie analysis
            cell_id = sam_file[:-11]
            output = os.path.join(output_dir, cell_id)

            # remove brie results directory if --overwrite-brie is provided:
            if os.path.exists(output) and args.overwrite_brie:
                sub.run(['rm', '-r', output])
                
            # create directory for brie outputs if needed
            if not os.path.exists(output):
                os.makedirs(output)
                # run brie analysis if output did not exist before
                sub.run(["brie",
                         "-o", output,
                         "-s", sam_file_location,
                         "-a", args.annotation_file,
                         "-f", args.factor_file]
                        + args.brie_arguments.split())
            # else, brie results will be taken from output directory
        
            output_dict[cell_id] = output

    ## compute correlation
    matrix_file = os.path.join(output_dir, 'WXmatrix.csv')
    extract_brie_psi_matrix(output_dict, matrix_file) # compute matrix of WX
    gene_corr_dict = compute_correlation(matrix_file, pseudotimes)

    # write the results in a file:
    correlation_file = os.path.join(out_dir, 'pseudotime-WX_correlation.tsv')
    with open(correlation_file, 'w') as f:
        for gene in gene_corr_dict: # for each exon inclusion transcript
            f.write(gene + '\t' + str(gene_corr_dict[gene]) + '\n')#write result

    # write pseudotime results in a file:
    pseudotime_file = os.path.join(out_dir,"pseudotimes.tsv")
    store_pseudotime(pseudotime_file, pseudotimes)


    ## run pseudotime brie analysis:
    pseudotime_auxiliary.main(["-o", output,
                               "-s", args.sam_dir,
                               "-a", args.annotation_file,
                               "-f", args.factor_file,
                               "--pseudotimes", pseudotime_file,
                               "--WX_matrix", matrix_file]
                              + args.brie_arguments.split())
    
    # if args.brie_pseudotime_path is not None:
        # sub.run(["python3.4", args.brie_pseudotime_path,
        #          "-o", output,
        #          "-s", args.sam_dir,
        #          "-a", args.annotation_file,
        #          "-f", args.factor_file,
        #          "--pseudotimes", pseudotime_file,
        #          "--WX_matrix", matrix_file]
        #         + args.brie_arguments.split())
        
        # os.system('longjob -28day -c "'
        #           + ' '.join(["python3.4", args.brie_pseudotime_path,
        #                       "-o", output,
        #                       "-s", args.sam_dir,
        #                       "-a", args.annotation_file,
        #                       "-f", args.factor_file,
        #                       "--pseudotimes", pseudotime_file,
        #                       "--WX_matrix", matrix_file]
        #                      + args.brie_arguments.split())
        #           + '"')

        # TODO: fix path (path must be accessible from this file)
        # >> inverse utils and main, fix path with args.brie_pseudotime_path
        # fix in script_pseudo
        
    return 0

if __name__ == "__main__":
    #store_pseudotime("/home/milan/prog/cours/internship/pseudotime/data/save/pseudotime_30.tsv", pseudotime("/home/milan/prog/cours/internship/pseudotime/data/save/matrix_of_counts_30.csv"))
    main()
