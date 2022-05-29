# This function is to count reads supporting each read category

import os
import sys
import time
import pysam
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

from ..version import __version__
from ..utils.count import get_smartseq_matrix, SE_effLen
from ..utils.count_droplet import get_droplet_matrix
from ..utils.io_utils import read_brieMM, read_gff, convert_to_annData


#TODO: something wrong as there is no reads for isoform 2. Please check!!

def smartseq_count(gff_file, samList_file, out_dir=None, nproc=1, 
        event_type='SE', verbose=False):
    """CLI for counting reads supporting isoforms
    """
    # Parameter check
    
    ## Load sam file list
    sam_table = np.loadtxt(samList_file, delimiter = None, dtype=str, ndmin = 2)
    print('[BRIE2] example head cells:')
    print(sam_table[:min(3, sam_table.shape[0])], "...")
    if sam_table.shape[1] == 1:
        sam_table = np.append(sam_table, 
            [['S%d' %x] for x in range(sam_table.shape[0])], axis=1)
        
    ## Check out_dir
    if out_dir is None:
        sam_dir = os.path.abspath(samList_file)
        out_dir = os.path.dirname(sam_dir) + "/brieCOUNT"
    else:
        out_dir = out_dir
        
    try:
        os.stat(os.path.abspath(out_dir))
    except:
        os.mkdir(os.path.abspath(out_dir))
        
    ## Load GFF file
    sys.stdout.write("\r[BRIE2] loading gene annotations ... ")
    sys.stdout.flush()
    genes = read_gff(gff_file)
    sys.stdout.write("\r[BRIE2] loading gene annotations ... Done.\n")
    sys.stdout.flush()
    
    # Running
    ## Output gene info
    gene_table = [["GeneID", "GeneName", "TranLens", "TranIDs"]]
    for g in genes:
        tran_ids, tran_lens = [], []
        for t in g.trans:
            tran_ids.append(t.tranID)
            tran_lens.append(str(t.tranL))
        out_list = [g.geneID, g.geneName, ",".join(tran_lens), 
                    ",".join(tran_ids)]
        gene_table.append(out_list)
        
    fid = open(out_dir + "/gene_note.tsv", "w")
    for out_list in gene_table:
        fid.writelines("\t".join(out_list) + "\n")
    fid.close()
        
    ## Output sam total counts
    reads_table = np.zeros(sam_table.shape[0])
    for i in range(sam_table.shape[0]):
        if not os.path.isfile(str(sam_table[i, 0])):
            print("Error: No such file\n    -- %s" %sam_table[i, 0])
            sys.exit(1)
        pysam_stats = pysam.idxstats(sam_table[i, 0])
        if type(pysam_stats) is not list:
            pysam_stats = pysam_stats.split("\n")
        for tp in pysam_stats: 
            tmp = tp.strip().split("\t")
            if len(tmp) >= 3:
                reads_table[i] += float(tmp[2])
                
    cell_table = [["samID", "samCOUNT"]]
    fid = open(out_dir + "/cell_note.tsv", "w")
    fid.writelines("samID\tsamCOUNT\n")
    for i in range(len(reads_table)):
        cell_table.append([sam_table[i, 1], reads_table[i]])
        fid.writelines("%s\t%d\n" %(sam_table[i, 1], reads_table[i]))
    fid.close()
    
    ## Generate isoform effective length matrix
    if event_type == 'SE':
        effLen_tensor = np.zeros((len(genes), 2, 3), dtype=np.float32)
        for ig, _gene in enumerate(genes):
            effLen_tensor[ig, :, :] = SE_effLen(_gene, rlen=76)
    else:
        # placeholder, not support yet
        effLen_tensor = np.zeros((len(genes), 1), dtype=np.float32)
    
    ## Load read counts
    print("[BRIE2] counting reads for %d genes in %d sam files with %d cores..." 
          %(len(genes), sam_table.shape[0], nproc))

    get_smartseq_matrix(genes, sam_table, out_dir, event_type="SE", 
        edge_hang=10, junc_hang=2, nproc=nproc, verbose=verbose)
    

    ## Don't save into h5ad if not SE
    if event_type != 'SE':
        sys.exit()
    
    ## Save data into h5ad
    sys.stdout.write("\r[BRIE2] saving matrix into h5ad ... ")
    sys.stdout.flush()

    Rmat_dict = read_brieMM(out_dir + "/read_count.mtx")
    adata = convert_to_annData(Rmat_dict=Rmat_dict, 
                               effLen_tensor=effLen_tensor, 
                               cell_note=np.array(cell_table, dtype='str'), 
                               gene_note=np.array(gene_table, dtype='str'))
    adata.uns['event_type'] = event_type
    adata.write_h5ad(out_dir + "/brie_count.h5ad")
    
    sys.stdout.write("\r[BRIE2] saving matrix into h5ad ... Done.\n")
    sys.stdout.flush()
    
    ## Save data into npz
    # np.savez(out_dir + "/brie_count.npz", 
    #         Rmat_dict=Rmat_dict, effLen_tensor=effLen_tensor, 
    #         cell_note=cell_table, gene_note=gene_table)


def droplet_count(gff_file, sam_file, barcode_file, out_dir=None, nproc=1, 
                  event_type='SE', CB_tag='CB', UMI_tag='UR', verbose=False):
    """CLI for counting reads supporting isoforms
    """
    ## TODO: Parameter check
    if sam_file == None:
        print("[BRIE2] Error: need --samFile for indexed & aliged sam/bam/cram file.")
        sys.exit(1)
    
    ## Load sam file list
    cell_list = np.loadtxt(barcode_file, delimiter = None, dtype=str, ndmin = 2)
    cell_list = cell_list[:, 0]
    print('[BRIE2] example head cells:')
    print(cell_list[:min(3, cell_list.shape[0])], "...")

    ## Check out_dir
    if out_dir is None:
        sam_dir = os.path.abspath(sam_file)
        out_dir = os.path.dirname(sam_dir) + "/brieCOUNT"
    else:
        out_dir = out_dir
        
    try:
        os.stat(os.path.abspath(out_dir))
    except:
        os.mkdir(os.path.abspath(out_dir))
        
    ## Load GFF file
    sys.stdout.write("\r[BRIE2] loading gene annotations ... ")
    sys.stdout.flush()
    genes = read_gff(gff_file)
    sys.stdout.write("\r[BRIE2] loading gene annotations ... Done.\n")
    sys.stdout.flush()
    
    ## Output gene info
    gene_table = [["GeneID", "GeneName", "TranLens", "TranIDs"]]
    for g in genes:
        tran_ids, tran_lens = [], []
        for t in g.trans:
            tran_ids.append(t.tranID)
            tran_lens.append(str(t.tranL))
        out_list = [g.geneID, g.geneName, ",".join(tran_lens), 
                    ",".join(tran_ids)]
        gene_table.append(out_list)
        
    fid = open(out_dir + "/gene_note.tsv", "w")
    for out_list in gene_table:
        fid.writelines("\t".join(out_list) + "\n")
    fid.close()
        
    ## Output sam total counts
    total_reads = 0
    if not os.path.isfile(sam_file):
        print("Error: No such file\n    -- %s" %sam_file)
        sys.exit(1)
    pysam_stats = pysam.idxstats(sam_file)
    if type(pysam_stats) is not list:
        pysam_stats = pysam_stats.split("\n")
    for tp in pysam_stats: 
        tmp = tp.strip().split("\t")
        if len(tmp) >= 3:
            total_reads += float(tmp[2])
                
    ### Output cell info
    fid = open(out_dir + "/cell_note.tsv", "w")
    fid.writelines("barcodes\n")
    for i in range(len(cell_list)):
        fid.writelines("%s\n" %(cell_list[i]))
    fid.close()
    
    ## Generate isoform effective length matrix
    if event_type == 'SE':
        effLen_tensor = np.zeros((len(genes), 2, 3), dtype=np.float32)
        for ig, _gene in enumerate(genes):
            effLen_tensor[ig, :, :] = SE_effLen(_gene, rlen=76)
    else:
        # placeholder, not support yet
        effLen_tensor = np.zeros((len(genes), 1), dtype=np.float32)
    
    ## Load read counts
    print("[BRIE2] counting reads for %d genes in %d cells with %d cores..." 
          %(len(genes), cell_list.shape[0], nproc))

    res = get_droplet_matrix(genes, sam_file, cell_list, out_dir, event_type, 
                             10, 2, CB_tag, UMI_tag, nproc, verbose)
    
    ## Don't save into h5ad if not SE
    if event_type != 'SE':
        sys.exit()
    
    ## Save data into h5ad
    sys.stdout.write("\r[BRIE2] saving matrix into h5ad ... ")
    sys.stdout.flush()

    Rmat_dict = read_brieMM(out_dir + "/read_count.mtx")
    cell_table = np.append(['barcodes'], cell_list).reshape(-1, 1)
    adata = convert_to_annData(Rmat_dict=Rmat_dict, 
                               effLen_tensor=effLen_tensor, 
                               cell_note=np.array(cell_table, dtype='str'), 
                               gene_note=np.array(gene_table, dtype='str'))
    adata.uns['event_type'] = event_type
    adata.uns['total_reads'] = total_reads
    adata.write_h5ad(out_dir + "/brie_count.h5ad")

    sys.stdout.write("\r[BRIE2] saving matrix into h5ad ... Done.\n")
    sys.stdout.flush()

    
def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--gffFile", "-a", dest="gff_file", default=None,
        help="GTF/GFF3 file for gene and transcript annotation")
    parser.add_option("--out_dir", "-o", dest="out_dir", default=None, 
        help="Full path of output directory [default: $samFile/brieCOUNT]")

    group0 = OptionGroup(parser, "SmartSeq-based input")
    group0.add_option("--samList", "-S", dest="samList_file", default=None,
        help=("A no-header tsv file listing sorted and indexed bam/sam/cram "
              "files. Columns: file path, cell id (optional)"))

    group1 = OptionGroup(parser, "Droplet-based input")
    group1.add_option("--samFile", "-s", dest="sam_file", default=None,
        help=("One indexed bam/sam/cram file"))
    group1.add_option("--barcodes", "-b", dest="barcodes_file", default=None,
        help=("A file containing cell barcodes without header"))
    group1.add_option("--cellTAG", dest="cell_tag", default="CB",
        help="Tag for cell barocdes [default: %default]")
    group1.add_option("--UMItag", dest="UMI_tag", default="UR",
        help="Tag for UMI barocdes [default: %default]")

    group2 = OptionGroup(parser, "Optional arguments")
    group2.add_option("--verbose", dest="verbose", default=False, 
        action="store_true", help="Print out detailed log info")    
    group2.add_option("--nproc", "-p", type="int", dest="nproc", default="4",
        help="Number of subprocesses [default: %default]")
    group2.add_option("--eventType", "-t", dest="event_type", default="SE",
        help="Type of splicing event for check. SE: skipping-exon; "
             "Any: no-checking [default: %default]")
    
    # group.add_option("--add_premRNA", action="store_true", dest="add_premRNA", 
    #     default=False, help="Add the pre-mRNA as a transcript")

    parser.add_option_group(group0)
    parser.add_option_group(group1)
    parser.add_option_group(group2)
    
    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to brie-count in BRIE v%s!\n" %(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)
        
    if options.gff_file == None:
        print("[BRIE2] Error: need --gffFile for gene annotation.")
        sys.exit(1)
    
    # bias mode is not supported yet
    add_premRNA = False
    
    if options.samList_file is not None:
        smartseq_count(options.gff_file, options.samList_file, options.out_dir, 
            options.nproc, options.event_type, options.verbose)
    else:
        droplet_count(options.gff_file, options.sam_file, options.barcodes_file, 
            options.out_dir, options.nproc, options.event_type, options.cell_tag,
            options.UMI_tag, options.verbose)


if __name__ == "__main__":
    main()
    
