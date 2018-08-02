"""compute pseudotime enhanced brie

   pseudotime and individual cell outputs of brie analysis are assumed to be
   stored in precoputed_brie_dir, supposed to contain a directory for each cell
   with cell id as a name, which each contain fractions.tsv output.
   Peudotime and brie output (WX) for each cell are stored in.
"""

import os
import sys
import time
import pysam # Working with BAM/CRAM/SAM-formatted files
import numpy as np
import multiprocessing
import csv # working with csv files
from optparse import OptionParser, OptionGroup

from utils.gtf_utils import loadgene # get list of Gene objects from gff/gtf
from utils.run_utils import set_info, map_data, save_data
from models.pseudotime_model_brie import brie_MH_Heuristic

# brie imports:
from brie import show_progress

PROCESSED = 0
TOTAL_GENE = 0
TOTAL_READ = []
START_TIME = time.time()

def parse_arguments(arguments=None):
    # parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="Annotation file for genes and transcripts in GTF or GFF3")
    parser.add_option("--sam_dir", "-s", dest="sam_dir", default=None,
        help=("Directory where are stored sorted and indexed bam/sam files."))
    parser.add_option("--out_dir", "-o", dest="out_dir", default="output", 
        help="Full path of output directory")
    parser.add_option("--factor_file", "-f", dest="factor_file", default=None,
        help=("Features in csv.gz file to predict isoform expression."))
    parser.add_option("--pseudotimes", default=None,
        help=("path to file with pseudotime for each cell in tsv format."))
    parser.add_option("--WX_matrix", default=None,
        help=("path to csv file with output WX from brie analysis."))

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default="4",
        help="Number of subprocesses [default: %default]")
    group.add_option("--weight_file", "-w", dest="weight_file", default=None, #! weights_dir
        help=("File with weights, an output of Brie."))
    group.add_option("--ftype", "-y", dest="ftype", default="Y",
        help="Type of function target: FPKM, Y, Psi [default: %default].")
    
    group.add_option("--fLen", type="float", nargs=2, dest="frag_leng",
        default=[None,None], help=("Two arguments for fragment length: "
        "mean and standard diveation, default: auto-detected"))
    group.add_option("--bias", nargs=3, dest="bias_args",
        default=["unif","None","None"], help=("Three argments for bias "
        "correction: BIAS_MODE,REF_FILE,BIAS_FILE(s). BIAS_MODE: unif, end5, "
        "end3, both. REF_FILE: the genome reference file in fasta format. "
        "BIAS_FILE(s): bias files from dice-bias, use '---' for time specific "
        "files, [default: unif None None]"))

    group.add_option("--sigma", dest="_sigma", default=None,
        help=("Sigma in Bayesian regression: the Gaussian standard deviation "
        "of residues [default: Auto]."))
    group.add_option("--lambda", dest="_lambda", default="0.1",
        help=("Lambda in Bayesian regression: the coeffiecient of "
        "L2 constrain on weights [default: %default]."))

    group.add_option("--mcmc", type="int", nargs=4, dest="mcmc_run",
        default=[500,5000,1000,50], help=("Four arguments for in MCMC "
        "iterations: save_sample,max_run,min_run,gap_run. Required: "
        "save_sample =< 3/4*mim_run. [default: 500 5000 1000 50]"))

    parser.add_option_group(group)

    if arguments is None:
        return parser.parse_args()
    else:
        return parser.parse_args(arguments)

def extract_info(WX_matrix_file, pseudotime_file):
    
    cell = {}
    return cell

def main(arguments=None):
    import warnings
    warnings.filterwarnings('error') # turn warnings into exceptions

    # if arguments is None:
    #     (options, args) = parse_arguments() # parse script's arguments
    # else: # if arguments are provided
    #     (options, args) = parse_arguments(arguments) # parse script's arguments
    (options, args) = parse_arguments(arguments) # parse script's arguments

    ###for each cell, run brie###

    ###appel de function generate data from brie_output_directory###

    if len(sys.argv[1:]) == 0:
        print("Welcome to pseudotime Brie!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.anno_file == None:
        print("[Brie] Error: need --anno_file for annotation.")
        sys.exit(1)
    else:
        sys.stdout.write("\r[Brie] loading annotation file... ")
        sys.stdout.flush()
        # list of genes (see Gene in .utils.gtf_utils):
        genes = loadgene(options.anno_file)
        sys.stdout.write("\r[Brie] loading annotation file... Done.\n")
        sys.stdout.flush()
        tran_len = [] # lenths of transcripts
        tran_ids = [] # ids of transcripts
        gene_ids = [] # ids of the gene of each transcript
        for g in genes: # for each gene, fill transcripts data
            for t in g.trans: # for each transcript (isoform)
                tran_len.append(t.tranL)
                tran_ids.append(t.tranID)
                gene_ids.append(g.geneID)
        # convert transcripts data as numpy.arrays:
        gene_ids = np.array(gene_ids)
        tran_ids = np.array(tran_ids)
        tran_len = np.array(tran_len)

        global TOTAL_GENE, TOTAL_TRAN # to modify TOTAL_GENE, TOTAL_TRAN
        TOTAL_GENE = len(genes) # number of genes
        TOTAL_TRAN = len(tran_ids) # number of transcripts

    # magical constants:
    two_isoform = True # currently, BRIE supports only the two isoforms scenario
    no_twice = False
    auto_min = 200
    mate_mode = "pair"
    add_premRNA = False
    print_detail = False

    nproc = options.nproc # number of subprocesses
    ftype = options.ftype
    FLmean, FLstd = options.frag_leng
    sample_num, M, initial, gap = options.mcmc_run

    # define features:
    if options.factor_file == None: # if there is no feature file (.csv.gz)
        # create random features (BRIE.NULL):
        feature_all = np.ones((len(tran_ids), 6))
        # range(1,6) to generate five uniformly distributed random features:
        feature_ids = ["random%d" %i for i in range(1,6)] + ["intercept"]
        if two_isoform: 
            idxF = np.arange(0, len(tran_ids), 2)#list of int from 0 to len(.) by 2
            feature_all[idxF+1,:] = None
        else:
            idxF = np.arange(0, len(tran_ids))
        feature_all[idxF,:-1] = np.random.rand(len(idxF), 5)
    else: # if there is a feature file:
        feature_all, feature_ids, idxF = map_data(options.factor_file,
            tran_ids, False) #options.feature_log

    # define weights options
    if options.weight_file == None:
        weights_in = None
    else:
        weights_in = np.loadtxt(options.weight_file, dtype="str", skiprows=1)
        weights_in = weights_in[:, 1].astype(float).reshape(-1)
    if options._sigma is None:
        _sigma = None
    else:
        _sigma = float(options._sigma)
    if options._lambda is None:
        _lambda = None
    else:
        _lambda = float(options._lambda)

    # define bias mode
    bias_mode, ref_file, bias_file = options.bias_args
    if bias_mode == "unif":
        ref_file = None
        bias_file = None
    elif ref_file is "None": 
        ref_file = None
        bias_file = None
        bias_mode = "unif"
        print("[Brie] No reference sequence, change to uniform mode.")
    elif bias_file is "None":
        ref_file = None
        bias_file = None
        bias_mode = "unif"
        print("[Brie] No bias parameter file, change to uniform mode.")
    else:
        bias_file = bias_file.split("---")

    # define sam_dir
    if options.sam_dir == None:
        print("[Brie] Error: need --sam_dir for indexed and aligned reads.")
        sys.exit(1)
    else:
        sam_dir = options.sam_dir

    # this dictionary will store every required information for each cell:
    cell_dict = {} # cell_dict[cell_id] represents cell of id cell_id
    ### for each cell
    for f in os.listdir(sam_dir): # for each sam_file in sam_dir:
        if (f[-11:] == ".sorted.bam" # if it is a sorted bam file
            or f[-11:] == ".sorted.sam"): # if it is a sorted sam file
            sam_file = os.path.join(sam_dir, f)
            # create output directory for each single brie analysis
            cell_id = f[:-11] # cell id, without file extention
            cell_dict[cell_id] = {}
        
            # count number of reads:
            global TOTAL_READ
            TOTAL_READ = 0
            if not os.path.isfile(sam_file):
                print("Error: No such file\n    -- %s" %sam_file)
                sys.exit(1)
            pysam_stats = pysam.idxstats(sam_file)
            if type(pysam_stats) is not list:
                pysam_stats = pysam_stats.split("\n")
            for tp in pysam_stats: 
                tmp = tp.strip().split("\t")
                if len(tmp) >= 3:
                    TOTAL_READ += float(tmp[2])
            # store TOTAL_READ in cell_dict:
            cell_dict[cell_id]["total_count"] = TOTAL_READ

            if options.out_dir is None: # if no output_directory is specified
                # create default output directory:
                out_dir = os.path.dirname(os.path.abspath(sam_file)) + "/pseudotime_brie_out/" + cell_id
            else:
                out_dir = os.path.join(options.out_dir, cell_id)
            try: # create output_dir if needed
                os.stat(os.path.abspath(out_dir))
            except:
                os.mkdir(os.path.abspath(out_dir))
            # store out_dir name in cell_dict:
            cell_dict[cell_id]["out_dir"] = out_dir
    
            # load reads
            print("[Brie] loading reads for %d genes with %d cores..." %(TOTAL_GENE, 
                nproc))
            global START_TIME
            START_TIME = time.time() # initialize time for BRIE analysis computing time

            R_all, len_iso_all, prob_iso_all = [], [], []
            if nproc <= 1: # in case of only one process
                for g in genes:
                    RV = set_info(g, sam_file, bias_mode, ref_file, bias_file, FLmean,
                        FLstd, mate_mode, auto_min)
                    show_progress(RV, total_gene=TOTAL_GENE,
                                  cell_number=len(cell_dict))
                    R_all.append(RV["Rmat"]) # R_mat: 2-D array_like, (N, K)
                    # N reads identities of belonging to K isoforms.
                    len_iso_all.append(RV["len_iso"])# prob_isos: 2-D array_like, (N, K)
                    # N reads probablity of belonging to K isoforms.
                    prob_iso_all.append(RV["prob_iso"]) # len_isos: 1-D array_like, (K,)
                    # The effective length of each isoform.

            else:
                pool = multiprocessing.Pool(processes=nproc)
                result = []
                for g in genes:
                    result.append(pool.apply_async(set_info, (g, sam_file, bias_mode,
                        ref_file, bias_file, FLmean, FLstd, mate_mode, auto_min), 
                        callback=show_progress))
                pool.close()
                pool.join()
                for res in result:
                    RV = res.get()
                    R_all.append(RV["Rmat"])
                    len_iso_all.append(RV["len_iso"])
                    prob_iso_all.append(RV["prob_iso"])

            # store relevant data in cell_dict
            cell_dict[cell_id]["R_mat"] = R_all
            cell_dict[cell_id]["len_isos"] = len_iso_all
            cell_dict[cell_id]["prob_isos"] = prob_iso_all

    ## exctract pseudotime and WX info from options.pseudotimes and WX_matrix
    # WX
    with open(options.WX_matrix, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            if row[0] in cell_dict: # if cell_id is in cell_dict (avoid header)
                cell_dict[row[0]]['WX'] = row[1:] # add corresponding WX
    # pseudotime:
    with open(options.pseudotimes, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0] in cell_dict: # if cell_id is in cell_dict (avoid header)
                cell_dict[row[0]]['t'] = row[1] # add corresponding WX
        

    print("\n[Brie] running Brie for %d isoforms on %d genes with %d cores..." 
          %(TOTAL_TRAN, TOTAL_GENE, nproc))
    
    results, W_all, sigma_ = brie_MH_Heuristic(cell_dict, feature_all, idxF,
                  weights_in=weights_in, _sigma=_sigma, _lambda=_lambda,
                  ftype=ftype, M=M, Mmin=initial, gap=gap, nproc=nproc)
    #brie_MH_Heuristic gives back an object like result["cell_id"]["Psi_all"]

    
    for _id in results: #for cell_id in results:
        # extract relevant variables
        out_dir = cell_dict[cell_id]["out_dir"]
        Psi_all = cell_dict[cell_id]["Psi_all"]
        Y_all = cell_dict[cell_id]["Y_all"]
        RPK_all = cell_dict[cell_id]["FPKM_all"]
        Cnt_all = cell_dict[cell_id]["Cnt_all"]
        ###out_dir = os.path.join(, "pseudotimes_fractions.tsv")
        save_data(out_dir, sample_num, gene_ids, tran_ids, tran_len, feature_all,
                  feature_ids, Psi_all, RPK_all, Cnt_all, W_all, sigma_) # checked



if __name__ == "__main__":
    main()
    
# add missing data from brie in cell !!
# check how will be stored output data
# read again this file and launch it
