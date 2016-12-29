# This is a main file to run the diceseq software, which will return the 
# isoform proportions ratio for each gene at all time points.

import os
import sys
import time
import pysam
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

# import pyximport; pyximport.install()
from utils.gtf_utils import loadgene
from utils.run_utils import set_info, map_data, save_data
from models.model_brie import brie_MH_Heuristic

PROCESSED = 0
TOTAL_GENE = 0
TOTAL_READ = []
START_TIME = time.time()
    
def show_progress(RV=None):
    global PROCESSED, TOTAL_GENE, START_TIME
    if RV is not None: 
        PROCESSED += 1
        bar_len = 20
        run_time = time.time() - START_TIME
        percents = 100.0 * PROCESSED / TOTAL_GENE
        filled_len = int(bar_len * percents / 100)
        bar = '=' * filled_len + '-' * (bar_len - filled_len)
        
        sys.stdout.write('\r[Brie] [%s] %.1f%% done in %.1f sec.' 
            % (bar, percents, run_time))
        sys.stdout.flush()
    return RV


def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--anno_file", "-a", dest="anno_file", default=None,
        help="Annotation file for genes and transcripts in GTF or GFF3")
    parser.add_option("--sam_file", "-s", dest="sam_file", default=None,
        help=("Sorted and indexed bam/sam files, use ',' for replicates "
        "e.g., rep1.sorted.bam,sam1_rep2.sorted.bam"))
    parser.add_option("--out_dir", "-o", dest="out_dir", default="output", 
        help="Full path of output directory")
    parser.add_option("--factor_file", "-f", dest="factor_file", default=None,
        help=("HDF5 file with features to predict isoform expression."))

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default="4",
        help="Number of subprocesses [default: %default]")
    group.add_option("--weight_file", "-w", dest="weight_file", default=None,
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

    # group.add_option("--add_premRNA", action="store_true", dest="add_premRNA", 
    #     default=False, help="Add the pre-mRNA as a transcript")
    # parser.add_option("--feature_log",action="store_true",dest="feature_log",
    #     default=False, help="Use log scale for all features.")
    # parser.add_option("--two_isoform",action="store_true",dest="two_isoform",
    #     default=False, help=("Only two isoforms for all genes. This is mainly "
    #     "for splicing events."))

    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to Brie!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.anno_file == None:
        print("[Brie] Error: need --anno_file for annotation.")
        sys.exit(1)
    else:
        sys.stdout.write("\r[Brie] loading annotation file... ")
        sys.stdout.flush()
        # anno = load_annotation(options.anno_file, options.anno_source)
        genes = loadgene(options.anno_file)
        sys.stdout.write("\r[Brie] loading annotation file... Done.\n")
        sys.stdout.flush()
        # genes = anno["genes"]
        tran_len = []
        tran_ids = []
        gene_ids = []
        for g in genes:
            for t in g.trans:
                tran_len.append(t.tranL)
                tran_ids.append(t.tranID)
                gene_ids.append(g.geneID)
        gene_ids = np.array(gene_ids)
        tran_ids = np.array(tran_ids)
        tran_len = np.array(tran_len)

        global TOTAL_GENE, TOTAL_TRAN
        TOTAL_GENE = len(genes)
        TOTAL_TRAN = len(tran_ids)

    if options.sam_file == None:
        print("[Brie] Error: need --sam_file for indexed and aliged reads.")
        sys.exit(1)
    else:
        global TOTAL_READ
        TOTAL_READ = 0
        sam_file = options.sam_file
        for ss in sam_file.split(","):
            if not os.path.isfile(ss):
                print("Error: No such file\n    -- %s" %ss)
                sys.exit(1)
            pysam_stats = pysam.idxstats(ss)
            if type(pysam_stats) is not list:
                pysam_stats = pysam_stats.split("\n")
            for tp in pysam_stats: 
                tmp = tp.strip().split("\t")
                if len(tmp) >= 3:
                    TOTAL_READ += float(tmp[2])

    two_isoform = True
    if options.factor_file == None:
        feature_all = np.ones((len(tran_ids), 6))
        feature_ids = ["random%d" %i for i in range(1,6)] + ["intercept"]
        if two_isoform: 
            idxF = np.arange(0, len(tran_ids), 2)
            feature_all[idxF+1,:] = None
        else:
            idxF = np.arange(0, len(tran_ids))
        feature_all[idxF,:-1] = np.random.rand(len(idxF), 5)
    else:
        feature_all, feature_ids, idxF = map_data(options.factor_file,
            tran_ids, False) #options.feature_log

    if options.out_dir is None:
        out_dir = os.path.dirname(os.path.abspath(ss)) + "/brie_out"
    else:
        out_dir = options.out_dir
    try:
        os.stat(os.path.abspath(out_dir))
    except:
        os.mkdir(os.path.abspath(out_dir))

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

    no_twice = False
    auto_min = 200
    mate_mode = "pair"
    two_isoform = True
    add_premRNA = False
    print_detail = False

    nproc = options.nproc
    ftype = options.ftype
    FLmean, FLstd = options.frag_leng
    sample_num, M, initial, gap = options.mcmc_run


    print("[Brie] loading reads for %d genes with %d cores..." %(TOTAL_GENE, 
        nproc))
    global START_TIME
    START_TIME = time.time()

    R_all, len_iso_all, prob_iso_all = [], [], []
    if nproc <= 1:
        for g in genes:
            RV = set_info(g, sam_file, bias_mode, ref_file, bias_file, FLmean,
                FLstd, mate_mode, auto_min)
            show_progress(RV)
            R_all.append(RV["Rmat"])
            len_iso_all.append(RV["len_iso"])
            prob_iso_all.append(RV["prob_iso"])

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

    
    print("\n[Brie] running Brie for %d isoforms on %d genes with %d cores..." 
        %(TOTAL_TRAN, TOTAL_GENE, nproc))

    Psi_all, Y_all, RPK_all, Cnt_all, W_all, sigma_ = brie_MH_Heuristic(R_all,
        len_iso_all, prob_iso_all, feature_all, idxF, weights_in=weights_in, 
        _sigma=_sigma, _lambda=_lambda, ftype=ftype, total_count=TOTAL_READ, 
        M=M,  Mmin=initial, gap=gap, nproc=nproc)

    save_data(out_dir, sample_num, gene_ids, tran_ids, tran_len, feature_all,
        feature_ids, Psi_all, RPK_all, Cnt_all, W_all, sigma_)



if __name__ == "__main__":
    main()
    
