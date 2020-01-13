# This function is to count reads supporting each read category

import os
import sys
import time
import pysam
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

from .utils.gtf_utils import loadgene
from .utils.run_utils import get_count_matrix

FID = None
PROCESSED = 0
TOTAL_GENE = 0
START_TIME = time.time()

#TODO: something wrong as there is no reads for isoform 2. Please check!!
    
def show_progress(RV=None):    
    global PROCESSED, TOTAL_GENE, START_TIME, FID
    
    if RV is None:
        return RV
    else:
        FID.writelines(RV)
    
    if RV is not None: 
        PROCESSED += 1
        bar_len = 20
        run_time = time.time() - START_TIME
        percents = 100.0 * PROCESSED / TOTAL_GENE
        filled_len = int(bar_len * percents / 100)
        bar = '=' * filled_len + '-' * (bar_len - filled_len)
        
        sys.stdout.write('\r[Brie] [%s] %.1f%% genes done in %.1f sec.' 
            % (bar, percents, run_time))
        sys.stdout.flush()
        
    return RV



def main():
    import warnings
    warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--gffFile", "-a", dest="gff_file", default=None,
        help="GTF/GFF3 file for gene and transcript annotation")
    parser.add_option("--samList", "-S", dest="samList_file", default=None,
        help=("A tsv file containing sorted and indexed bam/sam/cram files. "
              "No header line; file path and cell id (optional)"))
    parser.add_option("--out_dir", "-o", dest="out_dir", default=None, 
        help="Full path of output directory [default: $samList/brieOUT]")

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default="4",
        help="Number of subprocesses [default: %default]")
    group.add_option("--fLen", type="float", nargs=2, dest="frag_leng",
        default=[None,None], help=("Two arguments for fragment length: "
        "mean and standard diveation, default: auto-detected"))
    group.add_option("--bias", nargs=3, dest="bias_args",
        default=["unif","None","None"], help=("Three argments for bias "
        "correction: BIAS_MODE,REF_FILE,BIAS_FILE(s). BIAS_MODE: unif, end5, "
        "end3, both. REF_FILE: the genome reference file in fasta format. "
        "BIAS_FILE(s): bias files from dice-bias, use '---' for time specific "
        "files, [default: unif None None]"))
    group.add_option("--add_premRNA", action="store_true", dest="add_premRNA", 
        default=False, help="Add the pre-mRNA as a transcript")

    parser.add_option_group(group)
    
    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to Brie!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.samList_file == None:
        print("[Brie] Error: need --samList for indexed & aliged sam/bam/cram files.")
        sys.exit(1)
    else:
        sam_table = np.genfromtxt(options.samList_file, delimiter = None, dtype=str)
        sam_table = sam_table.reshape(sam_table.shape[0], -1)
        print(sam_table)
        
    if options.out_dir is None:
        sam_dir = os.path.abspath(samList_file)
        out_dir = os.path.dirname(sam_dir) + "/brieOUT"
    else:
        out_dir = options.out_dir
    try:
        os.stat(os.path.abspath(out_dir))
    except:
        os.mkdir(os.path.abspath(out_dir))
        
    if options.gff_file == None:
        print("[Brie] Error: need --gffFile for gene annotation.")
        sys.exit(1)
    else:
        sys.stdout.write("\r[Brie] loading gene annotations ... ")
        sys.stdout.flush()
        genes = loadgene(options.gff_file)
        sys.stdout.write("\r[Brie] loading gene annotations ... Done.\n")
        sys.stdout.flush()
    
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

    auto_min = 200
    mate_mode = "pair"
    add_premRNA = False

    nproc = options.nproc
    FLmean, FLstd = options.frag_leng
    
    ## Output gene info
    fid = open(out_dir + "/gene_note.tsv", "w")
    fid.writelines("GeneID\tGeneName\tTranLens\tTranIDs\n")
    for g in genes:
        tran_ids, tran_lens = [], []
        for t in g.trans:
            tran_ids.append(t.tranID)
            tran_lens.append(str(t.tranL))
        out_list = [g.geneID, g.geneName, ",".join(tran_lens), 
                    ",".join(tran_ids)]
        fid.writelines("\t".join(out_list) + "\n")
    fid.close()
    global TOTAL_GENE
    TOTAL_GENE = len(genes)
        
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
                
    fid = open(out_dir + "/cell_note.tsv", "w")
    fid.writelines("samID\tsamCOUNT\n")
    for i in range(len(reads_table)):
        fid.writelines("%s\t%d\n" %(sam_table[i, 1], reads_table[i]))
    fid.close()
    
    ## Load read counts
    print("[Brie] loading reads for %d genes in %d sam files with %d cores..." 
          %(TOTAL_GENE, sam_table.shape[0], nproc))
    
    global START_TIME, FID
    START_TIME = time.time()
    
    FID = open(options.out_dir + "/read_count.mtx", "w")
    
    if nproc <= 1:
        for g in range(len(genes)):
            RV = get_count_matrix(genes[g], g, gsam_table[:, 0], bias_mode, 
                                  ref_file, bias_file, FLmean, FLstd, mate_mode, 
                                  auto_min)
            show_progress(RV)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        result = []
        for g in range(len(genes)):
            result.append(pool.apply_async(get_count_matrix, (genes[g], g, 
                sam_table[:, 0], bias_mode, ref_file, bias_file, FLmean, 
                FLstd, mate_mode, auto_min), callback=show_progress))
        pool.close()
        pool.join()
    
    FID.close()
    print("")


if __name__ == "__main__":
    main()
    
