# This function is to count reads supporting each read category

import os
import sys
import time
import pysam
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

from .version import __version__
from .utils.gtf_utils import load_genes
from .utils.count import get_count_matrix, SE_probability

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
        
        sys.stdout.write('\r[BRIE2] [%s] %.1f%% genes done in %.1f sec.' 
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
    # group.add_option("--add_premRNA", action="store_true", dest="add_premRNA", 
    #     default=False, help="Add the pre-mRNA as a transcript")

    parser.add_option_group(group)
    
    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to BRIE2 v%s!\n" %(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.samList_file == None:
        print("[BRIE2] Error: need --samList for indexed & aliged sam/bam/cram files.")
        sys.exit(1)
    else:
        sam_table = np.genfromtxt(options.samList_file, delimiter = None, dtype=str)
        sam_table = sam_table.reshape(sam_table.shape[0], -1)
        print(sam_table)
        
    if options.out_dir is None:
        sam_dir = os.path.abspath(options.samList_file)
        out_dir = os.path.dirname(sam_dir) + "/brieOUT"
    else:
        out_dir = options.out_dir
    try:
        os.stat(os.path.abspath(out_dir))
    except:
        os.mkdir(os.path.abspath(out_dir))
        
    if options.gff_file == None:
        print("[BRIE2] Error: need --gffFile for gene annotation.")
        sys.exit(1)
    else:
        sys.stdout.write("\r[BRIE2] loading gene annotations ... ")
        sys.stdout.flush()
        genes = load_genes(options.gff_file)
        sys.stdout.write("\r[BRIE2] loading gene annotations ... Done.\n")
        sys.stdout.flush()
    
    # bias mode is not supported yet
    add_premRNA = False
    nproc = options.nproc

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
    print("[BRIE2] loading reads for %d genes in %d sam files with %d cores..." 
          %(TOTAL_GENE, sam_table.shape[0], nproc))
    
    global START_TIME, FID
    START_TIME = time.time()
    
    FID = open(options.out_dir + "/read_count.mtx", "w")
    FID.writelines("%" + "%MatrixMarket matrix coordinate integer general\n")
    FID.writelines("%d\t%d\t%d\n" %(TOTAL_GENE, sam_table.shape[0], 0))
    
    pool = multiprocessing.Pool(processes=nproc)
    result = []
    for g in range(len(genes)):
        result.append(pool.apply_async(get_count_matrix, (genes[g], g, 
            sam_table[:, 0], 10, 2), callback=show_progress))
    pool.close()
    pool.join()
    
    FID.close()
    print("")


if __name__ == "__main__":
    main()
    
