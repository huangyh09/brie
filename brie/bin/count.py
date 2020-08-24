# This function is to count reads supporting each read category

import os
import sys
import time
import pysam
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

from brie.version import __version__
from brie.utils.io_utils import read_brieMM, read_gff, convert_to_annData
from brie.utils.count import get_count_matrix, SE_effLen


FID = None
PROCESSED = 0
TOTAL_BAMs = 0
START_TIME = time.time()

#TODO: something wrong as there is no reads for isoform 2. Please check!!
    
def show_progress(RV=None):    
    global PROCESSED, TOTAL_BAMs, START_TIME, FID
    
    if RV is None:
        return RV
    else:
        FID.writelines(RV)
    
    if RV is not None: 
        PROCESSED += 1
        bar_len = 20
        run_time = time.time() - START_TIME
        percents = 100.0 * PROCESSED / TOTAL_BAMs
        filled_len = int(bar_len * percents / 100)
        bar = '=' * filled_len + '-' * (bar_len - filled_len)
        
        sys.stdout.write('\r[BRIE2] [%s] %.1f%% cells done in %.1f sec.' 
            % (bar, percents, run_time))
        sys.stdout.flush()
        
    return RV



def count(gff_file, samList_file, out_dir=None, nproc=1, add_premRNA=False):
    """CLI for counting reads supporting isoforms
    """
    # Parameter check
    
    ## Load sam file list
    sam_table = np.genfromtxt(samList_file, delimiter = None, dtype=str)
    sam_table = sam_table.reshape(sam_table.shape[0], -1)
    print(sam_table[:min(3, sam_table.shape[0])], "...")
        
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
    global TOTAL_BAMs
    TOTAL_BAMs = sam_table.shape[0]
    
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
    effLen_tensor = np.zeros((len(genes), 2, 3), dtype=np.float32)
    for ig, _gene in enumerate(genes):
        effLen_tensor[ig, :, :] = SE_effLen(_gene, rlen=76)
    
    ## Load read counts
    print("[BRIE2] loading reads for %d genes in %d sam files with %d cores..." 
          %(len(genes), sam_table.shape[0], nproc))
    
    global START_TIME, FID
    START_TIME = time.time()
    
    FID = open(out_dir + "/read_count.mtx", "w")
    FID.writelines("%" + "%MatrixMarket matrix coordinate integer general\n")
    FID.writelines("%d\t%d\t%d\n" %(sam_table.shape[0], len(genes), 0))
    
    if nproc <= 1:
        for s in range(len(sam_table[:, 0])):
            res = get_count_matrix(genes, sam_table[s, 0], s, 10, 2)
            show_progress(res)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        result = []
        for s in range(len(sam_table[:, 0])):
            result.append(pool.apply_async(get_count_matrix, (genes, 
                sam_table[s, 0], s, 10, 2), callback=show_progress))
        pool.close()
        pool.join()
    
    FID.close()
    print("")
    
    ## Save data into h5ad
    print("[BRIE2] save matrix into h5ad ...")
    Rmat_dict = read_brieMM(out_dir + "/read_count.mtx")
    adata = convert_to_annData(Rmat_dict=Rmat_dict, 
                               effLen_tensor=effLen_tensor, 
                               cell_note=np.array(cell_table, dtype='str'), 
                               gene_note=np.array(gene_table, dtype='str'))
    adata.write_h5ad(out_dir + "/brie_count.h5ad")
    
    ## Save data into npz
    # np.savez(out_dir + "/brie_count.npz", 
    #         Rmat_dict=Rmat_dict, effLen_tensor=effLen_tensor, 
    #         cell_note=cell_table, gene_note=gene_table)
    
    
def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--gffFile", "-a", dest="gff_file", default=None,
        help="GTF/GFF3 file for gene and transcript annotation")
    parser.add_option("--samList", "-S", dest="samList_file", default=None,
        help=("A tsv file containing sorted and indexed bam/sam/cram files. "
              "No header line; file path and cell id (optional)"))
    parser.add_option("--out_dir", "-o", dest="out_dir", default=None, 
        help="Full path of output directory [default: $samList/brieCOUNT]")

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default="4",
        help="Number of subprocesses [default: %default]")
    # group.add_option("--add_premRNA", action="store_true", dest="add_premRNA", 
    #     default=False, help="Add the pre-mRNA as a transcript")

    parser.add_option_group(group)
    
    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to brie-count in BRIE v%s!\n" %(__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.samList_file == None:
        print("[BRIE2] Error: need --samList for indexed & aliged sam/bam/cram files.")
        sys.exit(1)
        
    if options.gff_file == None:
        print("[BRIE2] Error: need --gffFile for gene annotation.")
        sys.exit(1)
    
    # bias mode is not supported yet
    add_premRNA = False
    
    count(options.gff_file, options.samList_file, options.out_dir, 
          options.nproc, add_premRNA)


if __name__ == "__main__":
    main()
    
