# This function is to quantify splcing isoforms and detect variable splicing 
# events associated with cell features

import os
import sys
import time
import pysam
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

import brie

def quant(h5ad_file, cell_file, gene_file=None, out_dir=None, nproc=1, 
          min_counts=50, min_counts_uniq=10, min_iter=5000, max_iter=20000):
    """CLI for quantifying splicing isoforms and detecting variable splicing 
    events associated with cell features
    
    Parameters
    ----------
    cell_file: str 
        file path for cell features with cell label and feature names. Can be 
        compressed with gzip
    gene_file: str
        file path for gene features with gene label and feature names. Can be 
        compressed with gzip
    """
    # Parameter check   
    if out_dir is None:
        print("No given out_dir, use the dir for h5ad file.")
        out_dir = os.path.dirname(os.path.abspath(h5ad_file)) + "/brieQuant"
    else:
        out_dir = out_dir
        
    try:
        os.stat(os.path.abspath(out_dir))
    except:
        os.mkdir(os.path.abspath(out_dir))
    
    ## Load anndata
    adata = brie.read_h5ad(h5ad_file)
    print(adata)
    
    ## Filter genes
    brie.pp.filter_genes(adata, min_counts=min_counts, 
                         min_counts_uniq=min_counts_uniq)
    
    ## Match cell featutures
    if cell_file is not None:
        dat_tmp = np.genfromtxt(gene_file, dtype="str", delimiter="\t")
        _idx = brie.match(adata.obs.index, dat_tmp[1:, 0])
        print(np.mean(adata.obs.index == dat_tmp[_idx+1, 0]))

        Xc = dat_tmp[_idx+1, :].astype(np.float32).transpose()
        Xc_ids = dat_tmp[0, 1:]
    else:
        Xc = None
    
    ## Match gene features
    if gene_file is not None:
        dat_tmp = np.genfromtxt(gene_file, dtype="str", delimiter=",")
        _idx = brie.match(adata.var.index, dat_tmp[1:, 0])
        print(np.mean(adata.var.index == dat_tmp[_idx+1, 0]))

        Xg = dat_tmp[_idx+1, 1:].astype(np.float32)
        Xg_ids = dat_tmp[0, 1:]
    else:
        Xg = None
    
    ## Test genes with each features
    model = brie.tl.fitBRIE(adata, Xc=Xc, Xg=Xg, LRT_index=None, 
                            min_iter=min_iter, max_iter=max_iter)
        
    adata.write_h5ad(out_dir + "/brie_quant.h5ad")
    
#     fid = open(out_dir + "/results_brie_detect.tsv", "w")
#     for out_list in gene_table:
#         fid.writelines("\t".join(out_list) + "\n")
#     fid.close()
        
    
def main():
    import warnings
    warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--h5adFile", "-a", dest="h5ad_file", default=None,
        help="AnnData of read counts in h5ad format.")
    parser.add_option("--cellFile", "-c", dest="cell_file", default=None,
        help=("File for cell features in tsv[.gz] with cell and feature ids."))
    parser.add_option("--geneFile", "-g", dest="gene_file", default=None, 
        help=("File for gene features in tsv[.gz] with gene and feature ids."))
    parser.add_option("--out_dir", "-o", dest="out_dir", default=None, 
        help="Full path of output directory [default: $h5adFile/brieDetect]")

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default="4",
        help="Number of subprocesses [default: %default]")
    group.add_option("--minCount", type="int", dest="min_count", default=50,
        help="Minimum total counts for fitltering genes [default: %default]")
    group.add_option("--minUniqCount", type="int", dest="min_uniq_count", default=10,
        help="Minimum unique counts for fitltering genes [default: %default]")
    group.add_option("--minIter", type="int", dest="min_iter", default=5000,
        help="Minimum number of iterations [default: %default]")
    group.add_option("--maxIter", type="int", dest="max_iter", default=20000,
        help="Maximum number of iterations [default: %default]")


    parser.add_option_group(group)
    
    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to brie-quant in BRIE v%s!\n" %(brie.__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.h5ad_file == None:
        print("[BRIE2] Error: need --h5adFile for count matrices in annData.")
        sys.exit(1)
        
    # run detection function
    quant(options.h5ad_file, options.cell_file, options.gene_file, 
          options.out_dir, options.nproc, options.min_count, 
          options.min_uniq_count, options.min_iter, options.max_iter)

if __name__ == "__main__":
    main()
    

