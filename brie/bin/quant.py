# This function is to quantify splcing isoforms and detect variable splicing 
# events associated with cell features

import os
import sys
import time
import numpy as np
from optparse import OptionParser, OptionGroup

import brie


def quant(in_file, cell_file=None, gene_file=None, out_file=None,
          LRT_index=[], layer_keys=['isoform1', 'isoform2', 'ambiguous'],
          intercept=None, intercept_mode='gene', nproc=1, min_counts=50, 
          min_counts_uniq=10, min_cells_uniq=30, min_MIF_uniq=0.001,
          min_iter=5000, max_iter=20000, MC_size=1, batch_size=500000,
          pseudo_count=0.01, base_mode='full'):
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
    if out_file is None:
        print("No given out_file, use the dir for input file.")
        out_file = os.path.dirname(os.path.abspath(in_file)) + "/brie_quant.h5ad"
    try:
        os.stat(os.path.dirname(os.path.abspath(out_file)))
    except:
        os.mkdir(os.path.dirname(os.path.abspath(out_file)))
    
    ## Load input data into anndata
    if in_file.endswith(".h5ad"):
        adata = brie.read_h5ad(in_file)
    if in_file.endswith(".npz"):
        adata = brie.read_npz(in_file)
    
    ## Match cell featutures
    if cell_file is not None:
        if cell_file.endswith('csv') or cell_file.endswith('csv.gz'):
            _delimeter = ","
        else:
            _delimeter = "\t"
        dat_tmp = np.genfromtxt(cell_file, dtype="str", delimiter=_delimeter)
        _idx = brie.match(adata.obs.index, dat_tmp[1:, 0]).astype(float)
        mm1 = _idx == _idx
        mm2 = _idx[mm1].astype(int)
        
        print("[BRIE2] %.1f%% cells are matched with features" 
              %(np.mean(mm1) * 100))

        Xc = dat_tmp[mm2+1, 1:].astype(np.float32)
        Xc_ids = dat_tmp[0, 1:]
        
        adata = adata[mm1, :]
    else:
        Xc = None
        Xc_ids = None
        
    ## Filter genes
    print("layers:", layer_keys)
    adata = brie.pp.filter_genes(adata, min_counts=min_counts,
                                 min_counts_uniq=min_counts_uniq, 
                                 min_cells_uniq=min_cells_uniq, 
                                 min_MIF_uniq=min_MIF_uniq,
                                 uniq_layers=layer_keys[:2],
                                 ambg_layers=layer_keys[2:], copy=True)
    
    ## Match gene features
    if gene_file is not None:
        if gene_file.endswith('csv') or gene_file.endswith('csv.gz'):
            _delimeter = ","
        else:
            _delimeter = "\t"
            
        dat_tmp = np.genfromtxt(gene_file, dtype="str", delimiter=_delimeter)
        _idx = brie.match(adata.var.index, dat_tmp[1:, 0]).astype(float)
        mm1 = _idx == _idx
        mm2 = _idx[mm1].astype(int)
        
        print("[BRIE2] %.1f%% genes are matched with features" 
              %(np.mean(mm1) * 100))

        Xg = dat_tmp[mm2+1, 1:].astype(np.float32)
        Xg_ids = dat_tmp[0, 1:]
        
        adata = adata[:, mm1]
    else:
        Xg = None
        Xg_ids = None
    
    print(adata)
    
    if 'unspliced' in adata.layers:
        tau_prior = [1, 1]
    else:
        tau_prior = [3, 27]
        
    ## Test genes with each cell features
    # model = brie.tl.fitBRIE(adata[:, :200])
    from brie.models import fitBRIE
    model = fitBRIE(adata, Xc=Xc, Xg=Xg, 
                    LRT_index=LRT_index, layer_keys=layer_keys, 
                    intercept=intercept, intercept_mode=intercept_mode,
                    min_iter=min_iter, max_iter=max_iter, 
                    MC_size=MC_size, batch_size=batch_size,
                    pseudo_count=pseudo_count, base_mode=base_mode,
                    tau_prior=tau_prior)
    
    adata.uns['brie_version'] = brie.__version__
    adata.uns['Xc_ids'] = Xc_ids
    adata.uns['Xg_ids'] = Xg_ids
    
    
    # Save adata
    adata.write_h5ad(out_file)
    
    # Save results tatble for identified splicing phenotypes
    out_table_file = ".".join(out_file.split('.')[:-1]) + '.brie_ident.tsv'
    df = brie.io.dump_results(adata)
    df.to_csv(out_table_file, sep='\t', header=True, 
              index=True, index_label='GeneID', float_format='%.3e')
    
    
def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--inFile", "-i", dest="in_file", default=None,
        help="Input read count matrices in AnnData h5ad or brie npz format.")
    parser.add_option("--cellFile", "-c", dest="cell_file", default=None,
        help=("File for cell features in tsv[.gz] with cell and feature ids."))
    parser.add_option("--geneFile", "-g", dest="gene_file", default=None, 
        help=("File for gene features in tsv[.gz] with gene and feature ids."))
    parser.add_option("--out_file", "-o", dest="out_file", default=None, 
        help="Full path of output file for annData in h5ad "
             "[default: $inFile/brie_quant.h5ad]")
    parser.add_option("--LRTindex", dest="LRT_index", default="None",
        help="Index (0-based) of cell features to test with LRT: All, None "
             "or comma separated integers [default: %default]")
    parser.add_option("--testBase", dest="test_base", default="full",
        help="Features in testing base model: full, null  [default: %default]")
    parser.add_option("--interceptMode", dest="intercept_mode", default="None",
        help="Intercept mode: gene, cell or None [default: %default]")
    parser.add_option("--layers", dest="layers", 
        default="isoform1,isoform2,ambiguous",
        help="Comma separated layers two or three for estimating Psi "
             "[default: %default]")
    
    group1 = OptionGroup(parser, "Gene filtering")
    group1.add_option("--minCount", type="int", dest="min_count", default=50,
        help="Minimum total counts for fitltering genes [default: %default]")
    group1.add_option("--minUniqCount", type="int", dest="min_uniq_count",
        default=10, help="Minimum unique counts for fitltering genes "
                         "[default: %default]")
    group1.add_option("--minCell", type="int", dest="min_cell", default=30,
        help="Minimum number of cells with unique count for fitltering genes "
             "[default: %default]")
    group1.add_option("--minMIF", type="float", dest="min_MIF", default=0.001,
        help="Minimum minor isoform frequency in unique count [default: %default]")
    
    group2 = OptionGroup(parser, "VI Optimization")
    group2.add_option("--MCsize", type="int", dest="MC_size", default=3,
        help="Sample size for Monte Carlo Expectation [default: %default]")
    group2.add_option("--minIter", type="int", dest="min_iter", default=5000,
        help="Minimum number of iterations [default: %default]")
    group2.add_option("--maxIter", type="int", dest="max_iter", default=20000,
        help="Maximum number of iterations [default: %default]")
    group2.add_option("--batchSize", type=int, dest="batch_size", default=500000, 
        help="Element size per batch: n_gene * total cell [default: %default]")
    group2.add_option("--pseudoCount", type=float, dest="pseudo_count", default=0.01, 
        help="Pseudo count to add on unique count matrices [default: %default]")
    # group.add_option("--nproc", "-p", type="int", dest="npoc", default="-1",
    #     help="Number of processes for computing [default: %default]")

    parser.add_option_group(group1)
    parser.add_option_group(group2)
    
    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to brie-quant in BRIE v%s!\n" %(brie.__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)
    if options.in_file == None:
        print("[BRIE2] Error: need --h5adFile for count matrices in annData.")
        sys.exit(1)
        
    if options.LRT_index.upper() == "NONE":
        LRT_index = []
    elif options.LRT_index.upper() == "ALL":
        LRT_index = None
    else:
        LRT_index = np.array(options.LRT_index.split(","), float).astype(int)
        
    intercept = None if options.intercept_mode.upper() in ["GENE", 'CELL'] else 0
    
    ## maximum number of threads (to fix)
    # if options.nproc != -1:
    #     tf.config.threading.set_inter_op_parallelism_threads(options.nproc)
    
    nproc = -1
            
    # run detection function
    quant(options.in_file, options.cell_file, options.gene_file, 
          options.out_file, LRT_index, options.layers.split(','),
          intercept, options.intercept_mode, nproc, options.min_count, 
          options.min_uniq_count, options.min_cell, options.min_MIF,
          options.min_iter, options.max_iter, options.MC_size, 
          options.batch_size, options.pseudo_count, options.test_base)

if __name__ == "__main__":
    main()
    
