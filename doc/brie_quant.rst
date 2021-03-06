==============
brie-quant CLI
==============

The ``brie-quant`` CLI (in brie>=2.0.0) uses the newly developed variational 
inference methods scalable to large data sets, which works both in CPU or 
`GPU <install.html#gpu-usage>`_ with the TensorFlow Backend. 
For using BRIE1 (<=0.2.4) with MCMC sampler, 
please refer to `BRIE1 <brie1.html>`_.

This command allows to quantify the splicing isoform proportion Psi and detect
variable splicing event along with cell level features, e.g., cell type, 
disease condition, development time.

As a Bayesian method, the key philosophy of BRIE is to combine likelihood (data 
driven) and prior (uninformative or informative). In BRIE2, a variety of prior
settings are supported, as follows.

Mode 1: None imputation
=======================

In this mode, the prior is uninformative logit-normal distribution with mean=0, 
and learned variance. Therefore, if a splicing event in a gene doesn't have any
read, it will return a posterior with Psi's mean=0.5 and 95% confidence interval 
around 0.95 (most case >0.9).

This setting is used if you have high covered data and you only want to 
calculate cells with sufficient reads for each interesting genes, e.g., by 
filtering out all genes with Psi_95CI > 0.3.

Otherwise, the 0.5 imputed genes will be confounded by the expression level, 
instead of the isoform proportion.

Example command line for mode 1:

.. code-block:: bash

  brie-quant -i out_dir/brie_count.h5ad -o out_dir/brie_quant_pure.h5ad --interceptMode None


Mode 2: Aggregated imputation
=============================

This mode requires argument ``--interceptMode gene``. It aims to learn a prior 
shared by all cells on each gene. The benefit for this mode is that dimension 
reduction can be performed, e.g., PCA and UMAP on splicing. As there are many 
splicing events that are not well covered, it has a high variance in the 
estimation, and is often suggested filtered out, which will cause missing values.
Based on the cell aggregated imputation, most dimension reduction methods can be
used, even it doesn't support missing values.

Example command line for mode 2:

.. code-block:: bash

  brie-quant -i out_dir/brie_count.h5ad -o out_dir/brie_quant_aggr.h5ad --interceptMode gene
  
  
Mode 3: Variable splicing detection
===================================

This mode requires argument ``-c`` for cell features and ``--LRTindex`` for the 
index (zero-based) of cell features to perform likelihood ratio test. Again we
suggest to keep the cell aggregation on each gene by ``--interceptMode gene``.

Then this mode will learn a prior from the given cell level features and perform
the second fit by leaving each feature out to calculate the ELBO gain, which 
can be further used as likelihood ratio test.

Example command line for mode 3:

.. code-block:: bash

  brie-quant -i out_dir/brie_count.h5ad -o out_dir/brie_quant_cell.h5ad \
      -c $DATA_DIR/cell_info.tsv --interceptMode gene --LRTindex=All

**Example**

As an example in the 
`MS data <brie2_msEAE.html#BRIE2-option-1:-differential-splicing-events>`_, 
we have cells labelled with 
two factors 1) multiple sclerosis & control - isEAE, 2) two mouse strains 
- isCD1. We want to identify alternative splicing events associated with the 
first factor multiple sclerosis, but also want to considerthe potential 
confounder in mouse strain, we could set the design matrix in ``cell_info.tsv`` 
file as follow, along with parameters ``--interceptMode gene --LRTindex 0``.

.. code-block:: bash

  samID   isEAE   isCD1
  SRR7102862      0       1
  SRR7103631      1       0
  SRR7104046      1       1
  SRR7105069      0       0


.. note::
   Be very careful on collinearity (i.e., redudance) of your design matrix in 
   ``cell_info.tsv`` and the constant intercept (if use 
   ``--interceptMode gene``), which is similar as DEG in edgeR or DESeq2.

   By default, BRIE2 uses the ``--testBase=full`` mode to detect differential 
   splicing by comparing all features + constant intercept versus leaving the 
   testing feature(s) out. In this setting, if your features can't have 
   collinearity with intercept, e.g., you may need to remove one cell type out 
   if it covers all you cells.
   
   Alternatively, you can change to another strategy ``--testBase=null`` by 
   comparing the testing feature(s) + intercept versus intercept only, 
   like the example data on
   `Dentate Gyrus <brie2_dentateGyrus.html#BRIE2’s-differential-momentum-genes-(DMGs)>`_.


Flexible settings
=================

There could be more flexible settings, for example only use gene features as in
BRIE1 by the following command:

.. code-block:: bash

  brie-quant -i out_dir/brie_count.h5ad -o out_dir/brie_quant_gene.h5ad \
      -g $DATA_DIR/gene_seq_features.tsv --interceptMode cell --LRTindex=All
      
      
Or use both gene features and cell features
      
.. code-block:: bash

  brie-quant -i out_dir/brie_count.h5ad -o out_dir/brie_quant_all.h5ad \
      -c $DATA_DIR/cell_info.tsv -g $DATA_DIR/gene_seq_features.tsv \
      --interceptMode gene --LRTindex=All
      

There are more parameters for setting (``brie-quant -h`` always give the version 
you are using):

.. code-block:: html

    Usage: brie-quant [options]

    Options:
      -h, --help            show this help message and exit
      -i IN_FILE, --inFile=IN_FILE
                            Input read count matrices in AnnData h5ad or brie npz
                            format.
      -c CELL_FILE, --cellFile=CELL_FILE
                            File for cell features in tsv[.gz] with cell and
                            feature ids.
      -g GENE_FILE, --geneFile=GENE_FILE
                            File for gene features in tsv[.gz] with gene and
                            feature ids.
      -o OUT_FILE, --out_file=OUT_FILE
                            Full path of output file for annData in h5ad [default:
                            $inFile/brie_quant.h5ad]
      --LRTindex=LRT_INDEX  Index (0-based) of cell features to test with LRT:
                            All, None or comma separated integers [default: None]
      --interceptMode=INTERCEPT_MODE
                            Intercept mode: gene, cell or None [default: None]
      --layers=LAYERS       Comma separated layers two or three for estimating Psi
                            [default: isoform1,isoform2,ambiguous]

      Gene filtering:
        --minCount=MIN_COUNT
                            Minimum total counts for fitltering genes [default:
                            50]
        --minUniqCount=MIN_UNIQ_COUNT
                            Minimum unique counts for fitltering genes [default:
                            10]
        --minCell=MIN_CELL  Minimum number of cells with unique count for
                            fitltering genes [default: 30]
        --minMIF=MIN_MIF    Minimum minor isoform frequency in unique count
                            [default: 0.001]

      VI Optimization:
        --MCsize=MC_SIZE    Sample size for Monte Carlo Expectation [default: 3]
        --minIter=MIN_ITER  Minimum number of iterations [default: 5000]
        --maxIter=MAX_ITER  Maximum number of iterations [default: 20000]
        --batchSize=BATCH_SIZE
                            Element size per batch: n_gene * total cell [default:
                            500000]
