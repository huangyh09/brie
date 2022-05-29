==============
brie-quant CLI
==============

The ``brie-quant`` CLI (in brie>=2.0.0) uses the newly developed variational 
inference methods scalable to large data sets, which works both in CPU or 
`GPU <install.html#gpu-usage>`_ with the TensorFlow Backend. 
For using BRIE1 (<=0.2.4) with an MCMC sampler, 
please refer to `BRIE1 <brie1.html>`_.

This command allows to quantify the splicing proportion Psi and detect
variable splicing events/genes along with cell-level features, e.g., cell type, 
disease condition, and development time. 
It can handle both alternative splicing isoforms or unspliced vs spliced ratios 
in a unified framework.

As a Bayesian method, the key philosophy of BRIE is to combine likelihood 
(data-driven) and prior (uninformative or informative). In BRIE2, a variety of 
prior settings are supported, as follows.

.. note::
   The mode2 (quant and diff) with cell features is the recommended option in 
   BRIE2, for both quantifying PSI and identifying differential momentum genes 
   for RNA velocity analysis and differential splicing events on alternative 
   splicing for either categorical or continuous covariates.

   For unspliced vs spliced in RNA velocity analysis, please add argument
   ``--layers spliced,unspliced``, as by default is set for alternative splicing.


Mode2-quant: Aggregated imputation
==================================

This mode requires argument ``--interceptMode gene``, which means fitting an 
offset for each gene (or splicing event). The offset denotes the mean of the 
prior distribution through aggregation, hence learns a prior 
shared by all cells on each gene. The benefit for this mode is that dimension 
reduction can be performed, e.g., PCA and UMAP on splicing PSI matrix. 
As there are many 
splicing events with low total reads, they will have a high variance in the 
estimation, and is often suggested filtered out, which will cause missing values.
Based on the cell aggregated imputation, most dimension reduction methods can be
used, even it doesn't support missing values.

Example command line for mode 2:

.. code-block:: bash

  brie-quant -i out_dir/brie_count.h5ad -o out_dir/brie_quant_aggr.h5ad --interceptMode gene
  
  
Mode2-diff: Variable splicing detection
=======================================

This mode requires argument ``-c`` for cell features and ``--LRTindex`` for the 
index (zero-based) of cell features to perform likelihood ratio test. Again we
suggest keeping the cell aggregation on each gene by ``--interceptMode gene``.

Then this mode will learn a prior from the given cell level features and perform
the second round fitting with leaving one feature out each time to calculate the 
ELBO gain, which can be used as a surrogate for Bayes factor.

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
first factor multiple sclerosis, but also want to consider the potential 
confounder in mouse strain, we could set the design matrix in ``cell_info.tsv`` 
file as follow, along with parameters ``--interceptMode gene --LRTindex 0``.

.. code-block:: bash

  samID   isEAE   isCD1
  SRR7102862      0       1
  SRR7103631      1       0
  SRR7104046      1       1
  SRR7105069      0       0


.. note::
   Be very careful on collinearity (i.e., redundancy) of your design matrix in 
   ``cell_info.tsv`` and the constant intercept (if use 
   ``--interceptMode gene``), which is similar to DEG in edgeR or DESeq2.

   By default, BRIE2 uses the ``--testBase=full`` mode to detect differential 
   splicing by comparing all features + constant intercept versus leaving the 
   testing feature(s) out. In this setting, if your features can't have 
   collinearity with intercept, e.g., you may need to remove one cell type out 
   if it covers all your cells.
   
   Alternatively, you can change to another strategy ``--testBase=null`` by 
   comparing the testing feature(s) + intercept versus intercept only, 
   like the example data on
   `Dentate Gyrus <brie2_dentateGyrus.html#BRIE2’s-differential-momentum-genes-(DMGs)>`_.



Mode 1: Imputation with gene features
=====================================

This mode is introduced in BRIE 1, where genomic sequences are leverage to 
learn a prior distribution of PSI. Then the predicted distribution of PSI from 
sequence features are combined with the likelihood obtained from the observed 
read counts. This framework provides a coherent way to combine the 
quantification from read counts and imputation from genomic features.

Initially, the inference in this mode is achieved by MCMC sample per cell
separated, which is generally slow for over hundreds of cells. BRIE2's 
variational framework keeps supporting the gene features and performs the 
inference for all cells in one go, though they are independent.

Example command line for Mode 1. We suggest use ``--interceptMode cell`` to 
learn an offset for each cell:

.. code-block:: bash

  brie-quant -i out_dir/brie_count.h5ad -o out_dir/brie_quant_gene.h5ad \
      -g $DATA_DIR/gene_feature.tsv --interceptMode cell


.. note::
   For the sake of convenience, we now recommend using Mode2-quant below to 
   perform imputation, which leverages the average PSI values in a cell 
   population to function as an informative prior.



Mode 0: None imputation
=======================

In this mode, the prior is an uninformative logit-normal distribution with mean=0, 
and learned variance. Therefore, if a splicing event in a gene doesn't have any
read, it will return a posterior with Psi's mean=0.5 and 95% confidence interval 
around 0.95 (most cases >0.9).

This setting is used if you have high covered data and you only want to 
calculate cells with sufficient reads for each interesting gene, e.g., by 
filtering out all genes with Psi_95CI > 0.3.

Otherwise, the 0.5 imputed genes will be confounded by the expression level, 
instead of the isoform proportion.

Example command line for mode 1:

.. code-block:: bash

  brie-quant -i out_dir/brie_count.h5ad -o out_dir/brie_quant_pure.h5ad --interceptMode None



All parameters
==============

There are more parameters for setting (``brie-quant -h`` always give the version 
you are using):

.. code-block:: html

    Usage: brie-quant [options]

    Options:
      -h, --help            show this help message and exit
      -i IN_FILE, --inFile=IN_FILE
                            Input read count matrices in AnnData h5ad or brie npz format.
      -c CELL_FILE, --cellFile=CELL_FILE
                            File for cell features in tsv[.gz] with cell and feature ids.
      -g GENE_FILE, --geneFile=GENE_FILE
                            File for gene features in tsv[.gz] with gene and feature ids.
      -o OUT_FILE, --out_file=OUT_FILE
                            Full path of output file for annData in h5ad [default:
                            $inFile/brie_quant.h5ad]
      --LRTindex=LRT_INDEX  Index (0-based) of cell features to test with LRT: All, None
                            or comma separated integers [default: None]
      --testBase=TEST_BASE  Features in testing base model: full, null  [default: full]
      --interceptMode=INTERCEPT_MODE
                            Intercept mode: gene, cell or None [default: None]
      --layers=LAYERS       Comma separated layers two or three for estimating Psi
                            [default: isoform1,isoform2,ambiguous]

      Gene filtering:
        --minCount=MIN_COUNT
                            Minimum total counts for fitltering genes [default: 50]
        --minUniqCount=MIN_UNIQ_COUNT
                            Minimum unique counts for fitltering genes [default: 10]
        --minCell=MIN_CELL  Minimum number of cells with unique count for fitltering genes
                            [default: 30]
        --minMIF=MIN_MIF    Minimum minor isoform frequency in unique count [default:
                            0.001]

      VI Optimization:
        --MCsize=MC_SIZE    Sample size for Monte Carlo Expectation [default: 3]
        --minIter=MIN_ITER  Minimum number of iterations [default: 5000]
        --maxIter=MAX_ITER  Maximum number of iterations [default: 20000]
        --batchSize=BATCH_SIZE
                            Element size per batch: n_gene * total cell [default: 500000]
        --pseudoCount=PSEUDO_COUNT
                            Pseudo count to add on unique count matrices [default: 0.01]
