==========
CLI Manual
==========

After properly installed BRIE (>=2.0.0) Python package, two CLI will be 
available directly in your Python path: ``brie-count``, ``brie-quant``. 

In this documentation, BRIE refers to BRIE2 (>=2.0.0). For using BRIE1 (<=0.2.2)
with MCMC sampler, please refer to BRIE1_. If you want to generate splicing
annotations for your data, e.g., a species different from human and mouse,
please use a separated package BRIE-kit_, which is developed in Python2. 

.. _BRIE1: https://brie.readthedocs.io/en/latest/brie1.html
.. _BRIE-kit: https://github.com/huangyh09/briekit/wiki


1. brie-count
=============

This CLI will return a count tensor for the number of reads aligned to four 
different categories in each splicing event and each cell:
1) Uniquely aligned to isoform1, e.g., in exon1-exon2 junction in SE event
2) Uniquely aligned to isoform2, e.g., in exon1-exon3 junction in SE event
3) Ambiguously aligned to isoform1 and isoform2, e.g., within exon1
0) Fetched in the region but not compatible with any of the above three types.
   We suggest ignore these reads.
   
As input, you need to generate the splicing event annotation. We have generated
data for human_ and mouse_. We suggest align RNA-seq reads to the according 
version of reference genome. Alternatively, you can use `briekit` package to 
generate.

.. _human: https://sourceforge.net/projects/brie-rna/files/annotation/human/gencode.v25/
.. _mouse: https://sourceforge.net/projects/brie-rna/files/annotation/mouse/gencode.vM12/


Then you fetch the counts on a list of bam files by the command line like this:

.. code-block:: bash

  brie-count -a AS_events/SE.gold.gtf -S sam_and_cellID.tsv -o out_dir -p 15

By default, you will have four output files in the out_dir: ``brie_count.h5ad``, 
``read_count.mtx.gz``, ``cell_note.tsv.gz``, and ``gene_note.tsv.gz``. The 
``brie_count.h5ad`` contains all information for downstream analysis, e.g., for
`brie-quant`.


There are more parameters for setting (``brie-count -h`` always give the version 
you are using):

.. code-block:: html

    Usage: brie-count [options]

    Options:
      -h, --help            show this help message and exit
      -a GFF_FILE, --gffFile=GFF_FILE
                            GTF/GFF3 file for gene and transcript annotation
      -S SAMLIST_FILE, --samList=SAMLIST_FILE
                            A tsv file containing sorted and indexed bam/sam/cram 
                            files. No header line; file path and cell id (optional)
      -o OUT_DIR, --out_dir=OUT_DIR
                            Full path of output directory [default: $samList/brieCOUNT]

      Optional arguments:
        -p NPROC, --nproc=NPROC
                            Number of subprocesses [default: 4]
  


2. brie-quant
=============


This command allows to quantify the splicing isoform proportion Psi and detect
variable splicing event along with cell level features, e.g., cell type, 
disease condition, development time.

As a Bayesian method, the key philosophy of BRIE is to combine likelihood (data 
driven) and prior (uninformative or informative). In BRIE2, a variety of prior
settings are supported, as follows.

2.1 Mode 1: None imputation
---------------------------

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


2.2 Mode 2: Aggregated imputation
---------------------------------

This mode requires argument `--interceptMode gene`. It aims to learn a prior 
shared by all cells on each gene. The benefit for this mode is that dimension 
reduction can be performed, e.g., PCA and UMAP on splicing. As there are many 
splicing events that are not well covered, it has a high variance in the 
estimation, and is often suggested filtered out, which will cause missing values.
Based on the cell aggregated imputation, most dimension reduction methods can be
used, even it doesn't support missing values.

Example command line for mode 2:

.. code-block:: bash

  brie-quant -i out_dir/brie_count.h5ad -o out_dir/brie_quant_aggr.h5ad --interceptMode gene
  
  
2.3 Mode 3: Variable splicing detection
---------------------------------------

This mode requires argument `-c` for cell features and `--LRTindex` for the 
index (zero-based) of cell features to perform likelihood ratio test. Again we
suggest to keep the cell aggregation on each gene by `--interceptMode gene`.

Then this mode will learn a prior from the given cell level features and perform
the second fit by leaving each feature out to calculate the EBLO gain, which 
can be further used as likelihood ratio test.

Example command line for mode 3:

.. code-block:: bash

  brie-quant -i out_dir/brie_count.h5ad -o out_dir/brie_quant_cell.h5ad \
      -c $DATA_DIR/cell_info.tsv --interceptMode gene --LRTindex=All


2.4 Flexible settings
---------------------

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
                            Input read count matrices in AnnData h5ad or brie npz format.
      -c CELL_FILE, --cellFile=CELL_FILE
                            File for cell features in tsv[.gz] with cell and feature ids.
      -g GENE_FILE, --geneFile=GENE_FILE
                            File for gene features in tsv[.gz] with gene and feature ids.
      -o OUT_FILE, --out_file=OUT_FILE
                            Full path of output file for annData in h5ad [default: $inFile/brie_quant.h5ad]
      --LRTindex=LRT_INDEX  Index (0-based) of cell features to test with LRT: 
                            All, None or comma separated integers [default: None]
      --interceptMode=INTERCEPT_MODE
                            Intercept mode: gene, cell or None [default: None]
      --layers=LAYERS       Comma separated layers two or three for estimating 
                            Psi [default: isoform1,isoform2,ambiguous]

      Optional arguments:
        --minCount=MIN_COUNT
                            Minimum total counts for fitltering genes [default: 50]
        --minUniqCount=MIN_UNIQ_COUNT
                            Minimum unique counts for fitltering genes [default: 10]
        --minCell=MIN_CELL  Minimum number of cells with unique count for fitltering genes [default: 30]
        --minIter=MIN_ITER  Minimum number of iterations [default: 5000]
        --maxIter=MAX_ITER  Maximum number of iterations [default: 20000]
        --batchSize=BATCH_SIZE
                            Element size per batch: n_gene * total cell [default: 500000]
