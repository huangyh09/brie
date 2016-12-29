BRIE: Bayesian Regression for Isoform Estimate
==============================================

About BRIE
----------

BRIE (Bayesian regression for isoform estimate) is a Bayesian method to estimate isoform proportions from RNA-seq data. Currently, BRIE could take sequence features to automatically learn informative prior of exon inclusion ratio in  exon-skippiing events. This informative prior is very important when limited data is available. In Bulk RNA-seq experiment, we could easily increase the amplification to get more sequencing reads to improve the accuracy of isoform estimate. However, in single cell RNA-seq (scRNA-seq) experiments, the initial molecular is very limited, which always results some genes with very low coverage or even drop-out. In scRNA-seq, the BRIE method, by integrating informative prior, e.g. learned from sequence feature, could provide accurate and reproducible estimates of splicing in single cells, as well as sensitive differential analyses.


BRIE provides following functions through command line:

1. ``brie``: Estimate isoform proportions and FPKM, and calculate weights for regulatory features.

2. ``brie-diff``: Calculate Bayes factor of differential splicing between two cells or two conditions. 

3. ``brie-event``: Extract the splicing events from gene annotation file in gff3/gtf formate. 

4. ``brie-factor``: Fetch genentic features from genome sequence reference file in fasta formate.


Quick Start
-----------

**Installation**: 

  - Download this repository
  - Type ``python setup.py install``; add ``--user`` if you don't have root permission and you don't use Anaconda_.

.. _Anaconda: https://www.continuum.io/anaconda-overview

**Arguments**

  - Type command line ``brie -h``




Detailed mannual
----------------

See the documentation_ on how to install, to use, to find the annotation data etc.

.. _documentation: http://brie-rna.sourceforge.net


References
----------

Coming soon