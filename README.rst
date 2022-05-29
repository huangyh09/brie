|PyPI| |Docs| |Build Status|

.. |PyPI| image:: https://img.shields.io/pypi/v/brie.svg
    :target: https://pypi.org/project/brie
.. |Docs| image:: https://readthedocs.org/projects/brie/badge/?version=latest
   :target: https://brie.readthedocs.io
.. |Build Status| image:: https://travis-ci.org/huangyh09/brie.svg?branch=master
   :target: https://travis-ci.org/huangyh09/brie


BRIE: Bayesian Regression for Isoform Estimate
==============================================

Top News
--------
* [29/05/2022] We have released v2.2 that fully supports counting droplet-based 
  data for both Skipping Exon events and other types of splcing events. See the
  `brie-count manual <https://brie.readthedocs.io/en/latest/brie_count.html>`_

* [29/05/2022] We have include small-sized test data sets (15MB) for both 
  smart-seq2 and 10x Genomics. See data in `brie-tutorials/tests repo 
  <https://github.com/huangyh09/brie-tutorials/tree/main/tests>`_


About BRIE
----------

Welcome to the new BRIE (>=2.0 or BRIE2), Bayesian Regression for Isoform 
Estimate, a scalable Bayesian method to accurately identify splicing phenotypes 
in single-cell RNA-seq experiments and quantify isoform proportions and their 
uncertainty.

BRIE2 supports the analysis of splicing processes at two molecular levels, 
either between alternative splicing isoforms or between unspliced and spliced 
RNAs. In either case, it returns cell-by-event or cell-by-gene matrices of PSI 
value and its 95% confidence interval (quantification) and the statistics for 
detecting DAS and DMG on each event or gene:

1. **Differential alternative splicing (DAS):** This task is to quantify the 
   proportions of alternative splicing isoforms and to detect DAS between groups
   of cells or along with a continuous covariate, e.g., pseudotime. 
   BRIE2 is designed for two-isoform splicing events with a focus on exon 
   skipping, but in principle also applicable for mutual exclusion, 
   intron-retaining, alternative poly-A site, 3' splice site and 5' splice site.

2. **Differential momentum genes (DMG):** This task is to quantify the 
   proportions of unspliced and spliced RNAs in each gene and each cell. 
   Similar to DAS, the DMG is a principled selection of genes that capture 
   heterogeneity in transcriptional kinetics between cell groups, e.g., cell 
   types, or continuous cell covariates, hence may enhance the RNA velocity 
   analyses by focusing on dynamics informed genes.


Installation
============

BRIE2 is available through PyPI_. To install, type the following command 
line, and add ``-U`` for upgrading:

.. code-block:: bash

  pip install -U brie

Alternatively, you can install from this GitHub repository for the latest (often 
development) version with the following command line

.. code-block:: bash

  pip install -U git+https://github.com/huangyh09/brie

In either case, add ``--user`` if you don't have the write permission for your 
Python environment.

For more instructions, see the installation_ manual.

.. _PyPI: https://pypi.org/project/brie
.. _installation: https://brie.readthedocs.io/en/latest/install.html


Manual and examples
===================

* The full manual is at https://brie.readthedocs.io 
* More examples and tutorials: https://github.com/huangyh09/brie-tutorials

In short, there are two steps for running BRIE2. 
First, obtain cell-by-gene or cell-by-event count matrices for each isoform. 
For the exon-skipping event, you can run ``brie-count``, which will return count 
matrices and hdf5 file for AnnData. 
For spliced and unspliced matrices, we listed a few options in the manual_.

Then you can use ``brie-quant`` to perform quantification of splicing ratio and 
detect differential alternative splicing or differential momentum genes. 

Type command line ``brie-count -h`` and ``brie-quant -h`` to see the full 
arguments.


.. _manual: https://brie.readthedocs.io/en/latest/quick_start.html#step1-read-counts


References
==========

* Yuanhua Huang and Guido Sanguinetti. `BRIE2: computational identification of 
  splicing phenotypes from single-cell transcriptomic experiments
  <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02461-5>`_.
  \ **Genome Biology**\, 2021; 22(1):251.

* Yuanhua Huang and Guido Sanguinetti. `BRIE: transcriptome-wide splicing 
  quantification in single cells 
  <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1248-5>`_. 
  \ **Genome Biology**\, 2017; 18(1):123.
