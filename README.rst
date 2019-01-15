|PyPI| |Docs| |Build Status|

.. |PyPI| image:: https://img.shields.io/pypi/v/brie.svg
    :target: https://pypi.org/project/brie
.. |Docs| image:: https://readthedocs.org/projects/brie/badge/?version=latest
   :target: https://brie.readthedocs.io
.. |Build Status| image:: https://travis-ci.org/huangyh09/brie.svg?branch=master
   :target: https://travis-ci.org/huangyh09/brie


BRIE: Bayesian Regression for Isoform Estimate
==============================================

About BRIE
----------

BRIE (Bayesian regression for isoform estimate) is a Bayesian method to 
estimate isoform proportions from RNA-seq data. Currently, BRIE could take 
sequence features to automatically learn informative prior of exon inclusion 
ratio in  exon-skippiing events. This informative prior is very important when 
limited data is available. In Bulk RNA-seq experiment, we could easily increase 
the amplification to get more sequencing reads to improve the accuracy of 
isoform estimate. However, in single cell RNA-seq (scRNA-seq) experiments, the 
initial molecular is very limited, which always results some genes with very 
low coverage or even drop-out. In scRNA-seq, the BRIE method, by integrating 
informative prior, e.g. learned from sequence feature, could provide accurate 
and reproducible estimates of splicing in single cells, as well as sensitive 
differential analyses.


BRIE provides following functions through command line:

1. ``brie``: Estimate isoform proportions and FPKM, and calculate weights for 
regulatory features.

2. ``brie-diff``: Calculate Bayes factor of differential splicing between 
multiple cells pair-wisely. 

Quick Start
-----------

**Installation**: 

- ``pip install brie``
- or download this repository, and type ``python setup.py install``; 
- add ``--user`` if you don't have root permission and you don't use Anaconda_.

.. _Anaconda: https://www.continuum.io/anaconda-overview

**Arguments**

- Type command line ``brie -h``


Detailed manual
---------------

See the documentation_ on how to install, to use, to find the annotation data 
etc.

.. _documentation: https://brie.readthedocs.io


References
----------

Yuanhua Huang and Guido Sanguinetti. `BRIE: transcriptome-wide splicing quantification in single cells <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1248-5>`_. 
\ **Genome Biology**\, 2017; 18(1):123.
