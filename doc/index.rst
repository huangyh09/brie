====
Home
====

.. :Author: Yuanhua Huang
.. :Version: 0.1.5
.. :Last viewed: Jun 3, 2018

About BRIE
==========

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

3. ``brie-event``: Extract the splicing events from gene annotation file in 
gff3/gtf formate. 

4. ``brie-factor``: Fetch genentic features from genome sequence reference file 
in fasta formate.


Contents
========

.. toctree::
   :maxdepth: 2

   install.rst
   manual.rst
   faq.rst
   release.rst



References
==========

Yuanhua Huang and Guido Sanguinetti. `BRIE: transcriptome-wide splicing 
quantification in single cells 
<https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1248-5>`_. 
\ **Genome Biology**\, 2017; 18(1):123.



.. seealso::

   Code: GitHub latest version
      https://github.com/huangyh09/brie

   Simulation wraps on GitHub
      https://github.com/huangyh09/brie/tree/master/simulator

   Data: splicing events annotations
      http://sourceforge.net/projects/brie-rna/files/annotation/

   Data: usage examples
      http://sourceforge.net/projects/brie-rna/files/examples/



