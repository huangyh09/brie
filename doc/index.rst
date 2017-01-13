====
Home
====

.. :Author: Yuanhua Huang
.. :Version: 0.2.4
.. :Last viewed: Dec 17, 2016

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

2. ``brie-diff``: Calculate Bayes factor of differential splicing between two 
cells or two conditions. 

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

Yuanhua Huang and Guido Sanguinetti. `Transcriptome-wide splicing quantification 
in single cells <http://biorxiv.org/content/early/2017/01/05/098517>`_. 
\ **bioRxiv**\, 2017, 0985172.



.. seealso::

   Code: GitHub latest version
      https://github.com/huangyh09/brie

   Data: splicing events annotations
      http://sourceforge.net/projects/brie-rna/files/annotation/

   Data: usage examples
      http://sourceforge.net/projects/brie-rna/files/examples/



