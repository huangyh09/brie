====
Home
====

.. :Author: Yuanhua Huang
.. :Version: 2.0.0
.. :Last viewed: Aug 23, 2020

About BRIE
==========

BRIE2 (Bayesian regression for isoform estimate, v2) is a scalable Bayesian 
method to robustly identify splicing phenotypes in single cells RNA-seq designs 
and accurately estimate isoform proportions and its uncertainty. 

BRIE2 supports isoform quantification for different needs:

1. likelihood based: without learning informative prior from any features

2. gene features: informative prior is learned from shared gene regulatory 
   features, e.g., sequences and RNA protein binding

3. cell features: informative prior is learned from shared cell processes. It 
   also allows to effectively detect splicing phenotypes by using Evidence Lower
   Bound gain, an approximate of Bayes factor.


Quick Resources
===============

**Code: GitHub latest version**
https://github.com/huangyh09/brie

**Simulation wraps on GitHub**
https://github.com/huangyh09/brie/tree/master/simulator

**Data: splicing events annotations**
http://sourceforge.net/projects/brie-rna/files/annotation/

**Data: usage examples**
http://sourceforge.net/projects/brie-rna/files/examples/

**All releases**
https://pypi.org/project/brie/#history

**Issue reports**
https://github.com/huangyh09/brie/issues



References
==========

Yuanhua Huang and Guido Sanguinetti. `BRIE: transcriptome-wide splicing 
quantification in single cells 
<https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1248-5>`_. 
\ **Genome Biology**\, 2017; 18(1):123.

