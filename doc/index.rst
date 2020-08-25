|PyPI| |Docs| |Build Status|

.. |PyPI| image:: https://img.shields.io/pypi/v/brie.svg
    :target: https://pypi.org/project/brie
.. |Docs| image:: https://readthedocs.org/projects/brie/badge/?version=latest
   :target: https://brie.readthedocs.io
.. |Build Status| image:: https://travis-ci.org/huangyh09/brie.svg?branch=master
   :target: https://travis-ci.org/huangyh09/brie
   
====
Home
====



About BRIE
==========

Welcome to the new BRIE2 (Bayesian regression for isoform estimate, v2), a 
scalable Bayesian method to robustly identify splicing phenotypes in single 
cells RNA-seq designs and accurately estimate isoform proportions and its 
uncertainty.

BRIE2 supports isoform quantification for different needs:

1. cell features: informative prior is learned from shared cell processes. It 
   also allows to effectively detect splicing phenotypes by using Evidence Lower
   Bound gain, an approximate of Bayes factor.
   
2. gene features: informative prior is learned from shared gene regulatory 
   features, e.g., sequences and RNA protein binding

3. no feature: use zero-mean logit-normal as uninformative prior, namely
   merely data deriven
   
Note, `BRIE1 CLI`_ is still available in this version but changed to `brie1` 
and `brie1-diff`.

.. _BRIE1 CLI: https://brie.readthedocs.io/en/latest/brie1.html


Questions or Bugs
=================
If you find any error or suspicious bug, we will appreciate your report.
Please write them in the github issues: 
https://github.com/huangyh09/brie/issues

If you have questions on using BRIE, feel free get in touch with us: 
yuanhua <at> hku.hk


Quick Resources
===============

**Code: GitHub latest version**
https://github.com/huangyh09/brie

**Data: splicing events annotations**
http://sourceforge.net/projects/brie-rna/files/annotation/

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



.. toctree::
   :caption: Main
   :maxdepth: 1
   :hidden:
   
   index
   install
   quick_start
   brie_count
   brie_quant
   release

.. toctree::
   :caption: Examples
   :maxdepth: 1
   :hidden:

   multiple_sclerosis
   Prior_distribution_BRIE2
