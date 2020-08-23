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

Installation
============

BRIE2 is available through PyPI_. To install, type the following command 
line, and add ``-U`` for upgrading:

.. code-block:: bash

  pip install -U brie

Alternatively, you can install from this GitHub repository for latest (often 
development) version by following command line

.. code-block:: bash

  pip install -U git+https://github.com/huangyh09/brie

In either case, add ``--user`` if you don't have the write permission for your 
Python environment.

For more instructions, see the installation_ manual.

.. _PyPI: https://pypi.org/project/brie
.. _installation: https://brie.readthedocs.io/en/latest/install.html


Manual and examples
===================

The full manual is at https://brie.readthedocs.io 
More examples and tutorials are coming soon.

In brief, you need to run `brie-count` first, which will return a count matrix
and hdf5 file for AnnData. Then you can use `brie-quant` to perform 
quantification in different settings. Type command line ``brie-count -h`` and 
``brie-quant -h`` to see the full arguments.


References
==========

Yuanhua Huang and Guido Sanguinetti. `BRIE: transcriptome-wide splicing 
quantification in single cells 
<https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1248-5>`_. 
\ **Genome Biology**\, 2017; 18(1):123.
