=============
Release notes
=============

Release v2.0.2 (24/08/2020)
===========================
* Fix minor bugs in brei-count and brie-quant cli for compatibility

Release v2.0.0 (23/08/2020)
===========================
* Change the whole BRIE model from MCMC sampler (v1) to variational inference (v2)
* Change the usage of each read to a summarised read counts for speedup
* Support splicing quantification without any feature or cell features or gene 
  features or both types of features.
* Support detection of variable splicing event associated with cell features
* Support acceleration with Graphic card (Nvidia GPU)
* Compatible with Scanpy analysis suite with a variety of plotting functions
* Restructure the whole package
* BRIE earlier version is still avaible with `brie1`


Release v0.2.4 (21/10/2019)
===========================
* fix a bug that fragment length longer than transcript length

Release v0.2.2 (15/01/2019)
===========================
* move __version__ into a separate file, so import brie is not required before
  installation.
* support cram file as input aligned reads

Release v0.2.0 (03/06/2018)
===========================
* split the pre-processing functions in BRIE to another separate package 
  BRIE-kit, as some functions in the pre-processing require Python2 environment.
  So, `BRIE` will keep two functions `brie` and `brie-diff`, while `BRIE-kit` 
  will have `briekit-event`, `briekit-event-filter`, and `briekit-factor`.
  See BRIE-KIT: https://github.com/huangyh09/briekit/wiki

Release v0.1.5 (03/06/2018)
===========================
* support gff in gz format
* add an example with 130 cells
* `brie-diff` supporting ranking genes

Release v0.1.4 (02/06/2018)
===========================
* turn off pylab
* fix a bug for function `rev_seq`, reported in issue #10
* update documentation

Release v0.1.3 (16/04/2017)
===========================
* ``brie-diff`` takes multiple cells, which could handle pair-wise comparisons 
  for 100 cells in ~10min with 30 CPUs; and 1000 cells within a few hours.
* Simulation wraps on Spanki are provided for simulating splicing events at 
  different coverages or drop-out probability and drop-out rate for single 
  cells: https://github.com/huangyh09/brie/tree/master/simulator

Release v0.1.2 (13/01/2017)
===========================
* support Python 3.x now
* do not depend on h5py anymore for hdf5 data storage.
* ``brie-factor`` returns xxx.csv.gz rather than xxx.h5
* ``brie`` returns sample.csv.gz rather than sample.h5
* ``brie-diff`` takes sample.csv.gz rather than sample.h5

Release v0.1.1 (02/01/2017)
===========================
* change licence to Apache License v2
* update ``brie-event-filter``

Release v0.1.0 (29/12/2016)
===========================
* Initial release of BRIE
