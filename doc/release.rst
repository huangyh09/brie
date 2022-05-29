=============
Release notes
=============

Release in future
=================

Release v2.2.0 (29/05/2022)
===========================
* **Big new feature** `brie-count` now supports cell barcodes based platforms, 
  including 10x Genomics
* provide small test files in the tests folder
* Update `brie-count` manual

Release v2.1.0 (22/03/2022)
===========================
* Enable `brie-count` to support other types of splicing events, not only SE
* Fix a minor bug in exteme scenario with missing a certain isoform for all genes

Release v2.0.6 (23/09/2021)
===========================
* Fix in format bug in brie-count
* Update brie.pl.volcano for using ELBO_gain as default y-axis
* Update manual according to the revised paper
* Add CLI ``brie``

Release v2.0.5 (04/11/2020)
===========================
* Support saving detection table to tsv file
* Add the exon start and stop positions in brie-count

Release v2.0.4 (26/09/2020)
===========================
* Tune the learning rate with multiple values
* For test model, the fitted parameters will be used as initials
* Support base model with full features or null feature
* For gene feature only, switch sigma into per cell base
* Add noise term in simulator
* A few minor bug fix
* Implement a Inv-Gamma prior distribution for sigma (not in use)

Release v2.0.3 (26/08/2020)
===========================
* Support to use minimum minor isoform frequency to fileter genes (default=0.001)
* Add pseudo count (default=0.01) for none-empty element in both unique counts 
  for more robust results
* Reduce the sample size for Monte Carlo Expectation (10 to 3) for computational
  efficiency
* Restructure the arguments in brie-quant
* Initialise the example notebook on multiple sclerosis data

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
