=======
History
=======


Release v0.1.3 (16/04/2017)
===========================
* ``brie-diff`` takes multiple cells, which could handle pair-wise comparisons 
  for 100 cells in ~10min with 30 CPUs; and 1000 cells within a few hours.
* Simulation wraps on Spanki are provided for simulating splicing events at 
  different coverages or drop-out probability and drop-out rate for single cells.
  https://github.com/huangyh09/brie/tree/master/simulator

Release v0.1.2 (13/01/2017)
===========================
* support Python 3.x now
* do not depend on h5py anymore for hdf5 data storage.
* ``brie-factor`` returns xxx.csv.gz rather than xxx.h5
* ``brie`` returns sample.csv.gz rather than sample.h5
* ``brie-diff`` takes sample.csv.gz rather than sample.h5

Release v0.1.1 (02/01/2017)
---------------------------
* change licence to Apache License v2
* update ``brie-event-filter``

Release v0.1.0 (29/12/2016)
---------------------------
* Initial release of BRIE
