Simulator for BRIE
==================

The all three simulation experiments in BRIE paper are based on Spanki_ RNA-seq 
reads simulator, which supports specific RPK (reads per kilo-base) value for 
each transcript.

Here, we provide two Python wraps that were used to perform the simulation 
experiments, and can be used to replicate them easily: 

1) simulation of exon-inclusion ratio at different coverages, and generating an
   auxiliary variable for learning informative prior; 

2) simulation of drop-out of given transcripts.


Simulation of exon-inclusion ratio
----------------------------------
In this simulator ``simuPSI.py``, we provide three modes for generating 
exon-inclusion ratio, (i.e., Percent of Spliced Inclusion, PSI), and an 
auxiliary variable for learning informative prior.

1. mode=LogitNormal

  This mode assumes that PSI follows a logitnormal_ distribution. The built-in 
  mean is 0, i.e., PSI=0.5, and the default standard deviation theta=3.0. 

2. mode=UniDiff1

  This mode assumes uniform distribution of PSI across all splicing events, with
  PSI values in order of [0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9].
  This mode can be used for generating reads for condition 1 in differential 
  splicing calling. 

3. mode=UniDiff2

  This mode also assumes uniform distribution of PSI, which has opposite PSI 
  values of UniDiff1, as [0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1]. It is for 
  reads in condition 2 in differential splicing calling.

In all above modes, ``rpk`` values for Spanki input will be generated in a file 
``tran_rpk.txt``, then the Spanki software will be run automatically, and reads 
in fastq file and output truth (slight different from input rpk values, 
``transcript_sims.txt``), will be generated, based on which an auxiliary variable 
will simulated in ``prior_fractions.txt``. Parameter ``priorR`` can be used to 
customize correlations between the auxiliary variable and the simulation truth.


Simulation of drop-out for single cell
--------------------------------------
In this simulator ``simuDropout.py``, we use the expression profile (FPKM or 
TPM) and drop-out probability profile for all transcript to simulate the RPK in 
single-cell RNA-seq reads. This simulator allows adjust the overall drop-out 
rate, but keep the similarity of the drop-out probability profile, as it 
achieved by add intercept to the drop-out probability in its logit space. There 
might be slight difference between the setting drop-out rate and the final 
output.

Here, as the expression profile is set esemble to the real data, the sequence 
features from the real data could be used to learn an informative prior.

This simulator for drop-out is not limited to the exon-skipping events, but also
can be used for the whole transcriptome in single cell, even one don't look into
splicing.


Installation
------------
You could follow the manul of Spanki_ to install its original version. But here, 
we suggest you to install the version we modified to support input `rpk` or `cov`
value in float format.

- download it here and unzip it:

  * GitHub repo: https://github.com/huangyh09/Spanki

  * Zip file: https://github.com/huangyh09/Spanki/archive/master.zip

- install via (in Spanki folder):

  * ``python setup.py install``, or 

  * ``python setup.py install --user`` if you don't have root permission.

  * check if installation is successful by typing ``spankisim_transcripts``

- run the simulator wrap and see help (go to the simulator folder)
  
  * ``python simuPSI.py -h``

  * ``python simuDropout.py -h``

**Notes** The fasta file downloaded from GENCODE is not compatible with Spanki,
because the pyfasta package used in Spanki require ``>chr1`` rather than 
``>chr1 extra information``. Thus, the fasta file needs a bit change, which can
be done with spanki_anno.py_. In addition, from our experience, we suggests 
generating reads in a Python2 environment, and specifiy versions of 
dependencies of Spanki, which you could do it as follows

.. code-block:: bash

  conda create --name brieSim python=2.7 numpy=1.13.0 
  conda activate brieSim

  pip install pysam===0.14.1
  conda install -c bioconda cufflinks=2.2.1
  conda install -c bioconda samtools=1.2.0

  git clone https://github.com/huangyh09/Spanki.git
  cd Spanki && python setup.py install && cd ..
  git clone https://github.com/huangyh09/brie.git
  cd brie && python setup.py install && cd ..

  python Spanki/extra/spanki_anno.py -f GRCm38.p4.genome.fa

Then you will have dependencies and ``GRCm38.p4.genome.spanki.fa``. 

**Also, we usually use the fastq file rather than the aligned bam file from Spanki.**


Examples
--------
In the BRIE paper, the simulation is performed as following command lines:

- PSI in LogitNormal for assessing quantification accuracy
  ::

    python simuPSI.py -a $anno_file -f $ref_file -o $out_dir --rpk $rpk --seed 0 --mode LogitNormal --theta=3.0 -m errorfree

- PSI in uniform for calling differential splicing
  ::

    python simuPSI.py -a $anno_file -f $ref_file -o $out_dir --rpk $rpk --seed 0 --mode UniDiff1 -m errorfree

  ::

    python simuPSI.py -a $anno_file -f $ref_file -o $out_dir --rpk $rpk --seed 0 --mode UniDiff2 -m errorfree

- PSI with drop-out for imputation
  ::

    python simuDropout.py -a $anno_file -f $ref_file -d $dice_file -o $out_dir --dropoutProb $prob_file --dropoutRate $rate -N 600000 -m errorfree

- demo files with direct Spanki be found in the examples: spanki_demo.sh_

.. _spanki_anno.py: https://github.com/huangyh09/Spanki/blob/master/extra/spanki_anno.py
.. _Spanki: http://www.cbcb.umd.edu/software/spanki/
.. _logitnormal: https://en.wikipedia.org/wiki/Logit-normal_distribution
.. _spanki_demo.sh: https://github.com/huangyh09/brie/tree/master/simulator/spanki_demo.sh