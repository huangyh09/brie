============
BRIE1 Manual
============

After properly installed BRIE Python package, two excutable binary files could 
be run from command line directly: ``brie1``, ``brie1-diff``. From version 0.2.0, 
all preprocessing are divided and moved into BRIE-kit_ package, which is aimed
to be used in Python2 only. 

.. _BRIE-kit: https://github.com/huangyh09/briekit/wiki

1. BRIE isoform estimate
========================

This is the main program to quitify the fraction of exon inclusion level. In 
order to automatically learn the informative prior, the predictive features are 
required. There are two ways to get the annotation and sequence features: 

1. use our processed annotation file and according sequence features, which you 
   can download from here_. Currently, we produced data for human_ and mouse_. 
   We suggest align RNA-seq reads to the according version of reference genome.

2. generate the annotation and fetch the sequence features with the help of 
   briekit-event_ and briekit-factor_ by yourself

.. _here: https://sourceforge.net/projects/brie-rna/files/annotation/
.. _human: https://sourceforge.net/projects/brie-rna/files/annotation/human/gencode.v25/
.. _mouse: https://sourceforge.net/projects/brie-rna/files/annotation/mouse/gencode.vM12/
.. _briekit-event: https://brie-rna.sourceforge.io/manual.html#splicing-events
.. _briekit-factor: https://brie-rna.sourceforge.io/manual.html#sequence-features


Then you could input the feature file obtained above, and run it like this:

::

  brie1 -a AS_events/SE.gold.gtf -s Cell1.sorted.bam -f mouse_features.csv.gz -o out_dir -p 15

By default, you will have three output files in the out_dir: ``fractions.tsv``, 
``weights.tsv`` and ``samples.csv.gz``. 

- In ``fractions.tsv``, there are 8 columns:

  * column 1: transcript id
  * column 2: gene id
  * column 3: transcript length
  * column 4: reads counts for whole events
  * column 5: FPKM for each isoform
  * column 6: fraction for each isoform, called Psi
  * column 7: lower bound of 95% confidence interval of isoform fraction
  * column 8: higher bound of 95% confidence interval of isoform fraction

- In ``weights.tsv``, there are the weights for the Bayesian regression, with 
  `#Feature+2` lines, involving each features, interpret and sigma (a hyperparameter). 
  There are two columns each line, including the label and the value.

- In ``sample.csv.gz``, there are the MCMC_ samples of posterior distribution of 
  Psi. These samples are used to detect the differential splicing.

.. _MCMC: https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo

There are more parameters for setting (``brie -h`` always give the version you 
are using)

.. code-block:: html

  Usage: brie1 [options]

  Options:
    -h, --help            show this help message and exit
    -a ANNO_FILE, --anno_file=ANNO_FILE
                          Annotation file for genes and transcripts in GTF or
                          GFF3
    -s SAM_FILE, --sam_file=SAM_FILE
                          Sorted and indexed bam/sam files, use ',' for
                          replicates e.g., rep1.sorted.bam,sam1_rep2.sorted.bam
    -o OUT_FILE, --out_file=OUT_FILE
                          Prefix of the output files with full path
    -f FACTOR_FILE, --factor_file=FACTOR_FILE
                          Features in csv.gz file to predict isoform expression.

    Optional arguments:
      -p NPROC, --nproc=NPROC
                          Number of subprocesses [default: 4]
      -w WEIGHT_FILE, --weight_file=WEIGHT_FILE
                          File with weights, an output of Brie.
      -y FTYPE, --ftype=FTYPE
                          Type of function target: FPKM, Y, Psi [default: Y].
      --fLen=FRAG_LENG    Two arguments for fragment length: mean and standard
                          diveation, default: auto-detected
      --bias=BIAS_ARGS    Three argments for bias correction:
                          BIAS_MODE,REF_FILE,BIAS_FILE(s). BIAS_MODE: unif,
                          end5, end3, both. REF_FILE: the genome reference file
                          in fasta format. BIAS_FILE(s): bias files from dice-
                          bias, use '---' for time specific files, [default:
                          unif None None]
      --sigma=_SIGMA      Sigma in Bayesian regression: the Gaussian standard
                          deviation of residues [default: Auto].
      --lambda=_LAMBDA    Lambda in Bayesian regression: the coeffiecient of L2
                          constrain on weights [default: 0.1].
      --mcmc=MCMC_RUN     Four arguments for in MCMC iterations:
                          save_sample,max_run,min_run,gap_run. Required:
                          save_sample =< 3/4*mim_run. [default: 500 5000 1000 50]

**Hyperparamers**

* ``sigma`` is the square rooted variance of Gaussian noise in Bayesian 
  regression. By default, it will learn it automatically. Alternatively, you 
  could set it with your experience, for example, 3 might be a good option. 
* ``lambda`` is the constrain on weights of Bayesian regression. 0.1 is good 
  option in ENCODE data.
* ``weight_file`` is fixed weights for Bayesian regression. Therefore, the 
  prior is predicted from the input weight file and its sequence features.
  


2. Differential splicing
========================

This command allows to detect differential splicing between many cells 
pair-wisely, including just two cells, by calculating Bayes factor. You could 
run it as follows:

For two cells (``-p 1 --minBF 0`` gives all events in the same order. Speed: 
10-20 second with 1 CPU)

::

  brie1-diff -i cell1/samples.csv.gz,cell2/samples.csv.gz -o c1_c2.diff.tsv -p 1 --minBF 0


For many cells (gives events with ``BF>10``. Speed: 100 cells in ~10min with 30 
CPUs)

::

  fileList=cell1/samples.csv.gz,cell2/samples.csv.gz,cell3/samples.csv.gz,cell4/samples.csv.gz

  brie1-diff -i $fileList -o c1_c4.diff.tsv

Then you will have two output files. The first one (in the format of xxx.diff.tsv) 
contains all Bayes factor passing the threshold, and it has with 15 columns:

* column1-2: transcript id and gene id
* column3-4: cell 1 and cell 2 names (the folder names)
* column5-6: prior of exon inclusion fraction for cell 1 and cell 2
* column7-8: posterior of exon inclusion fraction for cell 1 and cell 2
* column9-12: counts for inclusion and exclusion for cell1, and then cell 2
* column13-14: probability of prior and posterior diff<0.05
* column 15: Bayes factor

.. note::
  Bayes factor is different from p value in hypothesis test. A good threshold 
  could be ``Bayes factor > 10`` as differential splicing event between two 
  cells.

Also another file ranks these splicing events by the number of cell paris with
differential splicing. It has 4 columns: ``gene_id``, ``cell_pairs``, 
``mean_BF``, ``median_BF``.

There are more parameters for setting (``brie-diff -h`` always give the version 
you are using):

.. code-block:: html

  Usage: brie1-diff [options]

  Options:
  -h, --help            show this help message and exit
  -i IN_FILES, --inFiles=IN_FILES
                        Input files of Brie samples for multiple cells, comma
                        separated for each cell, e.g., cell1,cell2,cell3
  -o OUT_FILE, --outFile=OUT_FILE
                        Output file with full path

  Optional arguments:
    -p NPROC, --nproc=NPROC
                        Number of subprocesses [default: 4]
    -n BOOTSTRAP, --bootstrap=BOOTSTRAP
                        Number of bootstrap [default: 1000]
    --minBF=MINBF       Minimum BF for saving out, e.g., 3 or 10. If it is 0,
                        save all events [default: 10]





3. Examples
===========

One typical example on 130 mouse cells during gastrulation is in this folder, 
from which you will quantify the splicing with BRIE, identify the highly 
variable splicing events and visualise them with sashimi plot.
https://github.com/huangyh09/brie/tree/0.2.x/examples/gastrulation


There are some earlier examples: 
https://sourceforge.net/projects/brie-rna/files/examples/

- Example to quantify splicing with provided annotation (bash code and data): 
  brie-examples.zip_

- Example to quantify splicing with provided annotation (bash code): 
  brie_demo.sh_

- Example to generate splicing events and fetch sequence factors (bash codes): 
  anno_maker.sh_

.. _brie-examples.zip: http://ufpr.dl.sourceforge.net/project/brie-rna/examples/brie_quantify/brie-examples.zip
.. _brie_demo.sh: https://github.com/huangyh09/brie/blob/0.2.x/examples/brie_demo.sh
.. _anno_maker.sh: https://github.com/huangyh09/brie/blob/0.2.x/examples/anno_maker.sh

