======
Manual
======

.. contents:: Table of Contents
   :depth: 1
   :local:

After properly installed BRIE Python package, five excutable binary files could 
be run from command line directly: ``brie``, ``brie-diff``, ``brie-factor``, 
``brie-event``, and ``brie-event-filter``. Details about these utilities are 
shown below.


1. BRIE isoform estimate
========================

This is the main program to quitify the fraction of exon inclusion level. In 
order to automatically learn the informative prior, the predictive features are 
required. There are two ways to get the annotation and sequence features: 

1. use our processed annotation file and according sequence features, which you 
   can download from here_. Currently, we produced data for human_ and mouse_. 
   We suggest align RNA-seq reads to the according version of reference genome.

2. generate the annotation and fetch the sequence features with the help of 
   brie-event_ and brie-factor_ by yourself

.. _here: https://sourceforge.net/projects/brie-rna/files/annotation/
.. _human: https://sourceforge.net/projects/brie-rna/files/annotation/human/gencode.v25/
.. _mouse: https://sourceforge.net/projects/brie-rna/files/annotation/mouse/gencode.vM12/
.. _brie-event: https://brie-rna.sourceforge.io/manual.html#splicing-events
.. _brie-factor: https://brie-rna.sourceforge.io/manual.html#sequence-features


Then you could input the feature file obtained above, and run it like this:

::

  brie -a AS_events/SE.gold.gtf -s Cell1.sorted.bam -f mouse_features.csv.gz -o out_dir -p 15

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

  Usage: brie [options]

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

  brie-diff -i cell1/samples.csv.gz,cell2/samples.csv.gz -o c1_c2.diff.tsv -p 1 --minBF 0


For many cells (gives events with ``BF>10``. Speed: 100 cells in ~10min with 30 
CPUs)

::

  fileList=cell1/samples.csv.gz,cell2/samples.csv.gz,cell3/samples.csv.gz,cell4/samples.csv.gz

  brie-diff -i $fileList -o c1_c4.diff.tsv

Then you will have an output file with 15 columns:

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

There are more parameters for setting (``brie-diff -h`` always give the version 
you are using):

.. code-block:: html

  Usage: brie-diff [options]

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



3. Splicing events
==================

**Splicing events generating from full annotation**

This program is modified from Yarden Katz's Python package rnaseqlib_, with 
supporting different input annotation formats, e.g., gtf, gff3 and ucsc table.
For example, you could download a full annotation file for mouse from GENCODE_.
Then, you can generate the splicing event by the following command:

::

  brie-event -a gencode.vM12.annotation.gtf

.. _rnaseqlib: https://github.com/yarden/rnaseqlib
.. _GENCODE: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gtf.gz

Then in the same folder of the annotation file, there will a new folder 
``AS_events``, where the file for skipping-exon events, i.e., ``SE.gff3``
is located.

There are more parameters for setting (``brie-event -h`` always give 
the version you are using):

.. code-block:: html

  Usage: brie-event [options]

  Options:
    -h, --help            show this help message and exit
    -a ANNO_FILE, --anno_file=ANNO_FILE
                          The annotation files used in making the annotation.
                          You could input multiple files; use comma',' as
                          delimiter.
    --anno_type=ANNO_TYPE
                          The type of each annotation file. Use one for all
                          files or set for each file. Use comma ',' as
                          delimiter. You could choose 'ucsc', 'gtf', 'gff3'.
                          [default: gtf]
    -o OUTPUT_DIR, --output_dir=OUTPUT_DIR
                          Output directory.
    --flanking-rule=FLANKING_RULE
                          Rule to use when defining exon trios. E.g.
                          'commonshortest' to use the most common and shortest
                          regions are flanking exons to an alternative trio.
                          [default: commonshortest]
    --multi-iso           If passed, generates multi-isoform annotations. Off by
                          default.
    --genome-label=GENOME_LABEL
                          If given, used as label for genome in output files.
    --sanitize            If passed, sanitize the annotation. Off by default.



**Splicing events quality check**

As the annotation file is not perfect, there may be false splicing events 
generated from above command line. Therefore, we provide another function 
``brie-event-filter`` to only keep high-quality events, and use informative 
ids. Based on above ``SE.gff3``, we could select the gold-quality splicing 
event by following command line. Note, the reference genome sequence is also 
required, for example, mouse genome_ sequence here.

.. _genome : ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/GRCm38.p5.genome.fa.gz

::

  brie-event-filter -a AS_events/SE.gff3 -anno_ref gencode.vM12.annotation.gtf -r GRCm38.p5.genome.fa

Then you will find an output file as ``AS_events/SE.gold.gff3``, which only 
contains splicing events passing the following constrains:

* located on autosome and input chromosome
* not overlapped by any other AS-exon
* surrounding introns are no shorter than a fixed length, e.g., 100bp
* length of alternative exon regions, say, between 50 and 450bp
* with a minimum distance, say 500bp, from TSS or TTS
* surrounded by AG-GT, i.e., AG-AS.exon-GT

There are more parameters for setting (``brie-event-filter -h`` always give 
the version you are using):

.. code-block:: html

  Usage: brie-event-filter [options]

  Options:
    -h, --help            show this help message and exit
    -a ANNO_FILE, --anno_file=ANNO_FILE
                          The annotation file of SE events in gff3 format from
                          rnaseqlib.
    --anno_ref=ANNO_REF   The reference annotation file in gtf format.
    -r REFERENCE, --reference=REFERENCE
                          The genome reference sequence file in fasta format.
    -o OUT_FILE, --out_file=OUT_FILE
                          The prefix of out files.
    --as_exon_min=AS_EXON_MIN
                          the minimum length for the alternative splicing exon.
    --as_exon_max=AS_EXON_MAX
                          the maximum length for the alternative splicing exon.
    --as_exon_tss=AS_EXON_TSS
                          the minimum length for the alternative exon to TSS.
    --as_exon_tts=AS_EXON_TTS
                          the minimum length for the alternative exon to TTS.
    --add_chrom=ADD_CHROM
                          the extra chromosomes besides autosome, e.g.,
                          chrX,chrY,chrM



4. Sequence features
====================

With the splicing annotation file, a set of short sequence feature can be 
calculated by command line ``brie-factor``. Besides the annotation file, 
it also requires genome sequence file (the same as above), and a phast_ 
conservation file in bigWig_ format. For human and mouse, you could 
download it directly from UCSC browser: mm10.60way.phastCons.bw_ 
and hg38.phyloP100way.bw_. 

.. _phast: http://compgen.cshl.edu/phast/
.. _bigWig: https://genome.ucsc.edu/goldenpath/help/bigWig.html
.. _mm10.60way.phastCons.bw: http://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/
.. _hg38.phyloP100way.bw: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/

.. note::
  In order to fetch data from the bigWig file, we use a utility ``bigWigSummary``
  that is provided from UCSC. You could download the binary file for linux from 
  here: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary

Besides, we require that ``bigWigSummary`` can be access in the PATH environment. 
You can do it by the following command lines, and put into the ``.profile`` or 
``.bashrc`` file.

::

  chmod +x ~/ucsc/bigWigSummary
  export PATH="~/ucsc:$PATH"

Then, you could get the sequence features by ``brie-factor``, for example, 

::

  brie-factor -a AS_events/SE.gold.gtf -r GRCm38.p5.genome.fa -c mm10.60way.phastCons.bw -o mouse_features.csv -p 10

Then you will have the features stored in a ``mouse_features.csv.gz`` file, 
where #`factors` * #`gene_ids` features values are saved.
 
There are more parameters for setting (``brie-factor -h`` always give the 
version you are using):

.. code-block:: html

  Usage: brie-factor [options]

  Options:
    -h, --help            show this help message and exit
    -a ANNO_FILE, --anno_file=ANNO_FILE
                          Annotation file for genes and transcripts in GTF or
                          GFF3
    -r REF_SEQ, --ref_seq=REF_SEQ
                          Genome sequence reference in fasta file.
    -c PHAST_FILE, --phastCons=PHAST_FILE
                          PhastCons conservation scores in bigWig file.
    -o OUT_FILE, --out_file=OUT_FILE
                          Output in csv file, [default: brieFactor.cvs]

    Optional arguments:
      -p NPROC, --nproc=NPROC
                          Number of subprocesses [default: 4]
      --MSA5ss=MSA_5SS    Mutiple sequence alignment file for 5'splice-site. It
                          is from -4 to 7. As default, MSA is based on input 5
                          splice sites.
      --MSA3ss=MSA_3SS    Mutiple sequence alignment file for 3'splice-site. It
                          is from -16 to 4. As default, MSA is based on input 3
                          splice sites.



5. Preprocess
=============

5.1 reads alignment
-------------------

Usually, the initial RNA-seq reads is in fastq_ format, without information of 
where it comes from the genome location. BRIE, similar as DICEseq and MISO, it 
requires RNA-seq reads aligned to genome sequence. It should be in sam/bam 
format, after sorting and indexing.

There are quite a fewer aligner that allows mapping reads to genome reference 
with big gaps, mainly caused by splicing. For example, you could use STAR_ and 
HISAT_, which usually return good alignment quality.

You could run it like this (based on HISAT 0.1.5), which including alignment, 
sort and index:

::

  ($hisatDir/hisat -x $hisatRef -1 $fq_dir/"$file"_1.fq.gz -2 $fq_dir/"$file"_2.fq.gz --no-unal | samtools view -bS -> $out_dir/$file.bam) 2> $out_dir/$file.err
  samtools sort $out_dir/$file.bam $out_dir/$file.sorted
  samtools index $out_dir/$file.sorted.bam

.. _fastq: https://en.wikipedia.org/wiki/FASTQ_format
.. _STAR: https://code.google.com/p/rna-star/
.. _HISAT: https://ccb.jhu.edu/software/hisat/index.shtml


6. Examples
===========

There are some examples available here: 
https://sourceforge.net/projects/brie-rna/files/examples/

- Example to quantify splicing with provided annotation (bash code and data): 
  brie-examples.zip_

- Example to quantify splicing with provided annotation (bash code): 
  brie_demo.sh_

- Example to generate splicing events and fetch sequence factors (bash codes): 
  anno_maker.sh_

.. _brie-examples.zip: http://ufpr.dl.sourceforge.net/project/brie-rna/examples/brie_quantify/brie-examples.zip
.. _brie_demo.sh: https://github.com/huangyh09/brie/blob/master/example/brie_demo.sh
.. _anno_maker.sh: https://github.com/huangyh09/brie/blob/master/example/anno_maker.sh

