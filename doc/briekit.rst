===============
BRIE-kit manual
===============

In this manual, we will introduce how to install, run ``briekit-event``, 
``briekit-event-filter``, and ``briekit-factor``. You can find examples below 
or in bash file anno_human.sh_ and anno_mouse.sh_  in the example folder of the 
BRIE-kit repository.

.. _anno_human.sh: https://github.com/huangyh09/briekit/blob/master/example/anno_human.sh
.. _anno_mouse.sh: https://github.com/huangyh09/briekit/blob/master/example/anno_mouse.sh


Installation
============

We recommend using Python in Anaconda platform, which gives everything you need
in one folder, and you have the permission to change files even you're not root.

BRIE-kit is developed under Python 2.7 environment, and not full compatible 
with Python 3, **so please use it in Python2 environment**. We recommend you to 
create conda_ environment to get Python2.7 as following command lines. Of 
course, you can install Ananconda2 to get a default Python 2 environment, but 
we recommend the conda_ environment, as you only do this preprossing once.

.. code-block:: bash

  conda create -n briekit python=2.7 numpy=1.13.0 

  source activate briekit

Once you are in a Python 2 environment, there are usually two ways to isntall a
package: 

- Opt 1: Type in terminal: ``pip install briekit``. Add ``-U`` if you want to 
  upgrade your earlier installation.
- Opt 2: Download the GitHub repository_, and type ``python setup.py install``

Note, if you don't use Anaconda_  and don't have root permission, add 
``--user``, so you can install it in your folder.

Sometimes, you may need to install pysam separatly (hopefully not). In our test
pysam=0.10-0.14 works fine.


.. _conda: https://conda.io/docs/user-guide/tasks/manage-environments.html
.. _Anaconda: https://anaconda.org
.. _repository: https://github.com/huangyh09/briekit

1. Splicing events generation
=============================

.. note::
  This function is not compatible with Python 3.

``briekit-event`` for generating from full annotation. This program is modified 
from Yarden Katz's Python package rnaseqlib_, with supporting different input 
annotation formats, e.g., gtf, gff3 and ucsc table. For example, you could 
download a full annotation file for mouse from GENCODE_. Then, you can download
the gene annotation and generate the splicing event by the following command:

.. _rnaseqlib: https://github.com/yarden/rnaseqlib
.. _GENCODE: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gtf.gz


.. code-block:: bash

  cd $DATA_DIR

  # download gene annotation
  wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gff3.gz

  briekit-event -a $anno_ref -o $DATA_DIR/AS_events


Then in the `$DATA_DIR/AS_events` folder, the skipping-exon events, i.e., 
``SE.gff3.gz`` will be generated.

Check more arguments bu ``briekit-event -h``.

.. note::
  If the directory (BIN_DIR) of executable ``briekit-event`` is not in 
  PATH variable, use ``$BIN_DIR/briekit-event`` rather than ``briekit-event``. 
  The same for following functions.


2. Splicing events filtering
============================

As the annotation file is not perfect, there may be false splicing events 
generated from above command line. Therefore, it can be useful to add some 
quality control on these splicing events. Here, we provide another function 
``briekit-event-filter`` to only keep high-quality events, and also add 
informative ids (gene id / transcript id). Based on above ``SE.gff3.gz``, we 
could select the gold-quality splicing event by following command line. Note, 
the reference genome sequence is also required, for example, mouse genome_ 
sequence here.

.. _genome : ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/GRCm38.p5.genome.fa.gz


Here are examples to generate default filtering and lenient filtering events.

.. code-block:: bash

  cd $DATA_DIR

  # download genome reference
  wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/GRCm38.p5.genome.fa.gz
  gzip -d GRCm38.p5.genome.fa.gz

  # default filtering
  briekit-event-filter -a AS_events/SE.gff3.gz --anno_ref=gencode.vM12.annotation.gtf.gz -r GRCm38.p5.genome.fa

  # lenient filtering
  briekit-event-filter -a AS_events/SE.gff3.gz --anno_ref=gencode.vM12.annotation.gtf.gz \
  -r GRCm38.p5.genome.fa -o AS_events/SE.lenient.gtf --add_chrom chrX,chrY --as_exon_min 10 \
  --as_exon_max 100000000 --as_exon_tss 10 --as_exon_tts 10 --no_splice_site #--keep_overlap


Then you will find an output file as ``AS_events/SE.filtered.gff3.gz``, which 
only contains splicing events passing the following (default) constrains:

* located on autosome and input chromosome
* not overlapped by any other AS-exon
* surrounding introns are no shorter than a fixed length, e.g., 100bp
* length of alternative exon regions, say, between 50 and 450bp
* with a minimum distance, say 500bp, from TSS or TTS
* surrounded by AG-GT, i.e., AG-AS.exon-GT

Check more arguments for events filtering by ``briekit-event -h``.



3. Sequence features
====================

With the splicing annotation file, a set of short sequence feature can be 
calculated by command line ``briekit-factor``. Besides the annotation file, 
it also requires genome sequence file (the same as above), and a phast_ 
conservation file in bigWig_ format. For human and mouse, you could 
download it directly from UCSC browser: mm10.60way.phastCons.bw_ 
and hg38.phastCons100way.bw_. 

.. _phast: http://compgen.cshl.edu/phast/
.. _bigWig: https://genome.ucsc.edu/goldenpath/help/bigWig.html
.. _mm10.60way.phastCons.bw: http://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/
.. _hg38.phastCons100way.bw: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/

.. note::
  In order to fetch data from the bigWig file, we use a utility ``bigWigSummary``
  that is provided from UCSC. You could download the binary file for linux from 
  here: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary

Here is example to download the executable ``bigWigSummary`` to 
``/usr/local/bin``. Of course, you can download it to anywhere you want, e.g.,
the same directory to the data.

.. code-block:: bash

  cd /usr/local/bin
  wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary
  chmod +x bigWigSummary


In order to tell ``briekit-factor`` this directory, you could either add this
directory into PATH variable, or use it as an arguments of ``briekit-factor``
by ``--bigWigSummary /usr/local/bin/bigWigSummary``. If you prefer to add it 
to PATH, add this line to ``export PATH="~/ucsc:$PATH"`` into the ``.profile`` 
or ``.bashrc`` file.


Then, you could get the sequence features by ``briekit-factor``, for example, 

.. code-block:: bash

  cd $DATA_DIR

  #download phastCon file
  wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw

  briekit-factor -a AS_events/SE.filtered.gff3.gz -r GRCm38.p5.genome.fa -c mm10.60way.phastCons.bw -o mouse_features.csv -p 10 --bigWigSummary ./bigWigSummary

Then you will have the features stored in a ``mouse_features.csv.gz`` file, 
where #`factors` * #`gene_ids` features values are saved.

Check more arguments for fetch sequence features by ``briekit-factor -h``.


