==============
brie-count CLI
==============

BRIE (>=2.0.0) provides two CLI directly available in your Python path: 
``brie-count``, ``brie-quant``. 

From **v2.1**, the ``brie-count`` CLI supports any type of splicing events, 
with primarily focusing on exon-skipping events:

* exon skipping (SE): we provide a curated annotation for human_ and mouse_; the
  output is in h5ad format, seamlessly for downstream analysis with 
  ``brie-quant``
  
* other types: in generally, ``brie-count`` returns reads (or UMIs) counts for
  that are compatible to isoform groups (regardless of the number of exons it 
  contains). However, the downstream analysis (e.g., differentail splicing) 
  with ``brie-quant`` mainly supports two-isoform based events. Also, users may 
  need to curate the splicing annotation themselves (we may provide some 
  curation in future).

From **v2.2**, the ``brie-count`` CLI supports counting on both well-based 
platform (e.g., smart-seq2) where each cell has a separate bam (sam/cram) file 
and droplet-based platform (e.g., 10x Genomics) where all cells are pooled in 
one big bam file through cell barcodes.


Splicing annotations
====================
As input, it requires a list of annotated splicing events. We have generated
annotations for human_ and mouse_. If you are using it, please align RNA-seq 
reads to the according version (or close version) of reference genome. 
Alternatively, you can use `briekit`_ package to generate.

.. _human: https://sourceforge.net/projects/brie-rna/files/annotation/human/gencode.v25/
.. _mouse: https://sourceforge.net/projects/brie-rna/files/annotation/mouse/gencode.vM12/
.. _briekit: https://github.com/huangyh09/briekit/wiki

Read counting
=============

The ``brie-count`` CLI calculates a count tensor for the number of reads that 
are aligned to each splicing event and each cell, and stratifies four them into
four different categories (using SE as an example):

1. key 1: Uniquely aligned to isoform1, e.g., in exon1-exon2 junction in SE event
2. key 2: Uniquely aligned to isoform2, e.g., in exon1-exon3 junction in SE event
3. key 3: Ambiguously aligned to isoform1 and isoform2, e.g., within exon1
4. key 0: Partially aligned in the region but not compatible with any of the two 
   isoforms. We suggest ignoring these reads.

.. note::
   You can decode isoform compatibility from the key index by using the 
   `decimal-bo-binary 
   <https://docutils.sourceforge.io/docs/user/rst/quickref.html#internal-hyperlink-targets>`_
  
  For example, the key 5 (decimal) can be decoded to 101 (binary), which means
  compatible with both isoform 1 (the last) and 3 (the third from right).

Then you fetch the counts on a list of bam files by the command line like this:

.. code-block:: bash

  # for smart-seq
  brie-count -a AS_events/SE.gold.gtf -S sam_and_cellID.tsv -o out_dir -p 15

  # for droplet, e.g. 10x Genomics
  brie-count -a AS_events/SE.gold.gtf -s possorted.bam -b barcodes.tsv.gz -o out_dir -p 15


By default, you will have four output files in the out_dir: ``brie_count.h5ad``, 
``read_count.mtx.gz``, ``cell_note.tsv.gz``, and ``gene_note.tsv.gz``. The 
``brie_count.h5ad`` contains all information for downstream analysis, e.g., for
`brie-quant`.

Options
=======

There are more parameters for setting (``brie-count -h`` always give the version 
you are using):

.. code-block:: html

    Usage: brie-count [options]

    Options:
      -h, --help            show this help message and exit
      -a GFF_FILE, --gffFile=GFF_FILE
                            GTF/GFF3 file for gene and transcript annotation
      -o OUT_DIR, --out_dir=OUT_DIR
                            Full path of output directory [default: $samFile/brieCOUNT]

      SmartSeq-based input:
        -S SAMLIST_FILE, --samList=SAMLIST_FILE
                            A no-header tsv file listing sorted and indexed bam/sam/cram
                            files. Columns: file path, cell id (optional)

      Droplet-based input:
        -s SAM_FILE, --samFile=SAM_FILE
                            One indexed bam/sam/cram file
        -b BARCODES_FILE, --barcodes=BARCODES_FILE
                            A file containing cell barcodes without header
        --cellTAG=CELL_TAG  Tag for cell barocdes [default: CB]
        --UMItag=UMI_TAG    Tag for UMI barocdes [default: UR]

      Optional arguments:
        --verbose           Print out detailed log info
        -p NPROC, --nproc=NPROC
                            Number of subprocesses [default: 4]
        -t EVENT_TYPE, --eventType=EVENT_TYPE
                            Type of splicing event for check. SE: skipping-exon; Any: no-
                            checking [default: SE]
