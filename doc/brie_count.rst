==============
brie-count CLI
==============

BRIE (>=2.0.0) provids two CLI directly available in your Python path: 
``brie-count``, ``brie-quant``. 

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
four different categories in :

1. key 1: Uniquely aligned to isoform1, e.g., in exon1-exon2 junction in SE event
2. key 2: Uniquely aligned to isoform2, e.g., in exon1-exon3 junction in SE event
3. key 3: Ambiguously aligned to isoform1 and isoform2, e.g., within exon1
4. key 0: Partially aligned in the region but not compatible with any of the two 
   isoforms. We suggest ignoring these reads.

Then you fetch the counts on a list of bam files by the command line like this:

.. code-block:: bash

  brie-count -a AS_events/SE.gold.gtf -S sam_and_cellID.tsv -o out_dir -p 15

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
      -S SAMLIST_FILE, --samList=SAMLIST_FILE
                            A tsv file containing sorted and indexed bam/sam/cram 
                            files. No header line; file path and cell id (optional)
      -o OUT_DIR, --out_dir=OUT_DIR
                            Full path of output directory [default: $samList/brieCOUNT]

      Optional arguments:
        -p NPROC, --nproc=NPROC
                            Number of subprocesses [default: 4]
