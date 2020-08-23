==============
brie-count CLI
==============

After properly installed BRIE (>=2.0.0) Python package, two CLI will be 
available directly in your Python path: ``brie-count``, ``brie-quant``. 

In this documentation, BRIE refers to BRIE2 (>=2.0.0). For using BRIE1 (<=0.2.2)
with MCMC sampler, please refer to BRIE1_. If you want to generate splicing
annotations for your data, e.g., a species different from human and mouse,
please use a separated package BRIE-kit_, which is developed in Python2. 

.. _BRIE1: brie1.html
.. _BRIE-kit: https://github.com/huangyh09/briekit/wiki


Options
=======

This CLI will return a count tensor for the number of reads aligned to four 
different categories in each splicing event and each cell:
1) Uniquely aligned to isoform1, e.g., in exon1-exon2 junction in SE event
2) Uniquely aligned to isoform2, e.g., in exon1-exon3 junction in SE event
3) Ambiguously aligned to isoform1 and isoform2, e.g., within exon1
0) Fetched in the region but not compatible with any of the above three types.
   We suggest ignoring these reads.
   
As input, you need to generate the splicing event annotation. We have generated
data for human_ and mouse_. We suggest align RNA-seq reads to the according 
version of reference genome. Alternatively, you can use `briekit` package to 
generate.

.. _human: https://sourceforge.net/projects/brie-rna/files/annotation/human/gencode.v25/
.. _mouse: https://sourceforge.net/projects/brie-rna/files/annotation/mouse/gencode.vM12/


Then you fetch the counts on a list of bam files by the command line like this:

.. code-block:: bash

  brie-count -a AS_events/SE.gold.gtf -S sam_and_cellID.tsv -o out_dir -p 15

By default, you will have four output files in the out_dir: ``brie_count.h5ad``, 
``read_count.mtx.gz``, ``cell_note.tsv.gz``, and ``gene_note.tsv.gz``. The 
``brie_count.h5ad`` contains all information for downstream analysis, e.g., for
`brie-quant`.


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
