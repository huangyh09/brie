Alternative Splicing Events Generator
=====================================

This ``events_maker`` folder contains scripts to generate alternative splicing events, most of which are modified files from Yarden's package rnaseqlib_. It will probably be published as part of a new version of MISO_.

.. _rnaseqlib: https://github.com/yarden/rnaseqlib
.. _MISO: https://github.com/yarden/MISO

Features
--------

* Generate splicing events from annotation file(s):

  - Skipped exon (SE)

  - Retained intron (RI)

  - Mutually exclusive exon (MXE)

  - Alternative 5' splice site (A5SS) 

  - Alternative 3' splice site (A3SS) 

* Quality check for skipped exon events


Installation
------------

Get ``events_maker`` from the GitHub repository_, and run it locally (see usage). Easier installation will be provided later.

.. _repository: https://github.com/huangyh09/MISO/tree/fastmiso/misopy/events_maker

Note: Some files may requires `gffutils`_ and `pysam`_.

.. _gffutils: https://pythonhosted.org/gffutils/
.. _pysam: http://pysam.readthedocs.io


Usage of ``events_maker``
-------------------------

``events_maker`` can take a transcript-centric annotation to generate a exon-centric alternative events which can be used by MISO for quantitation using RNA-Seq. Here, you could input a single or multiple annotation files, with a same or different formats. The supported formats for the input annotation file include ``uscs`` table, ``gtf`` format and ``gff3`` format, and the output file will be in gff3 formats. 

To run it, you could use one of the following command lines.::

 $ $dir/gff_make_annotation.py my_anno.gtf gtf ./gff --flanking-rule commonshortest --genome-label hg38

 $ $dir/gff_make_annotation.py my_anno.gtf,my_anno2.gtf gtf ./gff --flanking-rule commonshortest --genome-label hg38

 $ $dir/gff_make_annotation.py my_anno.gtf,my_anno2.gff3,my_table.txt gtf,gff3,uscs ./gff --flanking-rule commonshortest --genome-label hg38

This will use the script ``gff_make_annotation.py`` to load the input annotation file(s) in given format(s) and output the set of GFF3 events into the directory ``./gff``. The ``--flanking-shortest`` parameter specifies how to pick the flanking exons of an alternative event when there are multiple options (e.g. which flanking exons to use for an alternatively skipped exon.) In the above call, we chose to the take the common shortest region as flanking exons. The ``--genome-label`` option specifies what the name of the GFF annotation should be. This call should generate an output directory with a GFF3 file for each of the event types: ::

  $ ls ./gff/commonshortest/
  A3SS.hg38.gff3  A5SS.hg38.gff3  MXE.hg38.gff3  RI.hg38.gff3  SE.hg38.gff3


**Example of annotation formats mentioned above:**

* ``ucsc`` (UCSC table): ::

    585     ENST00000619216.1       chr1    -       17368   17436   17368   17368   1       17368,  17436,  0       MIR6859-2       none    none    -1,
    585     ENST00000473358.1       chr1    +       29553   31097   29553   29553   3       29553,30563,30975,      30039,30667,31097,      0       MIR1302-11      none    none    -1,-1,-1,
    585     ENST00000469289.1       chr1    +       30266   31109   30266   30266   2       30266,30975,    30667,31109,    0       MIR1302-11      none    none    -1,-1,

* ``gtf`` (gtf_ format): ::
   
    XV  ensembl gene  93395 94402 . - . gene_id "YOL120C"; gene_name "RPL18A"; gene_biotype "protein_coding";
    XV  ensembl mRNA  93395 94402 . - . gene_id "YOL120C"; gene_name "RPL18A"; gene_biotype "protein_coding";
    XV  ensembl exon  94291 94402 . - . gene_id "YOL120C"; gene_name "RPL18A"; gene_biotype "protein_coding";
    XV  ensembl exon  93395 93843 . - . gene_id "YOL120C"; gene_name "RPL18A"; gene_biotype "protein_coding";
    XI  ensembl gene  431906  432720  . + . gene_id "YKL006W"; gene_name "RPL14A"; gene_biotype "protein_coding";
    XI  ensembl mRNA  431906  432720  . + . gene_id "YKL006W"; gene_name "RPL14A"; gene_biotype "protein_coding";
    XI  ensembl exon  431906  432034  . + . gene_id "YKL006W"; gene_name "RPL14A"; gene_biotype "protein_coding";
    XI  ensembl exon  432433  432720  . + . gene_id "YKL006W"; gene_name "RPL14A"; gene_biotype "protein_coding";

.. _gtf: http://www.ensembl.org/info/website/upload/gff.html


* ``gff3`` (gff3_ format): ::
   
    chr1	Gencode.v23	gene	184923	200322	.	-	0	ID=ENSG00000279457.2
    chr1	Gencode.v23	mRNA	184923	195411	.	-	0	ID=ENST00000623834.2;Parent=ENSG00000279457.2
    chr1	Gencode.v23	exon	184923	185350	.	-	0	ID=ENST00000623834.2.0;Parent=ENST00000623834.2
    chr1	Gencode.v23	exon	185491	185559	.	-	0	ID=ENST00000623834.2.1;Parent=ENST00000623834.2
    chr1	Gencode.v23	exon	187376	187577	.	-	0	ID=ENST00000623834.2.2;Parent=ENST00000623834.2
    chr1	Gencode.v23	exon	187755	187886	.	-	0	ID=ENST00000623834.2.3;Parent=ENST00000623834.2

.. _gff3: http://www.broadinstitute.org/annotation/argo/help/gff3.html


Notes
-----

* RI may need further check as I changed to use the old version.

* Events definition is mainly copied from Yarden's codes and hasn't been checked carefully. I mean the codes rather than the output files.

* Quality check will be provided, for example SE exon surrounded by AG-as_exon-GT.