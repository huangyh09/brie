#!/bin/bash
#a demo for running brie

# DATA_DIR=/afs/inf.ed.ac.uk/user/s13/s1333321/research

DATA_DIR=/home/yuanhua/test/brie-examples
BRIE_DIR=.

anno_ref=$DATA_DIR/anno/gencode.vM6.annotation.gtf
fasta_file=$DATA_DIR/anno/GRCm38.p4.genome.fa


### 1. Splicing events generation
# brie-event -a $anno_ref -o $DATA_DIR/anno/AS_events


### 1.1 Splicing events filtering
# python $BRIE_DIR/brie/events/event_filter.py -a $DATA_DIR/anno/AS_events/SE.gff3 --anno_ref $anno_ref -r $fasta_file -o $DATA_DIR/anno/gencode.vM6.SE 


### 2. Extacting sequence feature for BRIE
##download bigWigSummary binary file:
##http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary
##chmod +x bigWigSummary
##export PATH="/home/yuanhua/Tool/ucsc:$PATH"

anno_file=$DATA_DIR/anno/gencode.vM6.SE.gold.gtf
phast_file=$DATA_DIR/anno/mm10.60way.phastCons.bw
feature_file=$DATA_DIR/anno/mouse_factors.h5
# brie-factor -a $anno_file -r $fasta_file -c $phast_file -o $feature_file -p 4

### 3. BRIE quantification
sam1=$DATA_DIR/sam/E7.75_c1.sorted.bam
sam2=$DATA_DIR/sam/E7.75_c2.sorted.bam

# out_file=$DATA_DIR/out/E7.75_c1
# brie -a $anno_file -s $sam1 -f $feature_file -o $out_file -p 15

out_file=$DATA_DIR/out/E7.75_c2
brie -a $anno_file -s $sam2 -f $feature_file -o $out_file -p 15


### 4. BRIE differential splicing
# brie-diff -1 $DATA_DIR/out/E7.75_c1/samples.h5 -2 $DATA_DIR/out/E7.75_c2/samples.h5 -o $DATA_DIR/out/E7.75_c1_c2.diff.tsv
