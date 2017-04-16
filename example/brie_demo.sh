#!/bin/bash
#a demo for running brie

################################################################
#                                                              #
#     Part I. Isoform estimate and differential splicing       #
#                                                              #
################################################################

## Assume you already have exon-skipping events and sequence features
## Otherwise, see Part II for how to get it or downlowd from here:
## https://sourceforge.net/projects/brie-rna/files/annotation/

DATA_DIR=$PWD
DATA_DIR=/disk/scratch/yhuang/brei-examples

anno_file=$DATA_DIR/anno/gencode.vM12.SE.gold.gtf
feature_file=$DATA_DIR/anno/mouse_factors.SE.gold.csv.gz

### 1.1 BRIE quantification
sam1=$DATA_DIR/sam/E7.75_c1.sorted.bam
sam2=$DATA_DIR/sam/E7.75_c2.sorted.bam

out_file=$DATA_DIR/out/E7.75_c1
brie -a $anno_file -s $sam1 -f $feature_file -o $out_file -p 15

out_file=$DATA_DIR/out/E7.75_c2
brie -a $anno_file -s $sam2 -f $feature_file -o $out_file -p 15


# ### 1.2 BRIE differential splicing
fileList=$DATA_DIR/out/E7.75_c1/samples.csv.gz,$DATA_DIR/out/E7.75_c2/samples.csv.gz
brie-diff -i $fileList -o $DATA_DIR/out/E7.75_c1_c2.diff.tsv
