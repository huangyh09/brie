#!/bin/bash

# Top Note, besides the demo below, from v2.2 we now include test sets for both
# smart-seq2 and 10x Genomics data, in the brie-tutorials repo:
# https://github.com/huangyh09/brie-tutorials/tree/main/tests


#a demo for running brie2

## Assume you already have exon-skipping events and sequence features
## Otherwise, see Part II for how to get it or downlowd from here:
## https://sourceforge.net/projects/brie-rna/files/annotation/

SAM_DIR=$HOME/research/brie2/germ/sam
OUT_DIR=$HOME/research/brie2/test

mkdir $OUT_DIR
cd $OUT_DIR

### Download annotationa
wget http://ufpr.dl.sourceforge.net/project/brie-rna/annotation/mouse/gencode.vM17/SE.lenient.gff3.gz


ANNO=$OUT_DIR/SE.lenient.gff3.gz

### 0. make cell list table
rm $OUT_DIR/cell_table.tsv
for SAM in $SAM_DIR/*sorted.bam
do
    cellID=`basename $SAM | awk -F".sorted.bam" '{print $1}'`
    echo -e "$SAM\t$cellID" >> $OUT_DIR/cell_table.tsv
done

### 1. BRIE count
samList=$OUT_DIR/cell_table.tsv
brie-count -a $ANNO -S $samList -o $OUT_DIR -p 15


### 2. BRIE quant
brie-quant -i $OUT_DIR/brie_count.h5ad -o $OUT_DIR/brie_quant_gene.h5ad --interceptMode gene
