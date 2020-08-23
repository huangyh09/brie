#!/bin/bash

BAM_DIR=$1
suffix=".sorted.bam"

samList=$BAM_DIR/../cell_bams.tsv

## Make cell list
echo "" > $samList
for sample in `ls $BAM_DIR/*$suffix`
do
    NAME=$(basename $sample $suffix)
    echo $sample$'\t'$NAME >> $samList
done
