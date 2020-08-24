#!/bin/bash
#a demo for generate splicing events and sequence features

################################################################
#                                                              #
#     Part II. Define splicing events and fetch features       #
#                                                              #
################################################################


### Top notice: ###
# This requires briekit package: https://pypi.org/project/briekit/
# It requires to install in Python 2 environment
# Readmore: https://github.com/huangyh09/briekit/wiki

echo "Please use briekit package in Python 2 environment!"


DATA_DIR=/disk/scratch/yhuang/anno

### 1. download files
# wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/GRCm38.p5.genome.fa.gz -O $DATA_DIR/GRCm38.p5.genome.fa.gz
# wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gtf.gz -O $DATA_DIR/gencode.vM12.annotation.gtf.gz
# wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.p7.genome.fa.gz -O $DATA_DIR/GRCh38.p7.genome.fa.gz
# wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz -O $DATA_DIR/gencode.v25.annotation.gtf.gz
# gzip -d $DATA_DIR/*
# wget ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw -O $DATA_DIR/mm10.60way.phastCons.bw
# wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw -O $DATA_DIR/hg38.phastCons100way.bw

anno_ref1=$DATA_DIR/gencode.vM12.annotation.gtf
fasta1=$DATA_DIR/GRCm38.p5.genome.fa
anno_ref2=$DATA_DIR/gencode.v25.annotation.gtf
fasta2=$DATA_DIR/GRCh38.p7.genome.fa



### 2. Splicing events filtering
## https://brie-rna.sourceforge.io/manual.html#splicing-events

# briekit-event -a $anno_ref1 -o $DATA_DIR/mouse_AS
# briekit-event-filter -a $DATA_DIR/mouse_AS/SE.gff3 --anno_ref $anno_ref1 -r $fasta1
# briekit-event-filter -a $DATA_DIR/mouse_AS/SE.gff3 --anno_ref $anno_ref1 -r $fasta1 -o $DATA_DIR/mouse_AS/SE.extended --add_chrom chrX,chrY --as_exon_min 10 --as_exon_max 100000000 --as_exon_tss 10 --as_exon_tts 10 --no_splice_site #--keep_overlap

# briekit-event -a $anno_ref2 -o $DATA_DIR/human_AS
# briekit-event-filter -a $DATA_DIR/human_AS/SE.gff3 --anno_ref $anno_ref2 -r $fasta2
# briekit-event-filter -a $DATA_DIR/human_AS/SE.gff3 --anno_ref $anno_ref2 -r $fasta2 -o $DATA_DIR/human_AS/SE.extended --add_chrom chrX,chrY --as_exon_min 10 --as_exon_max 100000000 --as_exon_tss 10 --as_exon_tts 10 --no_splice_site #--keep_overlap


### 3. Extacting sequence feature for BRIE
## https://brie-rna.sourceforge.io/manual.html#sequence-features

anno_file=$DATA_DIR/anno/gencode.vM6.SE.gold.gtf

phast1=$DATA_DIR/mm10.60way.phastCons.bw
feature1=$DATA_DIR/mouse_AS/mouse_factors.SE.most.csv
feature2=$DATA_DIR/mouse_AS/mouse_factors.SE.gold.csv

briekit-factor -a $DATA_DIR/mouse_AS/SE.extended.gold.gtf -r $fasta1 -c $phast1 -o $feature1 -p 20
briekit-factor -a $DATA_DIR/mouse_AS/SE.gold.gtf -r $fasta1 -c $phast1 -o $feature2 -p 1


phast2=$DATA_DIR/hg38.phastCons100way.bw
feature3=$DATA_DIR/human_AS/human_factors.SE.most.csv
feature4=$DATA_DIR/human_AS/human_factors.SE.gold.csv

briekit-factor -a $DATA_DIR/human_AS/SE.extended.gold.gtf -r $fasta2 -c $phast2 -o $feature3 -p 20
briekit-factor -a $DATA_DIR/human_AS/SE.gold.gtf -r $fasta2 -c $phast2 -o $feature4 -p 20
