#! /bin/bash

ANNO_DIR=~/annotation
cd $ANNO_DIR

wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip -P $ANNO_DIR
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/GRCm38.p6.genome.fa.gz -P $ANNO_DIR/mouse
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/gencode.vM17.annotation.gff3.gz -P $ANNO_DIR/mouse
gzip -d $ANNO_DIR/mouse


### downloaded from GENCODE on 8th Feb 2018
fastaRef=$ANNO_DIR/mouse/GRCm38.p6.genome.fa
ERCCref=$ANNO_DIR/ERCC92.fa

hisatDir=/nfs/software/stegle/users/huangh/tool/hisat2-2.1.0
hisatRef=$ANNO_DIR/mouse/hisatRef/GRCm38.p6.ERCC92

### build the reference
mkdir $hisatRef
$hisatDir/hisat2-build $fastaRef,$ERCCref $hisatRef


