#! /bin/bash
# Tips: test the scripts block by block

SCR_DIR=`pwd`
DAT_DIR=~/splicing/germ
cd $DAT_DIR


#### Download ####
mkdir $DAT_DIR/fastq

i=1
while IFS=$'\t' read -r -a myArray
do
    test $i -eq 1 && ((i=i+1)) && continue
    echo "${myArray[1]}" "${myArray[2]}"
    wget "${myArray[2]}" -O $DAT_DIR/fastq/"${myArray[1]}".fq.gz
done < $SCR_DIR/E-MTAB-4079.130cells.txt


#### Alignment ####
ANNO_DIR=~/annotation
fastaRef=$ANNO_DIR/mouse/GRCm38.p6.genome.fa
ERCCref=$ANNO_DIR/ERCC92.fa

hisatDir=~/tool/hisat2-2.1.0
hisatRef=$ANNO_DIR/mouse/hisatRef/GRCm38.p6.ERCC92

build the reference
mkdir $hisatRef
$hisatDir/hisat2-build $fastaRef,$ERCCref $hisatRef

mkdir $DAT_DIR/sam $DAT_DIR/sam/_errs $DAT_DIR/sam/_logs

i=1
while IFS=$'\t' read -r -a myArray
do
    test $i -eq 1 && ((i=i+1)) && continue
    cell="${myArray[1]}"
    echo $cell
    
    fastq=$DAT_DIR/fastq/$cell.fq.gz
    ($hisatDir/hisat2 -x $hisatRef -U $fastq --no-unal -p 20 | samtools view -bS -> $DAT_DIR/sam/$cell.bam) 2> $DAT_DIR/sam/_errs/$cell.err
    samtools sort $DAT_DIR/sam/$cell.bam -o $DAT_DIR/sam/$cell.sorted.bam
    samtools index $DAT_DIR/sam/$cell.sorted.bam
done < $SCR_DIR/E-MTAB-4079.130cells.txt



#### BRIE quantification ####
SE_anno=$ANNO_DIR/mouse/AS_events/SE.filtered.gtf
FACTOR=$ANNO_DIR/mouse/mouse_features.csv.gz
i=1
while IFS=$'\t' read -r -a myArray
do
    test $i -eq 1 && ((i=i+1)) && continue
    cell="${myArray[1]}"
    echo $cell
    
    sam=$DAT_DIR/sam/$cell.sorted.bam
    brie -a $SE_anno -s $sam -f $FACTOR -o $DAT_DIR/brie/$cell -p 25
done < $SCR_DIR/E-MTAB-4079.130cells.txt



#### BRIE differential splicing ####
SE_anno=$ANNO_DIR/mouse/AS_events/SE.filtered.gtf.gz
FACTOR=$ANNO_DIR/mouse/mouse_features.csv.gz

fileList=$DAT_DIR/brie/E6.5_c21/samples.csv.gz
for i in `seq 22 100`
do
    fileList=$fileList,$DAT_DIR/brie/E6.5_c$i/samples.csv.gz
done
for i in `seq 51 100`
do
    fileList=$fileList,$DAT_DIR/brie/E7.75_c$i/samples.csv.gz
done
