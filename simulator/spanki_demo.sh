#!/bin/bash
#a demo for simulation with Spanki.

home_dir=/home/yuanhua/research

hisatDir=$home_dir/tool/hisat-0.1.6-beta
anno_dir=$home_dir/splicing/data/Annotation
hisatRef=$anno_dir/human/hisatRef/GRCh38.p2.genome
junc=$anno_dir/human/junc.v22.tsv


#### Generate reads ####
anno_file=$anno_dir/human/AS_event/SE.filtered.gtf
ref_file=$anno_dir/human/GRCh38.p2.genome.spanki.fa

dat_dir=$home_dir/splicing/data/brie/spanki


# Note, you need to generate the rpk_file yourself for your simulation.
# Here, simuPSI.py does the job for the simulation in BRIE.

rep=rep0
for rpk in rpk200 rpk400 rpk50 #rpk25 rpk100 
do
	echo processing $rpk
	rpk_file=$dat_dir/$rpk/$rpk.$rep.txt
	spankisim_transcripts -o $dat_dir/$rpk -g $anno_file -f $ref_file -bp 76 -frag 200 -ends 2 -t $rpk_file

	mv $dat_dir/$rpk/sim_1.fastq $dat_dir/$rpk/$rpk."$rep"_1.fq
	mv $dat_dir/$rpk/sim_2.fastq $dat_dir/$rpk/$rpk."$rep"_2.fq
	mv $dat_dir/$rpk/transcript_sims.txt $dat_dir/$rpk/$rpk."$rep".sims.txt
	gzip $dat_dir/$rpk/$rpk."$rep"_1.fq $dat_dir/$rpk/$rpk."$rep"_2.fq

	rm -rf $dat_dir/$rpk/tmp $dat_dir/$rpk/log
	rm $dat_dir/$rpk/sim.* $dat_dir/$rpk/junc*


	($hisatDir/hisat -x $hisatRef -1 $dat_dir/$rpk/$rpk."$rep"_1.fq.gz -2 $dat_dir/$rpk/$rpk."$rep"_2.fq.gz --known-splicesite-infile $junc --no-unal -p 20 | samtools view -bS -> $dat_dir/$rpk/$rpk.$rep.bam) 2> $dat_dir/$rpk/$rpk.$rep.err
	samtools sort $dat_dir/$rpk/$rpk.$rep.bam $dat_dir/$rpk/$rpk.$rep.sorted
	samtools index $dat_dir/$rpk/$rpk.$rep.sorted.bam

	rm $dat_dir/$rpk/$rpk.$rep.bam

	($hisatDir/hisat -x $hisatRef -U $dat_dir/$rpk/$rpk."$rep"_1.fq.gz,$dat_dir/$rpk/$rpk."$rep"_2.fq.gz --known-splicesite-infile $junc --no-unal -p 20 | samtools view -bS -> $dat_dir/$rpk/$rpk.$rep.U.bam) 2> $dat_dir/$rpk/$rpk.$rep.U.err
	samtools sort $dat_dir/$rpk/$rpk.$rep.U.bam $dat_dir/$rpk/$rpk.$rep.U.sorted
	samtools index $dat_dir/$rpk/$rpk.$rep.U.sorted.bam

	rm $dat_dir/$rpk/$rpk.$rep.U.bam
done

