#!/bin/bash

mkdir -p ../data/alignments   # create output directory, if it doesn't exist

for i in ../data/reads/*.gz; do
	base=$(basename $i .fastq.gz) # create a sample ID name from the input file name
	echo Mapping $base
	time bowtie2 -p 4 -x ../ref/NC_012967 -U $i `# run bowtie2 with 4 threads mapping reads to reference` \
		--rg-id $base --rg ID:$base --rg LB:Nextera --rg SM:$base --rg PL:ILLUMINA `# give a unique identifier and metadata` | \
		samtools view -Su - `# convert SAM format to BAM` | \
		samtools sort - ../data/alignments/$base  # sort output by chromosome position
	samtools index ../data/alignments/$base.bam   # create a lookup index for rapid file access
done

