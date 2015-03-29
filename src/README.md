# Identifying resistant mutants

For this exercise, we'll use the bowtie2 aligner, and the freebaes  variant calling. There are quite a range of tools for both jobs, and it might be interesting to give them a try. We will use VCFTools to filter the results. The short pipeline below can be easily modified for use with other packages.

## Aligning reads using bowtie2

```
bowtie2-build ../ref/NC_012967.fasta ../ref/NC_012967
```

This creates a lookup index that for the aligner.

Now we can map the reads. I actually don't recommend doing this, since it takes more than half an hour. Instead the results have been pre-computed for you in ```../data/alignments/*bam```. However, if you do want to run your own computations, run  ```./map.sh``` This script maps the reads, adds read groups (a kind of ID), sorts the output by chromosomal position and creates an index. It uses several programs to do this, pipint the ouput from one to the other.

Instead, let's take a look at the pre-computed files to see whether the alignments are any good. 

```
ls -lh ../data/reads/
ls -lh ../data/alignments/
```
samtools 
You see that there alignment files corresponding to every input file.

## Variant calling

There we use FreeBayes on a merged bam input created by samtools. Read groups tell freebayes what reads go with what sample.
```
samtools merge - ../data/alignments/*bam | freebayes --ploidy 1 --fasta-reference ../ref/NC_012967.fasta --stdin -v ../data/var/raw.vcf
```

## Variant filtering



## Visualizing results

We can view the results of our analysis in IGV.


samtools mpileup  -ugf ../ref/NC_012967.fasta ../data/alignments/*bam | bcftools view --min-ac 1 > raw.vcf

## Exercises

