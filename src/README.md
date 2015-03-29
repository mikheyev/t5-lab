# Identifying resistant mutants

For this exercise, we'll use the bowtie2 aligner, and the freebaes  variant calling. There are quite a range of tools for both jobs, and it might be interesting to give them a try. We will use VCFTools to filter the results. The short pipeline below can be easily modified for use with other packages.

The raw reads are in compressed [fastq format](http://en.wikipedia.org/wiki/FASTQ_format). Basically it contains the DNA sequence and a measure of how confident the sequencer is in each of the base calls.

## Aligning reads using bowtie2

```
bowtie2-build ../ref/NC_012967.fasta ../ref/NC_012967
```

This creates a lookup index that for the aligner.

Now we can map the reads. I actually don't recommend doing this, since it takes more than half an hour. Instead the results have been pre-computed for you in ```../data/alignments/*bam```. However, if you do want to run your own computations, run  ```./map.sh``` This script maps the reads, adds read groups (a kind of ID), sorts the output by chromosomal position and creates an index. It uses several programs to do this, piping the ouput from one to the other.



Instead, let's take a look at the pre-computed files to see whether the alignments are any good. 

```
ls -lh ../data/reads/
ls -lh ../data/alignments/
```
samtools 
You see that there alignment files corresponding to every input file.

## Variant calling

There we use FreeBayes on a list of bams produced by our previous analysis.

```
freebayes --ploidy 1 --fasta-reference ../ref/NC_012967.fasta --bam ../data/alignments/*bam -v ../data/var/raw.vcf
```

The outcome of this analysis is in VCF format

## Variant filtering

- Our bacteria are likely going to be a bit different from the actual sequenced reference. We can make our lives simpler, by removing sites where all the samples sites are non-reference. So, we keep any sample that has at least 1 non-reference allele, but less than all non-reference alleles.

-  We can remove low quality sites at the same time. Sites quality is coded by a logarithmic [phred quality score](http://en.wikipedia.org/wiki/Phred_quality_score). A score of 20 corresponds to 
Remove low quality sites and anything where all the samples are different from the reference
```
vcftools --vcf raw.vcf --minQ 20 --non-ref-ac 1 --max-non-ref-ac 4  --recode --out filtered
```

## Visualizing results

We can view the results of our analysis in IGV.

## Exercises

