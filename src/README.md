# Identifying resistant mutants

For this exercise, we'll use the bowtie2 aligner, and the FreeBayes variant calling package. We will use VCFTools to filter the results. The short pipeline below can be easily modified for use with other packages, of which ther are quite a range, and you are welcome to try some other ones for homework.

The raw reads are in compressed [fastq format](http://en.wikipedia.org/wiki/FASTQ_format). Basically it contains the DNA sequence and a measure of how confident the sequencer is in each of the base calls.

## A few quick notes about relative paths paths

- For the commands below you should be in the ```~/t5-data/src``` folder. 

- The ```..``` is shorthand for one directory above the current one. 

- Similarly the relative path ```.``` means in the current directory. On linux, when running software in the current directory you have to say it explicitly e.g., ```./foo```

- Relative paths allow you to write code that accesses files relative to the current directory. This is handy if your package, such as the t5-lab repository is placed in a different location other than the home directory (otherwise know as ```~``` or ```$HOME```).

- In the virtual machine your current path is displayed at the prompt after the colon (:)

## Aligning reads using bowtie2

We take data from [Jeong *et al* 2009](http://www.ncbi.nlm.nih.gov/pubmed/19786035), who sequenced a couple of *E. coli* B strains.

```
bowtie2-build ../ref/NC_012967.fasta ../ref/NC_012967
```

This creates a lookup index that for the aligner.

Now we can map the reads. I actually don't recommend doing this, since it takes more than half an hour. Instead the results have been pre-computed for you in ```../data/alignments/*bam```. However, first, look inside ```map.sh``` and figure out what's going on.

```map.sh``` makes extensive use of linux pipes, a way to pass seamlessly pass data between different programs. You should read about pipes [here](./pipes.md), before tackling the homework assignments.

However, if you do want to run your own computations, run  ```./map.sh``` This script maps the reads, adds read groups (sample ID, and meta-data, such as sequencing platform, *etc.*), sorts the output by chromosomal position and creates an index. It uses several programs to do this, piping the output from one to the other.

Instead, let's take a look at the pre-computed files to see whether the alignments are any good. 

```
ls -lh ../data/reads/
ls -lh ../data/alignments/
```
You see that there alignment files corresponding to every input file, each being roughly the same size as the input file.

You can also check to see the mapping percentages using ```samtools flagstat [filename]```

## Variant calling

There we use FreeBayes on a list of bams produced by our previous analysis. I also don't recommend running this command, as it takes a while to run. The output file ```../data/var/raw.vcf``` has been pre-computed.

```
freebayes --ploidy 1 --fasta-reference ../ref/NC_012967.fasta --bam ../data/alignments/*bam -v ../data/var/raw.vcf
```

The outcome of this analysis is in [VCF format](http://samtools.github.io/hts-specs/VCFv4.2.pdf). It is designed to be human-readable, but just barely. At its core, it contains evidence for variability at sites throughout the genome (column 6), and what the genotype of each sample is at that site (columns 9 and onward). See sections 1.3 and 1.4 of the VCF format specifications for more information.

## Variant filtering

Our bacteria are likely going to be a bit different from the actual sequenced reference. We can make our lives simpler, by removing sites where all the samples sites are non-reference. So, we keep any sample that has at least 1 non-reference allele, but less than all non-reference alleles.

We can remove low quality sites at the same time. Sites quality is coded by a logarithmic [phred quality score](http://en.wikipedia.org/wiki/Phred_quality_score). A score of 20 corresponds to 
Remove low quality sites and anything where all the samples are different from the reference
```
vcftools --vcf ../data/var/raw.vcf --minQ 20 --non-ref-ac 1 --max-non-ref-ac 4  --recode --out ../data/var/filtered
```

## Visualizing results in IGV

### Loading all of the data

- We can view the results of our analysis in IGV by running ```igv``` from the command line
- Use *Genomes... Load genome from file* to load ```ref/NC_012967.fasta```
- The load the annotation usnig *File... Load from file*  ```ref/NC_012967.gff```
- Then load all of the bam files from ```data/alignments``` and the filtered vcf file ```data/var/filtered.recode.vcf```

### Data intepretation

The T5 receptor is called [**fhuA**](https://www.wikigenes.org/e/gene/e/944856.html) and you can type this into the IGV location window to zoom to the region of particular interest.

#### In class exercises

1. Take a look at mutant #3. What is this mutation? What biological effect does it have?
- Take a look at the VCF output of FreeBayes on mutant #1, ignoring alignment data for the moment. How would you interpret this mutation?
- Now, take a look at the raw alignment data for mutant #1. Do you agree with the FreeBayes on the nature of the mutation?
- Characterize all four mutations, and explain how they contribute to T5 phage resistance.
	- Hint: Think about what the read mapper does and how this relates to coverage variation seen in the samples.
	- Hint: there are three different types of mutations.
	- Hint: You may wish to search the reference genome for motifs. You can do it at the UCSC genome browser using a search tool called [BLAT](http://microbes.ucsc.edu/cgi-bin/hgBlat?hgsid=411911&command=start)
- Why are there no mutations registered in the VCF file for mutants #2 and #4? 

#### Programming exercises for next week

1. If we didn't know where to look for mutations, say if we wanted to characterize a novel phage, we would have to look genome-wide. Parse the VCF to find out how many other mutations are specific to just one isolate.
- The transposes can show some sequence specificity. Are there systematic coverage biases? To check this, you can see if there the reads start coordinates on average the same across different samples. Use the output of ```samtools view``` on the bam files to see if the start positions are correlated across samples. Check [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf) for more details on the output. You will be looking at column 4, the start position.
	- Make sure you pipe the output of ```samtools view [filename]``` somewhere, or terminate it in some way. Otherwise, it will print a line for every read in the data set, which will take a while. One way to print just a few lines is by [piping](./pipes.md) the output to a command called ```head```, *e.g.*, ```samtools view foo.bam |head -50``` will display the first 50 lines in the file.
- What is the average depth of coverage per sequenced bacterial genome, i.e., what is that average number of time each base is sequenced?
- Bonus: re-write the analysis pipeline to use another aligner and base caller (*e.g.*, the [Stamy](http://www.well.ox.ac.uk/project-stampy) and [Platypus](http://www.well.ox.ac.uk/platypus) combo). Re-run the data analysis and compare VCF files. Do you get the same result?
