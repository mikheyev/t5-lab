# Analysis of an experimental evolution experiment using next-generation sequencing

The steps below should be done before class, since large files need to be downloaded and installed.

## Installation 

### Software installation

#### Virtual machine

- We prepared a Linux virtual machine (VM) with all the software and required for this analysis. The VM can be run in [VirtualBox](https://www.virtualbox.org/wiki/Downloads) which is available for pretty much every operating system.  
- Once VirtualBox is running, download and import the [VM](http://ecoevo.unit.oist.jp/class/t5-lab/T5.ova), using 
- You can then simply run the VM and proceed with the rest of the exercise.
     - The account name is *t5* and the password is *t5-lab*
- Make sure that you turn on bidirectional clipboard in VirtualBox (under Devices...)

#### Install your own software and data
Alternatively, you can install your own software and skip virtual machine. You will need:

- FreeBayes
- VCFTools
- samtools
- bowtie2
- IGV

Unless you are comfortable with the command line and know how to compile software, installing all of this software can be a bit tricky. However, it is a good learning experience to go through the installation. In bioinformatics you do a lot of this sort of thing. 

You can then install this repository, download the data there and unpack them:

	# get this repository
	git clone https://github.com/mikheyev/t5-lab.git   
	cd t5-lab 
	# get data (very slow!)
	wget http://ecoevo.unit.oist.jp/class/t5-lab/t5-data.tar
	# unpack data and remove archive
	tar -xf t5-data.tar && rm t5-data.tar

### Inside the virtual machine


Once the virtual machine is running, open the terminal and update the git repository, which will sync any changes (if needed)

	cd t5-lab
	git pull


### Making sure everything is there

The final directory structure should can be checked by running ```tree``` and should look like this:

```
.
├── README.md
├── data
│   ├── alignments
│   │   ├── mutant1_OIST-2015-03-28.bam
│   │   ├── mutant1_OIST-2015-03-28.bam.bai
│   │   ├── mutant2_OIST-2015-03-28.bam
│   │   ├── mutant2_OIST-2015-03-28.bam.bai
│   │   ├── mutant3_OIST-2015-03-28.bam
│   │   ├── mutant3_OIST-2015-03-28.bam.bai
│   │   ├── mutant4_OIST-2015-03-28.bam
│   │   ├── mutant4_OIST-2015-03-28.bam.bai
│   │   ├── reference_OIST-2015-03-28.bam
│   │   └── reference_OIST-2015-03-28.bam.bai
│   ├── reads
│   │   ├── mutant1_OIST-2015-03-28.fastq.gz
│   │   ├── mutant2_OIST-2015-03-28.fastq.gz
│   │   ├── mutant3_OIST-2015-03-28.fastq.gz
│   │   ├── mutant4_OIST-2015-03-28.fastq.gz
│   │   └── reference_OIST-2015-03-28.fastq.gz
│   └── var
│       └── raw.vcf
├── ref
│   ├── NC_012967.fasta
│   ├── NC_012967.fasta.fai
│   └── NC_012967.gff
└── src
    ├── README.md
    └── map.sh
```

The actual exercise is found in the [src](./src/) folder, so you can go there now and start the exercise.

