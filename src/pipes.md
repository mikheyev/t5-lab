# Piping in Linux

```map.sh``` uses linux pipes to pass the output of the mapper (bowtie2), which produces a file in [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf), to samtools, which converts the file to BAM format (a binary version of SAM), which in turn passes the output to another samtools module that sorts the reads by coordinates. These steps could easily be accomplished by creating a series of intermediate files, but pipes make for cleaner transitions. You also don't need to worry about the proliferation of intermediate files.

Pipes also allow you to redirect the output from the screen (known as *standard output*) to a file. For example:

	echo Hello, world!
	echo Hello, world! > hello.txt

The second command redirects the "Hello, world!" string to a file called ```hello.txt```. You can read its contents by ```cat hello.txt```

You can append text to files using pipes, too. *E.g.*, ```echo Hello, everyone >> hello.txt``` will add this line to the existing file. Using ```>``` would overwrite the contents of ```hello.txt```

```map.sh``` uses ```|``` to pass output from one program as input to another.  ```cat hello.txt | grep everyone``` prints the contents of ```hello.txt``` to a search utility ```grep``` that looks for a line containing the string ```everyone``` and prints it, if it is found.

## Random useful pipe examples
- You can read the contents of a long file one page at a time by passing them to ```less```
- You can read the end of the file using ```tail```, and an arbitrary line range in a file by using ```head``` and ```tail``` together
- You can select specific columns of an output: ```echo Hello, world! | cut -d " " -f2```
- **Advanced:** Linux and Mac OS come pre-installed with extremely powerful stream editors called ```awk``` and ```sed```. They can be used to perform complicated computations on streams. *E.g.*, you can quickly parse a stream and output only even numbers: ```seq 1 10 | awk '$1 % 2'``` or conduct complex substitutions ```echo Hello, world! | sed 's/world/Gary/'``` Mastering the full power of these editors takes time, but is an essential skill for bioinformatics.
