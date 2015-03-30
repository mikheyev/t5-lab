# Analysis of an experimental evolution experiment using next-generation sequencing

The steps below should be done before class, since large files need to be downloaded and installed.

## Installation 

### Software installation

#### Virtual machine

- We prepared a Linux virtual machine (VM) with all the software required for this analysis. The VM can be run in [VirtualBox](https://www.virtualbox.org/wiki/Downloads) which is available for pretty much every operating system.  
- Once VirtualBox is running, download and import the [VM](http://ecoevo.unit.oist.jp/class/t5-lab/T5.ova), using 
- You can then simply run the VM and proceed with the rest of the exercise.
     - The account name is *t5* and the password is *t5-lab*
- Make sure that you turn on bidirectional clipboard in VirtualBox (under Devices...)

#### Install your own data
Alternatively, you can install your own software and skip virtual machine. You will need:

- FreeBayes
- VCFTools
- samtools
- bowtie2

Unless you are comfortable with the command line and know how to compile software, installing all of this software can be a bit tricky. However, it is a good learning experience to go through the installation. In bioinformatics you do a lot of this sort of thing. 

### Inside the virtual machine

Or, if you installed the required software, on you computer.

#### Cloning the git repository

- Once the virtual machine is running, open the terminal and run ```git clone https://github.com/mikheyev/t5-lab.git```
- You can then go to the ```t5-lab``` folder created by git.

#### Download data 

The data are too big for git and need to be downloaded separately.

- Download data and references: ```wget http://ecoevo.unit.oist.jp/class/t5-lab/t5-data.tar.gz``` 
	- Decompress the files ```tar -xzf t5-data.tar.gz && rm t5-data.tar.gz```
		- This step takes a while, maybe a few hours.
- Alternatively, if you have already obtained the data in some way, you can create a shared folder in VirtualBox via *Devices... Shared Floder settings*. Call this folder ```data```
	- You can then create a folder in the VM, e.g., ```mkdir ~/share```
	- Then ```sudo mount -t vboxsf data ~/share```
		- This connects the folder on your computer to the VM.
	- You can then go to the ```t5-lab``` folder and run ```tar -xzf ~/share/t5-data.tar.gz``` to decompress the data.
- You should see two new folders called ```data``` and ```ref```
- The actual exercise is found in the [src](./src/) folder, so you can go there now and start the exercise.



