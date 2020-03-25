# Remove Host from Nanopore Reads Using Minimap2 and Samtools 

This tools takes basecalled fastq files from Nanopore sequencing, and removes reads mapping to a reference file from the fastq. The most obvious use is to remove the host contaminants from metagenomic sequencing.

It uses Minimap2 to map to the reference sequence and samtools to pull out the reads that did not map.

This workflow depends on Python 3.6.9 and the Gooey module (for the GUI interface) as well as Minimap2 and Samtools. This workflow includes a Conda recipe to install these environmental dependencies (Minimap2, Samtools, Python 3.6.9, and the Gooey module) in a easy and reproducible fashion using Miniconda/Anaconda.

A script for easy installation of Miniconda is also included.

The GUI version also requires direct access to a linux computer or X11 forward. 

## Summary - Installation 
1. Clone Repository 
2. Install Conda if not already in environment
3. Create conda environment
4. Acquire reference file - download the human_g1k_v37 reference file using the following command
```
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
```

## Summary - How to run after installation.
1. Activate conda environment - `conda activate py36_minimap2_samtools`
2. Run GUI version - `~/removehost/removehost_GUI.py`
3. Or run CLI version - `~/removehost/removehost.py -h`
4. Deactivate conda environment - `conda deactivate`

## Clone this code using GIT

### Install git for Debian systems using the following command (if necessary)
```
sudo apt update
sudo apt install git
```

## Installation directions 
These instructions install the code into your home path. Change the instructions if appropriate. 

### Clone the code from repository
```
cd ~
git clone https://github.com/jiangweiyao/removehost.git
```

### Install Miniconda, Download human_g1k_v37, and build conda environment
You can do all three things by using the following command. Then, you can skip the next 3 sections.
```
. ~/removehost/install_all.sh
```

### Install Miniconda (if no Conda is install on system). 
You can run the prepackaged script install_miniconda.sh to install into your home directory (recommended) by using the following command
```
. ~/removehost/install_miniconda.sh
```

Detailed instruction on the the Miniconda website if anything goes wrong:
https://conda.io/projects/conda/en/latest/user-guide/install/linux.html

### Clone the environment. Need to do once.

We use conda to create an environment (that we can activate and deactivate) to install our dependent software and resolve their dependencies. This environment is called "py36_minimap2_samtools". The following command assumes your environment file is in your home path. Modify as appropriate.

```
conda env create -f ~/removehost/environment.yml
```

The command to generate the environment originally is in the included `py36_minimap2_samtools.txt` file. 

## Download Reference Files
The human_g1k_v37 reference file can be downloaded from the link below:
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
with the following command
```
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
```

You can use any fasta or fasta.gz format reference file for this tool.


## Run the code.

Activating your environment makes the software you installed in that environment available for use. You will see "(py36_minimap2_samtools)" in front bash after activation.
```
conda activate py36_minimap2_samtools
```

Run the GUI version with the following command
```
~/removehost/removehost_GUI.py
```

Run the CLI version with the test files with the following command
```
~/removehost/removehost.py -i ~/removehost/test -r ~/removehost/human_g1k_v37.fasta.gz -o ~/dehost_output/test
```

Get help for the CLI version with the following command
```
~/removehost/removehost.py -h
```

CLI generic Usage
```
~/removehost/removehost.py -i <input directory contain fastq> -r <reference fasta/fasta.gz> -o <output location>
```


When you are finished running the workflow, exit out of your environment by running `conda deactivate`. Deactivating your environment exits out of your current environment and protects it from being modified by other programs. You can build as many environments as you want and enter and exit out of them. Each environment is separate from each other to prevent version or dependency clashes. The author recommands using Conda/Bioconda to manage your dependencies.

## Reference Files
The human_g1k_v37 reference file from NCBI is included with this repository. It was downloaded from: 
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz

You can use any fasta or fasta.gz format reference file for this tool. 

## What the code is actually doing

The code is running the following command using the reference file (human_g1k_v37.fasta.gz) and the fastq file (E0126.fastq).

```
minimap2 -ax map-ont human_g1k_v37.fasta.gz E0126.fastq -t 12 > E0126_g1k.sam
samtools view -u -f 4 E0126_g1k.sam > E0126_g1k_unmapped.sam
samtools bam2fq E0126_g1k_unmapped.sam > E0126_g1k_unmapped.fastq
```
