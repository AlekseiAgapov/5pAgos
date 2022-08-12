# 5pAgos

## What is this project for?
This is the code that was used for NGS analysis published in paper "Bacterial Argonaute nucleases with different modes of DNA targeting in vitro and in vivo" 
We sequenced small DNA molecules copurified with prokaryotic Argonaute proteins from *E. coli* cells.
This repository contains pipeline for phage genome assembly and scripts that align NGS reads to the reference DNA and prepares data for visualization.

## How to use it?

### Phage genome assembly
File p1_assembly.txt contains commands for phage genome assembly.

### smDNA analysis
Script reads_preprocessing.sh starts reads quality control with FastQC, then remove adaptors with cutadapt and makes quality control once again. For input one should specify (with `-d` argument) a path to a directory that contains single fq.gz file with reads. Output - directories with FastQC output and trimmed.fastq.gz file with the processed reads. **It is important to check the quality of the processed reads and manually adjust cutadapt arguments if needed.**

This repo contains 2 directories: <no_phage> and <with_phage>. Bash scripts pipeline.sh and pipeline_phage.sh in them launch the code. The choice depends on whether two (genome and plasmid) or three (genome, plasmid and phage) DNA molecules were present in the cell. The script has two arguments:
 
  `-p` is the number of threads to use for calculations (number of cores on the machine by default);
 
  `-d` is the path to the working directory, that should include trimmed.fastq.gz file (result of reads_preproessing.sh script) and a directory that contains two or three fasta files: genome.fa (the header should be >genome), plasmid.fa (the header should be >plasmid), and phage.fa (the header should be >phage) in case of <with_phage> pipeline.
 
The script will calculate alignments and create a set of directories in the working directory, where R files for plot drawing will be copied. Manual adjusting of parameters (such as axis limits and labels) in R scripts is required.
The result will be following plots:
- chromosome, plasmid and phage coverage with small DNA reads
- metaplot of small DNA reads coverage of regions around Chi-sites
- aligned reads logo
- GC-content along the guide length and in surrounding sequences of chromosomal DNA
- aligned reads length distribution

## Requirements
This script utilizes some commonly used programs for data analysis and NGS analysis:
- pandas library for Python https://github.com/pandas-dev/pandas
- FastQC https://github.com/s-andrews/FastQC
- cutadapt https://github.com/marcelm/cutadapt
- bowtie https://github.com/BenLangmead/bowtie
- bowtie2 https://github.com/mathworks/bowtie2
- samtools https://github.com/samtools/samtools
- bedtools https://github.com/arq5x/bedtools2
- ggplot2 library for R https://github.com/tidyverse/ggplot2
- ggseqlogo library for R https://github.com/omarwagih/ggseqlogo
- SPAdes https://github.com/ablab/spades

 ## Where to find the data that was processed with this code?
 The raw sequencing reads are deposited in SRA.
 - Small DNA reads are in the BioProject PRJNA827032.
 - Reads for P1 phage genome assembly are in the BioProject PRJNA827167.
