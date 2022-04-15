# 5pAgos

## What is this project for?
This is the code that was used for NGS analysis published in paper "Prokaryotic Argonaute nucleases with different patterns of bacterial and bacteriophage DNA targeting" 
We sequenced small DNA molecules copurified with prokaryotic Argonaute proteins from *E. coli* cells.
This script aligns the reads to the reference and prepares data for visualization.

## How to use it?
This repo contains 3 directories. <no_phage> and <with_phage> contain pipelines for analysis. Bash script in <pipeline> directory launches the code. The choice depends on whether two (genome and plasmid) or three (genome, plasmid and phage) DNA molecules were present in the cell. The script has two arguments:
  **-p** is the number of threads to use for calculations (number of cores on the machine by default);
  **-d** is the path to the working directory, that should include single raw reads <file_name>.fq.gz and a directory that contains two or three fasta files: genome.fa (the header should be >genome), plasmid.fa (the header should be >plasmid), and phage.fa (the header should be >phage) in case of <with_phage> pipeline.
The script will calculate alignments and create a set of directories in the working directory, where R files for plot drawing will be copied. Manual adjusting of parameters (such as axis limits and labels) in R scripts is required.
The result will be following plots:
- chromosome and phage genome coverage with small DNA reads
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
- samtools https://github.com/samtools/samtools
- bedtools https://github.com/arq5x/bedtools2
- ggplot2 library for R https://github.com/tidyverse/ggplot2
- ggseqlogo library for R https://github.com/omarwagih/ggseqlogo
