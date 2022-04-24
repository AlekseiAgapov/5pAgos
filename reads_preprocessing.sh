#!/usr/bin/bash

# The script takes -d argument, that specifies the path to the working directory that should contain single fq.gz file with reads.

working_directory="directory_that_clearly_could_not_exist_because_if_it_existed_it_wouldnt_be_named_so_weird"

while getopts d: flag
do
    case "${flag}" in
        d) working_directory=${OPTARG};;
    esac
done

if [ -d $working_directory ]; then
	echo "Working directory exists, starting the calculations."
else 
	echo "Working directory does not exist! Check -d argument."
	exit
fi

cd $working_directory

file_name=trimmed.fastq.gz
if [ -f $file_name ]; then
	echo "Trimmed.fastq file already exists. The script must have already been run. Exiting."
	exit
else 
	echo "Starting the FASTQ file preprocessing."
fi

# Making the fastqc report for the raw data
mkdir raw_fastqc_report
fastqc -o raw_fastqc_report *fastq.gz

# Trimming the adapters and filtering out reads that are less than 14 or more than 24 nt long.
cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 14 -M 24 -o trimmed.fastq.gz *fastq.gz

# Making the fastqc report for the processed data
mkdir proc_fastqc_report
fastqc -o proc_fastqc_report trimmed.fastq.gz

echo 'The FASTQ file is ready for alignment. Check FastQC report.'
