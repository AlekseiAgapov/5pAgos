#!/usr/bin/bash

# The script takes two arguments:
#-p is the number of threads to use (number of cores on the machine by default);
#-d is the path to the working directory that should contain trimmed.fastq.gz file with processed reads and the ./reference/ directory with FASTA files.

num_threads=$(nproc)
working_directory="directory_that_clearly_could_not_exist_because_if_it_existed_it_wouldnt_be_named_so_weird"

while getopts p:d: flag
do
    case "${flag}" in
        p) num_threads=${OPTARG};;
        d) working_directory=${OPTARG};;
    esac
done

if [ $num_threads -gt 0 ]; then
	echo "$num_threads threads"
else 
	echo "-p argument should be > 0."
	exit
fi


if [ -d $working_directory ]; then
	echo "Working directory exists, starting the calculations."
else 
	echo "Working directory does not exist! Check -d argument."
	exit
fi

cp -r ../pipeline/ $working_directory
cd $working_directory

# Make Python scripts executable.
chmod +x ./pipeline/filter_multimappers.py
chmod +x ./pipeline/insert_intervals_after_normalization.py
chmod +x ./pipeline/intervals_chi_logo_statistics.py
chmod +x ./pipeline/normalize_multimappers.py

# First, the script will create reference sequence files for alignment and visualisation. This script will look for a directory named "reference" and two files in it: "genome.fa" and "plasmid.fa". Attention: the name of the genome sequence must be ">genome" and the plasmid sequence - ">plasmid"!
# The script will combine the two fasta files in one and then make a bowtie index for it.

mkdir ref_tmp
cp ./reference/* ./ref_tmp/
cd ./ref_tmp/
touch ./ref.fa
cat ./genome.fa >> ./ref.fa
cat ./plasmid.fa >> ./ref.fa

bowtie-build --threads $num_threads --quiet ./ref.fa ./ref

echo "Created multifasta file with $(grep '>' ref.fa | wc -l) DNA molecules and then built a bowtie index for this reference."
cd ../alignment



# The following part of the script will launch the series of alignments using different bowtie options. First it alignes all reads to the reference with -k 10 option and filters out unaligned reads.
echo 'Starting the reads alignment to the reference'
bowtie -k 10 -v 0 -p $num_threads ../ref_tmp/ref ../trimmed.fastq.gz -S aligned.sam
samtools view -b -h aligned.sam | samtools sort -@ $num_threads > all_aligned.bam
rm aligned.sam
samtools view -F 4 all_aligned.bam | cut -f 1 | sort | uniq | wc -l > ../plasmid/total_reads_aligned.txt
rm all_aligned.bam
cd ../
rm -r pipeline/
rm -r ref_tmp/

echo 'Done!'
