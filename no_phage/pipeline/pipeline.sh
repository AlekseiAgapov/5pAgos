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
bowtie-build --threads $num_threads --quiet ./genome.fa ./genome

echo "Created multifasta file with $(grep '>' ref.fa | wc -l) DNA molecules and then built a bowtie index for this reference."
cd ../

# Creating folders and copying the necessary files
mkdir ./logo
cp ./pipeline/GC_content.R ./logo/
cp ./pipeline/ggseqlogo.R ./logo/
cp ./pipeline/lengths.R ./logo/

mkdir alignment

mkdir chi_metaplot
cp ./pipeline/chi_metaplot.R ./chi_metaplot/
cp ./pipeline/normalize_multimappers.py ./chi_metaplot/

mkdir coverage_1000
cp ./pipeline/coverage_1000.R ./coverage_1000/
cp ./pipeline/normalize_multimappers.py ./coverage_1000/
cp ./pipeline/insert_intervals_after_normalization.py ./coverage_1000/

mkdir coverage_10000
cp ./pipeline/coverage_10000.R ./coverage_10000/
cp ./pipeline/normalize_multimappers.py ./coverage_10000/
cp ./pipeline/insert_intervals_after_normalization.py ./coverage_10000/

mkdir plasmid
cp ./pipeline/coverage_plasmid.R ./plasmid/


# Align the reads only to genome in order to calculate GC-content.
echo 'Aligning the reads to the genome only'
cd ./alignment/
bowtie -k 1 -v 0 -p $num_threads ../ref_tmp/genome ../trimmed.fastq.gz --al aligned_genome.fastq -S aligned.sam # this creates a FASTQ file that contains only the reads that mapped to the reference
rm aligned.sam
# The following part of the script will launch the series of alignments using different bowtie options. First it alignes all reads to the reference with -k 10 option and filters out unaligned reads.
echo 'Starting the reads alignment to the reference'
bowtie -k 10 -v 0 -p $num_threads ../ref_tmp/ref ../trimmed.fastq.gz --al aligned.fastq -S aligned.sam # this creates a FASTQ file that contains only the reads that mapped to the reference
samtools view -b -h aligned.sam | samtools sort -@ $num_threads > all_aligned.bam
rm aligned.sam
samtools view -F 4 all_aligned.bam | cut -f 1 | sort | uniq | wc -l > ../plasmid/total_reads_aligned.txt
samtools view -F 16 -b all_aligned.bam > ../plasmid/all_aligned_plus.bam
samtools view -f 16 -b all_aligned.bam > ../plasmid/all_aligned_minus.bam
rm all_aligned*
cp aligned.fastq ../logo/
# Now we align only the reads that are mapped uniquely to the reference (-m 1 option). 
bowtie -m 1 -v 0 -p $num_threads ../ref_tmp/ref ./aligned.fastq -S uniq_aligned.sam
grep -v "@" uniq_aligned.sam | cut -f 2,3,4,10 > ../logo/aligned.tsv
samtools fastq -@ $num_threads -F 4 uniq_aligned.sam > selected.fastq # The reads that were marked as aligned are mapped uniquely. This command saves these reads in a separate FASTQ file.
samtools fastq -@ $num_threads -f 4 uniq_aligned.sam > multimappers.fastq # The reads that were not marked as aligned are multimappers. This command saves these reads in a separate FASTQ file.
# And then we take multimappers and align them to the reference with -a option (all alignments will be reported).
bowtie -a --best -strata -v 0 -p $num_threads ../ref_tmp/ref ./multimappers.fastq -S multimappers.sam
samtools view multimappers.sam | grep -P "\tgenome\t" | cut -f 1 | sort | uniq > multi_genome.txt # this line saves the names of the reads mapped to the genome sequence in a TXT file
samtools view multimappers.sam | grep -P "\tplasmid\t" | cut -f 1 | sort | uniq > multi_plasmid.txt # this line saves the names of the reads mapped to the plasmid sequence in a TXT file
rm multimappers.sam
../pipeline/filter_multimappers.py
echo 'Filtered out the reads that map to more than one chromosome.'

# Aligning the selected reads to the reference
echo 'Making the final BAM file...'
cat uniq_multi.fastq >> selected.fastq
bowtie -a --best -strata -v 0 -p $num_threads ../ref_tmp/ref ./selected.fastq -S final.sam
samtools view -b -h final.sam | samtools sort -@ $num_threads > final_sorted.bam
rm *.sam
echo 'The final_sorted.bam was created'

# Calculating alignment statisitcs and preparing files for further analysis
echo 'Calculating alignment statisitcs'
samtools view final_sorted.bam | cut -f 1 | sort | uniq -c | sed 's/^[ ]*//' | sed 's/ /\t/' > counts.tsv # This line makes a TXT file that contains the information of how many sites each read is mapped to.
samtools view final_sorted.bam | cut -f 1 | sort | uniq | wc -l > total_reads_aligned.txt # Creates a file with the number of reads mapped to the reference
cp total_reads_aligned.txt ../coverage_1000/
cp total_reads_aligned.txt ../coverage_10000/
cp total_reads_aligned.txt ../chi_metaplot/
samtools view final_sorted.bam | grep -P  '\tgenome\t' | cut -f 1 | sort | uniq | wc -l > ./aligned_on_genome.txt # Creates a file with the number of reads mapped to the genome
samtools view final_sorted.bam | grep -P  '\tplasmid\t' | cut -f 1 | sort | uniq | wc -l > ./aligned_on_plasmid.txt # Creates a file with the number of reads mapped to the plasmid
# Making BAM files for plus and minus strands
samtools view -F 16 -b final_sorted.bam > plus.bam
samtools view -f 16 -b final_sorted.bam > minus.bam
python -W ignore ../pipeline/intervals_chi_logo_statistics.py  # Prepares tables for read length distribution and logo, GC-content around mapped reads, chi-metaplot and calculates alignment statistics.
echo 'Alignment statistics is calculated.'
rm *.fastq
rm ../logo/aligned.fastq

echo 'Working on plasmid coverage.'
cd ../plasmid/
bedtools genomecov -d -ibam all_aligned_plus.bam | grep "plasmid" > plus_coverage.tsv
bedtools genomecov -d -ibam all_aligned_minus.bam | grep "plasmid" > minus_coverage.tsv
python -W ignore ../pipeline/plasmid_coverage.py # Makes TSV file with coverage of plus and minus strand in plasmid
rm plus_coverage.tsv
rm minus_coverage.tsv
rm *bam
echo 'Plasmid coverage is calculated.'


echo 'Calculating coverage in 1000-nt intervals.'
cd ../coverage_1000
# Intersect all aligned reads with genome intervals.
bedtools intersect -a intervals.bed -b ../alignment/final_sorted.bam -wa -wb -F 0.5 > intersected.tsv
./normalize_multimappers.py > normalized.tsv
./insert_intervals_after_normalization.py > coverage.tsv
# Intersect reads mapped to the genome in positive orientation with genome intervals. 
bedtools intersect -a intervals.bed -b ../alignment/plus.bam -wa -wb -F 0.5 > intersected.tsv
./normalize_multimappers.py > normalized.tsv
./insert_intervals_after_normalization.py > plus_coverage.tsv
# Intersect reads mapped to the genome in negative orientation with genome intervals. 
bedtools intersect -a intervals.bed -b ../alignment/minus.bam -wa -wb -F 0.5 > intersected.tsv
./normalize_multimappers.py > normalized.tsv
./insert_intervals_after_normalization.py > minus_coverage.tsv
rm normalized.tsv
rm *.py


echo 'Calculating coverage in 10000-nt intervals.'
cd ../coverage_10000
# Intersect all aligned reads with genome intervals. 
bedtools intersect -a intervals.bed -b ../alignment/final_sorted.bam -wa -wb -F 0.5 > intersected.tsv
./normalize_multimappers.py > normalized.tsv
./insert_intervals_after_normalization.py > coverage.tsv
# Intersect reads mapped to the genome in positive orientation with genome intervals. 
bedtools intersect -a intervals.bed -b ../alignment/plus.bam -wa -wb -F 0.5 > intersected.tsv
./normalize_multimappers.py > normalized.tsv
./insert_intervals_after_normalization.py > plus_coverage.tsv
# Intersect reads mapped to the genome in negative orientation with genome intervals. 
bedtools intersect -a intervals.bed -b ../alignment/minus.bam -wa -wb -F 0.5 > intersected.tsv
./normalize_multimappers.py > normalized.tsv
./insert_intervals_after_normalization.py > minus_coverage.tsv
rm normalized.tsv
rm *.py


echo 'Calculating coverage around chi-sites.'
cd ../chi_metaplot
# Intersect reads mapped to the genome in positive orientation with intervals around chi-sites on positive strand.
bedtools intersect -a plus_intervals.bed -b ../alignment/plus.bam -wa -wb -F 0.51 > intersected.tsv
./normalize_multimappers.py > plus_plus.tsv
# Intersect reads mapped to the genome in negative orientation with intervals around chi-sites on positive strand.
bedtools intersect -a plus_intervals.bed -b ../alignment/minus.bam -wa -wb -F 0.51 > intersected.tsv
./normalize_multimappers.py > minus_plus.tsv
# Intersect reads mapped to the genome in negative orientation with intervals around chi-sites on negative strand.
bedtools intersect -a minus_intervals.bed -b ../alignment/minus.bam -wa -wb -F 0.51 > intersected.tsv
./normalize_multimappers.py > minus_minus.tsv
# Intersect reads mapped to the genome in positive orientation with intervals around chi-sites on negative strand.
bedtools intersect -a minus_intervals.bed -b ../alignment/plus.bam -wa -wb -F 0.51 > intersected.tsv
./normalize_multimappers.py > plus_minus.tsv
rm *.py

cd ../
rm -r pipeline
rm -r ref_tmp

echo 'Done!'
