#!/usr/bin/bash

# Please insert the number of threads you would like to use.
num_threads=8

# First, we create reference sequence files for alignment and visualisation. This script will look for a directory named "reference" and two files in it: "genome.fa" and "plasmid.fa". Attention: the name of the plasmid sequence must be ">plasmid" and the phage must be ">phage".
# The script will combine the three fasta files in one and then make a bowtie index for it.

cd ../
mkdir ref_tmp
cp ./reference/* ./ref_tmp/
cd ./ref_tmp/
touch ./ref.fa
cat genome.fa > ./ref.fa
cat plasmid.fa >> ./ref.fa
cat phage.fa >> ./ref.fa

bowtie-build --threads $num_threads --quiet ./ref.fa ./ref

echo 'Created multifasta file with the number of DNA molecules:'
grep '>' ref.fa | wc -l
echo 'and then built a bowtie index for this reference.'
cd ../

# Making the fastqc report for the raw data
mkdir raw_fastqc_report
fastqc -o raw_fastqc_report *fastq.gz

# Trimming the adapters and filtering out reads that are less than 14 or more than 24 nt long.
cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 14 -M 24 -o trimmed.fastq *fastq.gz

# Making the fastqc report for the processed data
mkdir proc_fastqc_report
fastqc -o proc_fastqc_report trimmed.fastq

echo 'The FASTQ file is ready for alignment. Check FastQC report.'

# Creating folders and copying the necessary files
mkdir ./logo
cp ./pipeline/GC_content_phage.R ./logo/
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

mkdir phage_coverage
cp ./pipeline/phage_coverage.R ./phage_coverage/
cp ./pipeline/normalize_multimappers.py ./phage_coverage/
cp ./pipeline/insert_intervals_after_normalization.py ./phage_coverage/


# The following part of the script will launch the series of alignments using different bowtie options. First it alignes all reads to the reference with -k 1 option and filters out unaligned reads.
echo 'Starting the reads alignment to the reference'
cd ./alignment/
bowtie -k 1 -v 0 -p $num_threads ../ref_tmp/ref ../trimmed.fastq -S all_aligned.sam
gzip ../trimmed.fastq
samtools fastq -@ $num_threads -F 4 ./all_aligned.sam > ./aligned.fastq # this creates a FASTQ file that contains only the reads that mapped to the reference
cp aligned.fastq ../logo/
# Now we align only the reads that are mapped uniquely to the reference (-m 1 option). 
bowtie -m 1 -v 0 -p $num_threads ../ref_tmp/ref ./aligned.fastq -S uniq_aligned.sam
grep -v "@" uniq_aligned.sam | cut -f 2,3,4,10 > ../logo/aligned.tsv
samtools fastq -@ $num_threads -F 4 uniq_aligned.sam > selected.fastq # The reads that were marked as aligned are mapped uniquely. This command saves these reads in a separate FASTQ file.
samtools fastq -@ $num_threads -f 4 uniq_aligned.sam > multimappers.fastq # The reads that were not marked as aligned are multimappers. This command saves these reads in a separate FASTQ file.
# And then we take multimappers and align them to the reference with -a option (all alignments will be reported).
bowtie -a --best -strata -v 0 -p $num_threads ../ref_tmp/ref ./multimappers.fastq -S multimappers.sam
samtools view multimappers.sam | grep -P "\tNC_" | cut -f 1 | sort | uniq > multi_genome.txt # this line saves the names of the reads mapped to the genome sequence in a TXT file
samtools view multimappers.sam | grep -P "\tplasmid\t" | cut -f 1 | sort | uniq > multi_plasmid.txt # this line saves the names of the reads mapped to the plasmid sequence in a TXT file
samtools view multimappers.sam | grep -P "\tphage\t" | cut -f 1 | sort | uniq > multi_phage.txt # this line saves the names of the reads mapped to the phage genome sequence in a TXT file
rm multimappers.sam
../pipeline/filter_multimappers_phage.py
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
bedtools genomecov -d -ibam final_sorted.bam | grep "NC_" | awk '{sum += $3} END {print sum/NR}' > average_cov.txt # calculates average coverage depth of the genome
samtools view final_sorted.bam | cut -f 1 | sort | uniq -c | sed 's/^[ ]*//' | sed 's/ /\t/' > counts.tsv # This line makes a TXT file that contains the information about how many sites each read is mapped to.
samtools view final_sorted.bam | cut -f 1 | sort | uniq | wc -l > total_reads_aligned.txt # Creates a file with the number of reads mapped to the reference
cp total_reads_aligned.txt ../coverage_1000/
cp total_reads_aligned.txt ../coverage_10000/
cp total_reads_aligned.txt ../chi_metaplot/
cp total_reads_aligned.txt ../phage_coverage
samtools view final_sorted.bam | grep -P  '\tNC_' | cut -f 1 | sort | uniq | wc -l > ./aligned_on_genome.txt # Creates a file with the number of reads mapped to the genome
samtools view final_sorted.bam | grep -P  '\tplasmid\t' | cut -f 1 | sort | uniq | wc -l > ./aligned_on_plasmid.txt # Creates a file with the number of reads mapped to the plasmid
samtools view final_sorted.bam | grep -P  '\tphage\t' | cut -f 1 | sort | uniq | wc -l > ./aligned_on_phage.txt # Creates a file with the number of reads mapped to the phage genome
# Making BAM files for plus and minus strands
samtools view -F 16 -b final_sorted.bam > plus.bam
samtools view -f 16 -b final_sorted.bam > minus.bam
python -W ignore ../pipeline/intervals_chi_logo_statistics_phage.py  # Prepares tables for read length distribution and logo, GC-content around mapped reads, chi-metaplot and calculates alignment statistics.
echo 'Alignment statistics is calculated.'


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


echo 'Calculating phage genome coverage.'
cd ../phage_coverage
# Intersect all aligned reads with phage genome intervals. 
bedtools intersect -a intervals.bed -b ../alignment/final_sorted.bam -wa -wb -F 0.5 > intersected.tsv
./normalize_multimappers.py > normalized.tsv
./insert_intervals_after_normalization.py > coverage.tsv
# Intersect reads mapped to the genome in positive orientation with phage genome intervals. 
bedtools intersect -a intervals.bed -b ../alignment/plus.bam -wa -wb -F 0.5 > intersected.tsv
./normalize_multimappers.py > normalized.tsv
./insert_intervals_after_normalization.py > plus_coverage.tsv
# Intersect reads mapped to the genome in negative orientation with phage genome intervals. 
bedtools intersect -a intervals.bed -b ../alignment/minus.bam -wa -wb -F 0.5 > intersected.tsv
./normalize_multimappers.py > normalized.tsv
./insert_intervals_after_normalization.py > minus_coverage.tsv

echo 'Done!'
