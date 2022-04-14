#!/usr/bin/python3
# This script takes three TXT files with read names and filters out the reads present in more than one file. 
# Then it takes FASTQ file and returnes FASTQ only with reads present in only one file.
# Input - three TXT files with multimappers from cell genome, phage genome and plasmid, one FASTQ file with multimappers.
# Output - FASTQ file with reads that are present only in one of the three TXT files.

import pandas as pd

# read the file with multimappers mapped to genome
genome = pd.read_csv('../alignment/multi_genome.txt', sep='\t', header=None)
genome.columns = ['read_name']
genome['mapped_to_genome'] = 1

# read the file with multimappers mapped to plasmid
plasmid = pd.read_csv('../alignment/multi_plasmid.txt', sep='\t', header=None)
plasmid.columns = ['read_name']
plasmid['mapped_to_plasmid'] = 1

# read the file with multimappers mapped to phage
phage = pd.read_csv('../alignment/multi_phage.txt', sep='\t', header=None)
phage.columns = ['read_name']
phage['mapped_to_phage'] = 1

# merge three dataframes
merged_df = genome.merge(plasmid, how='outer', on='read_name').merge(phage, how='outer', on='read_name')
merged_df = merged_df.fillna(0)

# select reads that mapped to genome only or plasmid only
selected_reads = merged_df.query('(mapped_to_genome + mapped_to_plasmid + mapped_to_phage) == 1')['read_name']
selected_reads = '@' + selected_reads

# read FASTQ file with all multimappers
multi_fastq = pd.read_csv('../alignment/multimappers.fastq', header=None)
multi_fastq.columns = ['column']
index_selected_names = list(multi_fastq.query('column in @selected_reads').index)

# save the reads that are not cross-chromomsomal multimappers in a new file
index_selected_fastq = []
for element in index_selected_names:
    index_selected_fastq.append(element)
    index_selected_fastq.append(element + 1)
    index_selected_fastq.append(element + 2)
    index_selected_fastq.append(element + 3)

selected_fastq = multi_fastq.query('index in @index_selected_fastq')
selected_fastq.to_csv('../alignment/uniq_multi.fastq', header=False, index=False, sep='\t')


