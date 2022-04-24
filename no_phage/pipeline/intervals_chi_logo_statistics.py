#!/usr/bin/python3
# This script is a part of a pipeline for calculation and visualization of illumina reads parameters and distribution along the reference DNA. It does several things:
# 1. Creates BED files with equal intervals of length 1000 and 10000 nt for the genome.
# 2. Calculates GC-content for each DNA from reference genome.
# 3. Calculates reads length distribution.
# 4. Prepares TXT file for reads logo calculation.
# 5. Calculates GC-content for each position in the window -15 : +35 nucleotides from 5'-end of the mapped read.
# 6. Finds chi-sites in the genome and prepares a BED file with 500-nt intervals around each chi-site.
# 7. Calculates alignment statistics and write it in a TSV file.

import pandas as pd

# read the genome and the plasmid files and count their lengths
genome = pd.read_csv('../reference/genome.fa', sep='\t')
genome.columns = ['fasta']
genome = str(genome['fasta'].sum()).upper()
genome_length = len(genome)

plasmid = pd.read_csv('../reference/plasmid.fa', sep='\t')
plasmid.columns = ['fasta']
plasmid = str(plasmid['fasta'].sum()).upper()
plasmid_length = len(plasmid)

print('The FASTA files have been read.')


def create_bed(interval, dna, dna_name):
  
    '''
    This function creates a BED file with equal intervals.
    You should specify the length of interval (int) and provide two strings with dna sequence and name.
    '''
    dna_length = len(dna)
    end = []
    name = []
    
    coordinate = interval
    interval_name = 0
    while coordinate < dna_length:
        end.append(coordinate)
        name.append(interval_name)
        coordinate += interval
        interval_name += 1

    bed = pd.DataFrame({'end': end, 'name': name})
    bed['start'] = bed['end'] - interval
    bed['chr'] = dna_name
    bed = bed[['chr', 'start', 'end', 'name']]
    bed.loc[len(bed)] = [dna_name, coordinate - interval, dna_length, interval_name]
    file_name = '../coverage_' + str(interval) + '/intervals.bed'

    bed.to_csv(file_name, header=False, index=False, sep='\t')
    

# Apply the function to create two BED files that split genome into equal intervals of desired length.
create_bed(1000, genome, 'genome')
create_bed(10000, genome, 'genome')

print('BED files for whole genome coverage calculation have been created.')

# Calculate GC-content in dna files.

def GC_content(seq):
    '''
    This function calculates and returns percent of "G" or "C" in the input string.
    '''
    GC_number = 0
    for N in seq:
        if N == 'C' or N == 'G':
            GC_number += 1
    return GC_number / len(seq) * 100

# Save plasmid GC-content in a separate TXT file.
with open('../logo/total_plasmid_GC.txt', 'w') as ouf:
    ouf.write(str(GC_content(plasmid)))
    ouf.write('\n')

# Save genomic GC-content in a separate TXT file.
with open('../logo/total_genome_GC.txt', 'w') as ouf:
    ouf.write(str(GC_content(genome)))
    ouf.write('\n')


# Read FASTQ file with the reads aligned to the reference DNA. Save only reads sequences and calculate their length.
reads = pd.read_csv('../logo/aligned.fastq', sep='\t', header=None)
reads.columns = ['fastq']
reads = reads.query('(index+3)%4 == 0')
reads['length'] = reads.apply(lambda x: len(x['fastq']),axis=1)


# Save table with reads length distribution in a separate file.
lengths = reads['length'].value_counts().reset_index().sort_values(by='index').reset_index(drop=True)
lengths.columns = ['length', 'number']
lengths['percent'] = lengths['number'] / (lengths['number'].sum()) * 100
lengths.to_csv('../logo/lengths.tsv', sep='\t', header=True, index=False)

print('Ready to visualize guide length distribution.')


# Select the reads that are 17 nt or longer.
reads_17 = reads.query('length >= 17')['fastq']

# Cut these reads at 17 nt.

reads_17_list = []

for read in reads_17:
    reads_17_list.append(read[:17])
    
reads_17 = pd.DataFrame(reads_17_list)

# Save these sequences in a separate file.
reads_17.to_csv('../logo/logo.txt', index=False, header=False)

print('Ready to calculate guide logo.')

# Read a preformated alignment file
reads = pd.read_csv('../logo/aligned.tsv', sep='\t', header=None)
reads.columns = ['flag', 'chr', 'position', 'sequence']
reads['length'] = reads.apply(lambda x: len(x['sequence']),axis=1)

# Create tables for each chromosome and orientation. Calculate start position (5') for each read.
plus_genome = reads.query('flag == 0 and chr == "genome"')
plus_genome['start_position'] = plus_genome['position']

plus_plasmid = reads.query('flag == 0 and chr == "plasmid"')
plus_plasmid['start_position'] = plus_plasmid['position']

minus_genome = reads.query('flag == 16 and chr == "genome"')
minus_genome['start_position'] = minus_genome['position'] + minus_genome['length'] - 1

minus_plasmid = reads.query('flag == 16 and chr == "plasmid"')
minus_plasmid['start_position'] = minus_plasmid['position'] + minus_plasmid['length'] - 1


# Collect sequences from -15 to + 30 positions from the 5'-end of each read in lists.

plus_genome_logo = []

for start_position in plus_genome['start_position']:
    if start_position >= 16 and start_position <= (genome_length - 30):
        plus_genome_logo.append(genome[start_position - 16 : start_position + 29])

plus_plasmid_logo = []

for start_position in plus_plasmid['start_position']:
    if start_position >= 16 and start_position <= (plasmid_length - 30):
        plus_plasmid_logo.append(plasmid[start_position - 16 : start_position + 29])
        

def reverse_complement(seq):
    '''
    This function takes a string with dna sequence and returns a string with reverse complement sequence.
    '''
    reverse_seq = seq[::-1]
    dictionary = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    reverse_complement_seq = str()
    for letter in reverse_seq:
        reverse_complement_seq += dictionary[letter]
    return reverse_complement_seq
    
minus_genome_logo = []

for start_position in minus_genome['start_position']:
    if start_position >= 30 and start_position <= (genome_length - 16):
        minus_genome_logo.append(reverse_complement(genome[start_position - 30 : start_position + 15]))
        
minus_plasmid_logo = []

for start_position in minus_plasmid['start_position']:
    if start_position >= 30 and start_position <= (plasmid_length - 16):
        minus_plasmid_logo.append(reverse_complement(plasmid[start_position - 30 : start_position + 15]))

# Calculate GC content in each position and write the data in one file.
genome_logo = plus_genome_logo + minus_genome_logo
genome_GC = []

for nucleotide in range(len(genome_logo[0])):
    position = str()
    for read in genome_logo:
        position += read[nucleotide]
    genome_GC.append(GC_content(position))
    
plasmid_logo = plus_plasmid_logo + minus_plasmid_logo
plasmid_GC = []

for nucleotide in range(len(plasmid_logo[0])):
    position = str()
    for read in plasmid_logo:
        position += read[nucleotide]
    plasmid_GC.append(GC_content(position))

positions = []
for i in range(-15, 0):
    positions.append(i)
for i in range(1, 31):
    positions.append(i)

GC = pd.DataFrame({'position': positions, 'genome': genome_GC, 'plasmid': plasmid_GC})
GC.to_csv('../logo/GC_content.tsv', index=False, header=True, sep='\t')

print('Ready to visualize GC-content.')

# Define chi-sequence.
chi_sequence = "GCTGGTGG"

# Calculate all positions of chi-sequences in the positive strand and save them to TXT file.
plus_strand_chi = []

for coordinate in range(len(genome) - 8):
	if genome[coordinate:(coordinate + 8)] == chi_sequence:
		plus_strand_chi.append(coordinate)

plus_df = pd.DataFrame(plus_strand_chi)
plus_df.to_csv('../coverage_1000/plus_chi.txt', header=False, index=False, sep='\t')


# Calculate all positions of chi-sequences in the negative strand and save them to TXT file.
minus_strand_chi = []

for coordinate in range(8, len(genome)):
	if genome[coordinate:(coordinate + 8)] == reverse_complement(chi_sequence):
		minus_strand_chi.append(coordinate)

minus_df = pd.DataFrame(minus_strand_chi)
minus_df.to_csv('../coverage_1000/minus_chi.txt', header=False, index=False, sep='\t')


# Select only those chi-sites that are not in the interval 1200000-1700000 nt.

def site_slice(coordinate_list, start, end):
    '''
    This function takes a list of coordinate and two numbers: start and end.
    It returns another list with those coordinates that fit in the interval [start:end].
    '''
    result_list = []
    for coordinate in coordinate_list:
        if (10000 < coordinate < start) or (end < coordinate < genome_length - 10000):
            result_list.append(coordinate)
            
    return result_list

plus_strand = site_slice(plus_strand_chi, 1200000, 1700000)
minus_strand = site_slice(minus_strand_chi, 1200000, 1700000) 


# Create BED files with intervals around chi-sites (interval length - 500 nt, number of intervals - 40).

# BED file for positive-oriented chi-sites.
chr_name = []
interval_start = []
interval_end = []
interval_name = []

for coordinate in plus_strand:
    coordinate = coordinate - 10000
    step = 500
    interval_number = 1
    while interval_number <= 40:
        chr_name.append('genome')
        interval_start.append(coordinate)
        interval_end.append(coordinate + step)
        interval_name.append(interval_number)
        coordinate += step
        interval_number += 1

df = pd.DataFrame({'chr_name': chr_name, 'interval_start': interval_start, 'interval_end': interval_end, 'interval_name': interval_name})
df.to_csv('../chi_metaplot/plus_intervals.bed', header=False, index=False, sep='\t')

# BED file for negative-oriented chi-sites.
chr_name = []
interval_start = []
interval_end = []
interval_name = []

for coordinate in minus_strand:
    coordinate = coordinate + 9500
    step = 500
    interval_number = 1
    while interval_number <= 40:
        chr_name.append('genome')
        interval_start.append(coordinate)
        interval_end.append(coordinate + step)
        interval_name.append(interval_number)
        coordinate -= step
        interval_number += 1

df = pd.DataFrame({'chr_name': chr_name, 'interval_start': interval_start, 'interval_end': interval_end, 'interval_name': interval_name})
df.to_csv('../chi_metaplot/minus_intervals.bed', header=False, index=False, sep='\t')

print('Chi-sites have been found. BED files are ready.')

plasmid_copy_number = 12 # plasmid copy number

with open('../alignment/aligned_on_plasmid.txt', 'r', encoding='utf-8') as inf:
	for line in inf:
		aligned_on_plasmid = int(line.strip())

with open('../alignment/aligned_on_genome.txt', 'r', encoding='utf-8') as inf:
	for line in inf:
		aligned_on_genome = int(line.strip())

expected = (aligned_on_plasmid + aligned_on_genome) * plasmid_length * plasmid_copy_number / (genome_length + plasmid_length * plasmid_copy_number)

with open('../alignment/statistics.tsv', 'w') as ouf:
    ouf.write('genome length')
    ouf.write('\t')
    ouf.write(str(genome_length))
    ouf.write('\n')
    ouf.write('plasmid length')
    ouf.write('\t')
    ouf.write(str(plasmid_length))
    ouf.write('\n')
    ouf.write('plasmid copy number')
    ouf.write('\t')
    ouf.write(str(plasmid_copy_number))
    ouf.write('\n')
    ouf.write('total number of aligned reads')
    ouf.write('\t')
    ouf.write(str(aligned_on_plasmid + aligned_on_genome))
    ouf.write('\n')
    ouf.write('reads mapped to genome')
    ouf.write('\t')
    ouf.write(str(aligned_on_genome))
    ouf.write('\n')
    ouf.write('reads mapped to plasmid')
    ouf.write('\t')
    ouf.write(str(aligned_on_plasmid))
    ouf.write('\n')
    ouf.write('expected number of reads mapped to plasmid')
    ouf.write('\t')
    ouf.write(str(round(expected)))
    ouf.write('\n')
    ouf.write('real/expected ratio')
    ouf.write('\t')
    ouf.write(str(round((aligned_on_plasmid / expected), 1)))
    ouf.write('\n')

      
