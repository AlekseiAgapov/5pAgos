#!/usr/bin/python3
# This script takes for input intersected.tsv file with reads intersected with intervals and counts.tsv file with information about number of sites where each read is mapped to.
# Next, each read receives a count as 1/(number of sites). It's always 1 for uniquely mapped reads and a number between > 0 and < 1 for a multimapper.
# The script prints as output a tab-separated table with two columns: interval name and number of aligned reads.

import pandas as pd

# Read the file with intersection.
intersected = pd.read_csv('./intersected.tsv', sep='\t', header=None)
intersected = intersected.loc[:, (1,2,3,7)]
intersected.columns = ['start', 'end', 'interval_name', 'read']
intersected['interval_length'] = intersected['end'] - intersected['start']

# Read the file with information about the number of sites where a read is mapped to.
counts = pd.read_csv('../alignment/counts.tsv', header=None, sep='\t')
counts.columns = ['count', 'read']

# Merge the two tables and calculate a value of each read.
intersected = intersected.merge(counts, how='left', on='read')
intersected['value'] = 1 / intersected['count']

# Calculate the sum of mapped reads values for each interval.
result = intersected.pivot_table(index=['start', 'end', 'interval_name', 'interval_length'], values='value', aggfunc='sum').reset_index()
result.columns = ['start', 'end', 'interval_name', 'interval_length', 'coverage']

# Make three lists from a DataFrame.
intervals = list(result['interval_name'])
length = list(result['interval_length'])
coverage = list(result['coverage'])

# Print out the result.

print('interval_name' + '\t' + 'interval_length' + '\t' + 'coverage')

for row in range(len(intervals)):
    print(intervals[row], end='\t')
    print(length[row], end='\t')
    print(coverage[row])
