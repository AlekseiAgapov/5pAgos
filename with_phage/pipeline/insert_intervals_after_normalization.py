#!/usr/bin/python3
# This script takes for input normalized.tsv file and intervals.bed file. If any intervals are skipped in the first file, they are taken from the second and count is filled with 0.

import pandas as pd

# Read the file with normalized interval coverage.
normalized = pd.read_csv('./normalized.tsv', sep='\t')

# Read the BED file with original intervals.
intervals = pd.read_csv('./intervals.bed', sep='\t', header=None)
intervals.columns = ['chr', 'start', 'end', 'interval_name']

# Merge the two dataframes and fill NA with 0.
intervals = intervals.merge(normalized, how='left', on='interval_name')
intervals = intervals.fillna(0)
result = intervals[['interval_name', 'interval_length', 'coverage']]

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
