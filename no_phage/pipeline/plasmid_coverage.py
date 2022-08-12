#!/usr/bin/python3

import pandas as pd

window = 200
step = 10

plasmid = pd.read_csv('../reference/plasmid.fa', sep='\t')
plasmid.columns = ['fasta']
plasmid = str(plasmid['fasta'].sum()).upper()
plasmid_length = len(plasmid)

coordinate = []
start = []
end = []

c = 0
while c <= plasmid_length:
    coordinate.append(c)
    start.append(c - window/2)
    end.append(c + window/2)
    c += step
    
intervals = pd.DataFrame({'coordinate' : coordinate, 'start': start, 'end': end})

plus = pd.read_csv('../plasmid/plus_coverage.tsv', sep='\t', header=None)
plus.columns = ['chr', 'coordinate', 'coverage']
plus = plus [['coordinate', 'coverage']]

upper_coverage = list(plus.tail(int(window/2))['coverage'])

upper_coordinate = []
for i in range (-int(window/2), 0):
    upper_coordinate.append(i)

upper_plus = pd.DataFrame({'coordinate' : upper_coordinate, 'coverage' : upper_coverage})

lower_coverage = list(plus.head(int(window/2))['coverage'])

lower_coordinate = []
for i in range (int(plus.tail(1)['coordinate']) + 1, int(plus.tail(1)['coordinate']) + int(window/2+1)):
    lower_coordinate.append(i)
    
lower_plus = pd.DataFrame({'coordinate' : lower_coordinate, 'coverage' : lower_coverage})

plus = pd.concat([upper_plus, plus, lower_plus], ignore_index=True, sort=False)

minus = pd.read_csv('../plasmid/minus_coverage.tsv', sep='\t', header=None)
minus.columns = ['chr', 'coordinate', 'coverage']
minus = minus[['coordinate', 'coverage']]

upper_coverage = list(minus.tail(int(window/2))['coverage'])

upper_coordinate = []
for i in range (-int(window/2), 0):
    upper_coordinate.append(i)
    
upper_minus = pd.DataFrame({'coordinate' : upper_coordinate, 'coverage' : upper_coverage})

lower_coverage = list(minus.head(int(window/2))['coverage'])

lower_coordinate = []
for i in range (int(minus.tail(1)['coordinate']) + 1, int(minus.tail(1)['coordinate']) + int(window/2+1)):
    lower_coordinate.append(i)

lower_minus = pd.DataFrame({'coordinate' : lower_coordinate, 'coverage' : lower_coverage})

minus = pd.concat([upper_minus, minus, lower_minus], ignore_index=True, sort=False)

def calculate_plus_coverage(row):
    '''
    This function calculates coverage of positive strand in sliding windows.
    '''
    cov = plus.query('coordinate >= @row["start"] and coordinate < @row["end"]')['coverage'].mean()
    return cov

intervals['plus_coverage'] = intervals.apply(calculate_plus_coverage, axis=1)

def calculate_minus_coverage(row):
    '''
    This function calculates coverage of negative strand in sliding windows.
    '''
    cov = minus.query('coordinate >= @row["start"] and coordinate < @row["end"]')['coverage'].mean()
    return cov

intervals['minus_coverage'] = intervals.apply(calculate_minus_coverage, axis=1)

intervals.to_csv('../plasmid/coverage.tsv', sep='\t')
