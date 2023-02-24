#!/usr/bin/env python3
import gzip
import sys

# Initiate
bed = sys.argv[1]
avgcov = float(sys.argv[2])
n, total = 0, 0
cov = 1.2 * avgcov

# Read file and process the data
with gzip.open(bed, 'rt') as f:
    for line in f:
        x = line.split()
        if float(x[4]) > cov: n += 1
        total += 1 

# Output uniformity 
unf = n/total*100
print("%.1f" % unf)
