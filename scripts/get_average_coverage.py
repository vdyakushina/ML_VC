#!/usr/bin/env python3
import gzip
import sys

# Initiate
file = sys.argv[1]
cov = 0
length = 0

# Read file and process the data
with gzip.open(file, 'rt') as f:
    for line in f:
        x = line.split()
        l = int(x[2]) - int(x[1])
        cov += l * float(x[4])
        length += l 

# Output average coverage
print("%.2f" % (cov/length))
