#!/usr/bin/env python3
# Parser sample.mosdepth.region.dist.txt for breadth10x and breadth20x 
import sys
import re

# Initiate
file = sys.argv[1]

# Parser the input file
with open(file, 'r') as f:
    for line in f:
        line = line.rstrip()
        if re.search("^total\t10\t", line):
            breadth10 = float(line.split()[2])
        if re.search("^total\t20\t", line):
            breadth20 = float(line.split()[2])

# Output the values found
print("%.2f %.2f" % (breadth10, breadth20))