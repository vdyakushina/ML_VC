#!/usr/bin/env python3
# Estimates the average coverage and uniformity 
# from the file generated with 'samtools coverage' command 

import sys
import pandas as pd

# Initiate variables
file=sys.argv[1]

# Load data
df = pd.read_csv(file, sep = "\t", header = None, prefix = "v")
# Get average coverage
avgcov = df.v6.mean()
# Get uniformity
n = sum(df["v6"] >= avgcov * 0.2)
m = df.shape[0]
print("%i %.2f " % (int(avgcov), n/m))
