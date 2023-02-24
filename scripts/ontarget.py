#!/usr/bin/env python3

# Estimate the coverage on target 

import sys, subprocess, pandas as pd

# Initiate
bam, target = sys.argv[1:3]
resources = '/home/onco-admin/ATLAS_software/aod-pipe/Pipe/AODTKAbatch/pipeline/scripts/resources.csv'
samtools='samtools'

# Get resources
def read_resources(csv):
	df = pd.read_csv(csv, header = None, usecols = [0, 1], names = ["id", "path"])
	return df.set_index('id').T.to_dict('records')[0]

path = read_resources(resources)

# Estimate percentage of reads mapped on target 
total = subprocess.check_output(
	"%s view -q 1 -c %s" % (samtools, bam),
	shell = True)
on = subprocess.check_output(
	"%s view -q 1 -c -L %s %s" % (samtools, path[target], bam),
	shell = True)
ontarget = round(int(on)/int(total), 2)

# Output
print({'ontarget': ontarget})
