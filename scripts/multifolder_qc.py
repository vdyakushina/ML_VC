#!/usr/bin/env python3
# Finds all bam files in subfolders, get sample ids, and
# process qc.

import os
import re
import subprocess
from pathlib import Path
from shutil import move
from func import read_resources

# Initiate
meta = {}
resources = "/pipeline/scripts/resources.csv"
rs = read_resources(resources)

# Get working directory 
wd = os.getcwd()

# Find all files in subfolders of the current folder
for path in Path().glob('**/*.bam'):
	bam = str(path.resolve())
	# Get basename
	fn = path.name
	# Create regular expression
	p = re.compile('[a-z0-9]*', re.I)
	# Match the pattern at the begining of the basename
	res = p.match(fn)
	if res == None: continue
	sample = res.group()
	if sample == "Undetermined": continue
	# Construct meta 
	meta[sample] = bam

for sample in meta.keys():
	bam =  meta[sample]
	print(f'Processing {sample}')
	# Step in the folder with sample
	os.chdir(sample)
	subprocess.check_call("%s/qc.py %s %s" % (rs['scripts'], bam, sample),
		shell = True)
	os.chdir(wd)
