#!/usr/bin/env python3
# Finds all fastq files in subfolders, get sample ids, and aligns the reads. 

import sys
import re
import csv
import subprocess
from pathlib import Path
from shutil import move
from func import read_resources

# Initiate
meta = {}
resources = "/pipeline/scripts/resources.csv"
rs = read_resources(resources)

# Find all fastq.gz files in subfolders of the current folder
for path in Path().glob('**/*fastq.gz'):
	p1 = str(path.resolve())
	# Get basename
	fq = path.name
	# Create regular expression
	p = re.compile('[a-z0-9]*', re.I)
	# Match the pattern at the begining of the basename
	res = p.match(fq)
	if res == None: continue
	sample = res.group()
	if sample == "Undetermined": continue
	# Construct meta 
	if sample not in meta: meta[sample] = []
	meta[sample].append(p1)

for sample in meta.keys():
	fq =  sorted(meta[sample])
	# Skip if the bam file already exist
	if Path(f'{sample}/{sample}.bam').is_file(): continue
	print(f'Processing {sample}')
	subprocess.check_call("%s/align.sh %s %s %s" % (rs['scripts'], fq[0], fq[1], sample),
		shell = True)
	# Move bam files
	for ext in ['.bam', '.bam.bai']:
		move(sample+ext, sample)
