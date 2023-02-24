#!/usr/bin/env python3
# Copy files from one folder to another. 
# Arrange the files at new location  according to the file structure recommended
# Write meta file with from/to pathes 

import sys
import re
import csv
from pathlib import Path
from shutil import copyfile

# Initiate
(in_folder, out_folder) = sys.argv[1:3]
meta_file = "meta.csv"
meta = {}

# Show input
print(f"Input folder: {in_folder}\nOutput folder: {out_folder}")

# Make output folder
Path(out_folder).mkdir(
	parents = True, exist_ok = True)

# Glob all subfolders with fastq.gz files in the input folder
for path in Path(in_folder).glob('**/*fastq.gz'):
	p1 = str(path)
	# Get basename
	fq = path.name
	# Create regular expression
	p = re.compile('[a-z0-9_]*_L001', re.I)
	# Match the pattern at the begining of the basename
	res = p.match(fq)
	if res == None: continue
	sample = res.group()
	if sample == "Undetermined": continue
	# Construct new path
	p2 = f'{out_folder}/{sample}/raw/{fq}'
	if sample not in meta: meta[sample] = []
	meta[sample].append([p1, p2])

# Open meta file
with open(f'{out_folder}/{meta_file}', 'w') as mf:
	w = csv.writer(mf)
	# Loop via all samples and copy the files 
	for sample in meta.keys():
		# Make folder for raw data
		raw = f'{out_folder}/{sample}/raw'
		Path(raw).mkdir(
			parents = True, exist_ok = True)
		# Write meta and copy files 
		for p in meta[sample]:
			copyfile(p[0], p[1])
			w.writerow(p)
