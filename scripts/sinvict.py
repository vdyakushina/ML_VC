#!/usr/bin/env python3
# A wrapper to run sinvict tool
# Article: https://pubmed.ncbi.nlm.nih.gov/27531099/
# Source: https://github.com/sfu-compbio/sinvict/

import argparse
import os
import sys
from subprocess import run
from tempfile import TemporaryDirectory
from func import read_resources
from func import ts
from func import check_log

# Initiate the variables
resources = "/pipeline/scripts/resources.csv"
path = read_resources(resources)
pref="sinvict"
cpu = 8

# Process arguments
p = argparse.ArgumentParser(description='SNV and indels calling with SinVICT')
p.add_argument('bam', help="path/to/bam file")
p.add_argument('target_id', help="target id")
p.add_argument('folder', help="folder to save output")
p.add_argument('sample', help="sample id")

args = p.parse_args()
bam = args.bam

target = path.get(args.target_id)
if not target:
	print(f'{args.target_id} is unknown target id')
	sys.exit(1)

sample = args.sample
folder = args.folder
ref = path['ref']

# Check input
if not os.path.isfile(bam): sys.exit(f'Error: {bam1} doesn\'t exist!')
if not os.path.isfile(target): sys.exit(f'Error: {target} doesn\'t exist!')
if not os.path.isfile(ref): sys.exit(f'Error: {ref} doesn\'t exist!')

# Create the output folder
if not os.path.isdir(folder):
	os.makedirs(folder)

# Generate a readcount file in the format required by SinVICT
print(ts(), "Generate readcount")
run("%s -w10 -f %s -l %s %s > %s/%s.rc" %
	(path['readcount'], ref, target, bam, folder, sample), shell = True)

# Run SinVICT
print(ts(), "Run SinVICT")
run("%s -t %s -o %s" % (path['sinvict'], folder, folder), shell = True)

