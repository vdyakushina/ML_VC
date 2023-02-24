#!/usr/bin/env python3
# Realign locally around the indels
# Source: https://qcb.ucla.edu/wp-content/uploads/sites/14/2016/03/GATKwr12-3-IndelRealignment.pdf

import sys
import os
from subprocess import run
from func import read_resources

# Initiate
(bam, pref) = sys.argv[1:3]
res = '/pipeline/scripts/resources.csv'

# Check input
if not os.path.isfile(bam):
	print(f'{bam} doesn\'t exist!')
	sys.exit(1)

path = read_resources(res)

# Create target intervals list 
out = pref + ".realigner.intervals"
run("%s -jar %s -T RealignerTargetCreator -R %s -I %s -o %s" %
(path['java8'], path['gatk3.8'], '/pipeline/data/reference/genome_hg19/hg19.fa', bam, out), shell = True)

# Realign reads
inp = out
out = pref + ".realigned.bam"
run("%s -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s --maxConsensuses 990000 --maxReadsForConsensuses 990000 --maxReadsForRealignment 990000 --maxReadsInMemory 990000" %
(path['java8'], path['gatk3.8'], '/pipeline/data/reference/genome_hg19/hg19.fa', bam, inp, out), shell = True)

# Clean
try:
	os.remove(inp)
except:
	print("Error deleting ", inp)

