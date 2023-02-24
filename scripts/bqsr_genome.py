#!/usr/bin/env python3
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-
# Run Base Quality Score Recalibration (BQSR)

import sys
from subprocess import run
from func import read_resources
from pathlib import Path

(bam, pref) = sys.argv[1:3]

# Initiate
res = "/pipeline/scripts/resources.csv"
path = read_resources(res)

# Check input
if not Path(bam).is_file():
	print(f'{bam} doesn\'t exist')
	sys.exit(1)

# Create the recalibration table
out = pref + ".recall.table"
run("%s BaseRecalibrator -I %s -R %s --known-sites %s -O %s" %
	(path['gatk'], bam, '/pipeline/data/reference/genome_hg19/hg19.fa', path['dbsnp'], out), shell = True)

# Apply base quality score recalibration
inp = out
out = pref + ".bqsr.bam"
run("%s ApplyBQSR -I %s -R %s --bqsr-recal-file %s -O %s" %
	(path['gatk'], bam, '/pipeline/data/reference/genome_hg19/hg19.fa', inp, out), shell = True)


