#!/usr/bin/env python3
# Run mutect2 in tumor only mode 
import sys
import os
from subprocess import run
from func import read_resources

# Intitate
if len(sys.argv) == 3:
	(bam, sample) = sys.argv[1:3]
	mode = "nobam"
elif len(sys.argv) == 4:
	(bam, sample, mode) = sys.argv[1:4]
else:
	print("Missing arguments...")
	sys.exit(1)

res = "/pipeline/scripts/resources.csv"
path = read_resources(res)

# Check input
if not os.path.isfile(bam):
	print(bam, "doesn't exist")
	sys.exit(1)

print(f'bam: {bam}\nsample: {sample}\nmode: {mode}')

# Run mutect2 in tumor only mode
if mode == 'bam':
	# Run mutect2 with bam file creation
	print("Run mutect2 with bam file creation")
	output = sample + ".mutect2.vcf.gz"
	bamout = sample + ".mutect2.bam"
	run("%s --java-options \"-Xmx8G\" Mutect2 \
		-R %s \
		-I %s \
		-O %s \
		--max-reads-per-alignment-start 0 \
		-bamout %s" %
		(path['gatk'], path['ref'], bam, output, bamout), shell = True)
else:
	# Run mutect2 without creating bam file
	print("Run mutect2 without bam file creation")
	output = sample + ".mutect2.vcf.gz"
	run("%s --java-options \"-Xmx8G\" Mutect2 \
		-R %s \
		-I %s \
		--max-reads-per-alignment-start 0 \
		-O %s" %
		(path['gatk'], path['ref'], bam, output), shell = True)

