#!/usr/bin/env python3.8
# Generate pseudo normal sample by random sampling of control tumor samples 
import sys
import os
import csv
from subprocess import run
from func import read_resources
import pprint
import random
import subprocess
import statistics


pp = pprint.PrettyPrinter(indent=4)

# Intitate
(meta, sample) = sys.argv[1:3]

print(meta) # input list of control samples (/home/gkhvorykh/samples/meta.csv)
print(sample) # output .bam file

# read input list of control files

samples_test = []
samples_control = []
with open(meta, newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
	for row in spamreader:
		if row[1] == 'test' :
			samples_test.append(row[0])
		else :
			samples_control.append(row[0])
print(len(samples_test))
samples_control = random.sample(samples_control, 20)

res = "/pipeline/scripts/resources.csv"
path = read_resources(res)

# count reads in all control samples
total_reads = 0
control = dict()
for i in samples_control:
	bam = path['samples'] + '/' + i + '/raw/' + i + ".bam"
	result = subprocess.run([path['samtools'], 'view', '-c', '-L', path['AODABCV1'], bam], stdout=subprocess.PIPE)
	control[i] = int(result.stdout.decode('utf-8'))

target_read_count = int(statistics.median(list(control.values()))/len(samples_control))
t = 0
input_bam = []
seed = int(random.uniform(1000000, 9999999))

# subsample reads from target control bam files
for i in samples_control:
	t += 1
	bam = path['samples'] + '/' + i + '/raw/' + i + ".bam"
	fraction = min(control[i], target_read_count)/control[i]
	fraction = str(int(random.uniform(0, 100)) + fraction)
	output_bam = sample + '.' + str(seed) + '.' + str(t)
	f = open(output_bam, "w")
	subprocess.run([path['samtools'], 'view', '-L', path['AODABCV1'], '-s', fraction, bam, '-b', '-h'], stdout=f)
	input_bam.append(output_bam)

# merge resulting bam files into the final output
subprocess.run(path['samtools'] + ' merge -f ' + sample + ' ' + (' '.join(input_bam)), shell=True)
subprocess.run('rm ' + ' ' + (' '.join(input_bam)), shell=True)
