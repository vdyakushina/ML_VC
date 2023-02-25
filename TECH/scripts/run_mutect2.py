#!/usr/bin/env python3
# Wrapper script to process the whole pipeline

from subprocess import run
from pathlib import Path
from func import read_resources
import sys
import os.path
import csv

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
target = 'AODABCV1'

samples=[]
with open(path['meta'], newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
	for row in spamreader:
		samples.append(row[0])

for sample in samples:
	bam='/home/gkhvorykh/samples/%s/analysis/provisional.bam' % sample
	out='/home/gkhvorykh/samples/%s' % sample
	print(f'sample: {sample}, bam: {bam}, out:{out}')
	run("%s/mutect2_to.py %s %s %s %s" % (path['scripts'], bam, out + "/", path['ref'], target), shell = True)