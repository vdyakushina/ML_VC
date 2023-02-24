#!/usr/bin/env python3
# Wrapper script to process the Interpretation

from subprocess import run
from pathlib import Path
from func import read_resources
import sys
import os.path
import csv
import multiprocessing as mp

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
target = 'AODABCV1'

samples=[]
with open(path['meta'], newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
	for row in spamreader:
		samples.append(row[0])

def run_Interpretation_script(sample):
	input="/home/gkhvorykh/samples/%s/variant.mutect2_to.vcf.gz" % (sample)
#	if not os.path.exists("/home/gkhvorykh/samples/%s/analysis/variant.raw.SNV.tsv" % sample):
#		run("echo \"%s\" >> /home/gkhvorykh/samples/samples_not_ready_for_Interpretation.txt" % sample, shell = True)
#	else:
	run("%s/Interpretation.py %s" % (path['scripts'], input), shell = True)
#		run("mv /home/gkhvorykh/samples/%s/analysis/Interpretation_result.tsv /home/gkhvorykh/samples/%s/analysis/Interpretation_result_%s.tsv" % (sample, sample, type), shell = True)
		
for sample in samples:
	run_Interpretation_script(sample)
 
