#!/usr/bin/env python3
# Wrapper script to process the whole pipeline

from subprocess import run
from pathlib import Path
from func import read_resources
import sys
import os.path
import csv
import multiprocessing as mp
import pandas as pd

res = os.getcwd()+"/resources.csv"
path = read_resources(res)

out_folder='/home/gkhvorykh/train_positive/positive_fastq/'

# Create dictionary {VarName:samples} from File1.tsv (file with positive samples)
var_dict={}

with open('/home/gkhvorykh/train_positive/File1_Max.tsv', 'r') as inf:
	for line in inf:
		line=line.strip()
		if len(line.split('\t')) > 1:
			var_dict[line.split('\t')[0]]=line.split('\t')[1].split(',')


# Run generate_fastq_positive.py for each variant from var_dict
for variant in var_dict.keys():
	samples=var_dict[variant]
	print(variant)
	print(samples)
	for sample in samples:
		run("./generate_fastq_positive.py '%s' /home/gkhvorykh/samples/%s/raw/%s.bam %s" % (variant, sample, sample, out_folder), shell = True)

