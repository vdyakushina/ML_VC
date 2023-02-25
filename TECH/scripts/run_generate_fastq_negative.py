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

# Create dictionary {VarName:samples} from File2.tsv (file with neagtive samples)
var_dict={}
with open('/home/gkhvorykh/train_negative/File2_Max.tsv', 'r') as inf:
	for line in inf:
		line=line.strip()
		if len(line.split('\t')) > 1:
			var_dict[line.split('\t')[0]]=line.split('\t')[1].split(',')

# Run generate_fastq_negative.py for each variant from var_dict

def run_func(variant, out_folder):
	samples=var_dict[variant]
	print(variant)
	print(samples)
	for sample in samples:
		run("./generate_fastq_negative.py '%s' %s %s" % (variant, sample, out_folder), shell = True)


pool = mp.Pool(10)
pool.starmap(run_func, [(variant, '/home/gkhvorykh/train_negative/negative_fastq/') for variant in var_dict.keys()])
pool.close()