#!/usr/bin/env python3

from subprocess import run
from pathlib import Path
from func import read_resources
import sys
import os.path
import csv
import multiprocessing as mp
import numpy as np

res = os.getcwd()+"/resources.csv"
path = read_resources(res)

vars_exist=[]
with open('/home/gkhvorykh/train_true_samples/train_true_results_SNV.tsv', 'r') as inf:
	next(inf)
	for line in inf:
		line=line.strip()
		var=line.split('\t')[0]
		vars_exist.append(var)

with open('/home/gkhvorykh/train_true_samples/train_true_results_INDEL.tsv', 'r') as inf:
	next(inf)
	for line in inf:
		line=line.strip()
		var=line.split('\t')[0]
		vars_exist.append(var)

vars=[]
with open('/home/gkhvorykh/train_true_samples/meta_1.tsv', 'r') as inf:
	for line in inf:
		line=line.strip()
		var=line.split('\t')[0]
		if var not in vars_exist:
			vars.append(var)

with open('/home/gkhvorykh/train_true_samples/meta_05.tsv', 'r') as inf:
	for line in inf:
		line=line.strip()
		var=line.split('\t')[0]
		if var not in vars_exist:
			vars.append(var)

with open('/home/gkhvorykh/train_true_samples/meta_5.tsv', 'r') as inf:
	for line in inf:
		line=line.strip()
		var=line.split('\t')[0]
		if var not in vars_exist:
			vars.append(var)

with open('/home/gkhvorykh/train_true_samples/meta_25.tsv', 'r') as inf:
	for line in inf:
		line=line.strip()
		var=line.split('\t')[0]
		if var not in vars_exist:
			vars.append(var)

with open('/home/gkhvorykh/train_true_samples/meta_ba.tsv', 'r') as inf:
	for line in inf:
		line=line.strip()
		var=line.split('\t')[0]
		if var not in vars_exist:
			vars.append(var)
'''
def run_func(variantf):
	variant=variantf.split('_')[0]+':'+variantf.split('_')[1]+'>'+variantf.split('_')[2]
	f=variantf.split('_')[-2]
	depth=variantf.split('_')[-1]
	run("./generate_fastq_1j_temp.py '%s' %s %s" % (variant, f, depth), shell = True)


pool = mp.Pool(4)
pool.map(run_func, [variantf for variantf in vars])
pool.close()
'''