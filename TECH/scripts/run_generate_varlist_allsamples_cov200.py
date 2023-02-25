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
from func import read_resources
from subprocess import run

res = os.getcwd()+"/resources.csv"
path = read_resources(res)

variants=[]
mafs={}
with open('/pipeline/TECH/data/positive_candidats_INDEls.txt', 'r') as file:
	lines = csv.reader(file, delimiter='\t')
	for line in lines:
		variants.append(line[0])
		mafs[line[0]]=line[2]

variants=variants[1:]

def run_func(variant, infolder, out_folder):
	print(variant)
	print(mafs[variant])
	run("/pipeline/TECH/scripts/generate_varslist_allsamples_cov200.py \"%s\" %s %s %s" % (variant, infolder, out_folder, mafs[variant]), shell = True)

pool = mp.Pool(10)
pool.starmap(run_func, [(variant, '/home/gkhvorykh/samples/', '/home/gkhvorykh/samples/varians_all_samples/') for variant in variants])
pool.close()


out_folder='/home/gkhvorykh/samples/varians_all_samples/'
for variant in variants:
	if os.path.exists(out_folder+variant):
		var_f=pd.read_csv('%s%s' % (out_folder, variant), header=None).rename({0:variant}, axis=1)
		if not os.path.exists(out_folder+'all_variants_VAF_samples_expanded_INDEL.xls'):
			var_f.to_excel('%sall_variants_VAF_samples_expanded_INDEL.xls' % out_folder, index=False)
		else:
			xls=pd.read_excel('%sall_variants_VAF_samples_expanded_INDEL.xls' % out_folder)
			xls[variant]=var_f
			xls.to_excel('%sall_variants_VAF_samples_expanded_INDEL.xls' % out_folder, index=False)
	else:
		print('not ready')
