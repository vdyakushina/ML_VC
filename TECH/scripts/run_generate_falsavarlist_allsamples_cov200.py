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
meta=path['meta']
gnomad_genome=path['gnomad_genome_abc']
gnomad_exome=path['gnomad_exome_abc']

gnomad={}
with open('%s' % gnomad_exome, 'r') as in_g:
	next(in_g)
	for line in in_g:
		line=line.strip()
		gnomad[''.join([':'.join(['chr'+str(line.split('\t')[0]), line.split('\t')[1]]), '>'.join([line.split('\t')[3], line.split('\t')[4]])])]=line.split('\t')[5]


with open('%s' % gnomad_genome, 'r') as in_g:
	next(in_g)
	for line in in_g:
		line=line.strip()
		varg=str(''.join([':'.join(['chr'+str(line.split('\t')[0]), line.split('\t')[1]]), '>'.join([line.split('\t')[3], line.split('\t')[4]])]))
		if varg in gnomad:
			varg_old=gnomad[varg]
			gnomad[varg]=max(float(line.split('\t')[5]), float(varg_old))
		else:
			gnomad[varg]=float(line.split('\t')[5])


variants=[]
mafs={}

with open('/home/gkhvorykh/train_positive/File1_Max.tsv', 'r') as file:
	#next(file)
	lines = csv.reader(file, delimiter='\t')
	for line in lines:
		variants.append(line[0])
		try:
			mafs[line[0]]=gnomad[line[0]]
		except KeyError:
			mafs[line[0]]=0
variants=variants[1:]

def run_func(variant, infolder, out_folder):
	print(variant)
	print(mafs[variant])
	if not os.path.exists(out_folder+variant):
		run("./generate_varslist_allsamples_cov200.py \"%s\" %s %s %s" % (variant, infolder, out_folder, mafs[variant]), shell = True)

'''
pool = mp.Pool(1)
pool.starmap(run_func, [(variant, '/home/gkhvorykh/samples/', '/home/gkhvorykh/train_false/') for variant in variants])
pool.close()
'''
out_folder='/home/gkhvorykh/train_false/'
for variant in variants:
	if os.path.exists(out_folder+variant):
		try:
			var_f=pd.read_csv('%s%s' % (out_folder, variant), header=None).rename({0:variant}, axis=1)
			if not os.path.exists(out_folder+'all_false_SNV_VAF_samples.xls'):
				var_f.to_excel('%sall_false_SNV_VAF_samples.xls' % out_folder, index=False)
			else:
				xls=pd.read_excel('%sall_false_SNV_VAF_samples.xls' % out_folder)
				xls[variant]=var_f
				xls.to_excel('%sall_false_SNV_VAF_samples.xls' % out_folder, index=False, engine='xlsxwriter')
		except pd.errors.EmptyDataError:
			with open('%s/empty.tsv' % out_folder, 'a') as empty_out:
				empty_out.write(variant+'\n')
	else:
		print('not ready')