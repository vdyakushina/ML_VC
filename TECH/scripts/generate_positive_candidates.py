#!/usr/bin/env python3.8
from subprocess import run
from subprocess import PIPE
import re
import sys
import os
from func import read_resources
import multiprocessing as mp


(infolder, out) = sys.argv[1:3]

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
meta=path['meta']
gnomad_genome=path['gnomad_genome_abc']
gnomad_exome=path['gnomad_exome_abc']

samples=[]
with open('%s' % meta, 'r') as in_s:
	for line in in_s:
		samples.append(line.split(',')[0])

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


old=[]
with open('/pipeline/TECH/data/vcf_to_txt/vcf_to_txt.maf_panel.maf_gnomad.short_variants.full.tsv', 'r') as in_o:
	next(in_o)
	for line in in_o:
		line=line.strip()
		varo=line.split('\t')[0]
		old.append(varo)


vars_count={}
vars_gnomad={}
for sample in samples:
	with open('%s/%s/analysis/variant.raw.INDEL.tsv' % (infolder, sample), 'r') as in_t:
		next(in_t)
		for line in in_t:
			line=line.strip()
			var=str(line.split('\t')[0])
			try:
				if (float(line.split('\t')[13])>=0.2 and float(line.split('\t')[14])>=0.5 and float(line.split('\t')[17])>=500) and var not in old:
					if var in vars_count:
						vars_count[var]+=1
					else:
						vars_count[var]=1
						if var in gnomad:
							vars_gnomad[var]=str(gnomad[var])
						else:
							vars_gnomad[var]='N/A'
			except (ValueError):
				continue


with open ('%s' % out, 'w') as out:
	out.write('VarName' + '\t' + 'Number of samples' + '\t'+ 'MAF_gnomAD(max(exome,genome))' +'\n')
	for var in vars_count:
		out.write('\t'.join([str(var), str(vars_count[var])])+'\t'+str(vars_gnomad[var])+'\n')