#!/usr/bin/env python3.8
from subprocess import run
from subprocess import PIPE
import csv
import re
import sys
import os
from func import read_resources
import multiprocessing as mp
import pandas as pd

try:
	(variant, target_sample, infolder, out_folder, maf) = sys.argv[1:6]
except ValueError:
	res = os.getcwd()+"/resources.csv"
	path = read_resources(res)
	meta=path['meta']
	gnomad_genome=path['gnomad_genome_abc']
	gnomad_exome=path['gnomad_exome_abc']
	(variant, target_sample, infolder, out_folder) = sys.argv[1:5]
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

	if variant in gnomad:
		maf=gnomad[variant]
	else:
		maf=0.0

if not os.path.exists(out_folder):
  run("mkdir %s" % out_folder, shell=True)

chr=variant.split(':')[0]
start=re.split(r"[A-Z]", variant)[0].split(':')[1]
(refl,altl)=re.findall(r"[A-Z]+", variant)

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
samtools = path['samtools']
sga = path['sga']
ref = path['ref']

# Generate vcf with target variant for SGA
run("echo \"%s\t%s\t.\t%s\t%s\t.\t.\t.\t.\t.\" > %s/%s:%s%s\">\"%s.vcf" % (chr,start,refl,altl,out_folder,chr,start,refl,altl), shell=True)

samples=[]
with open(path['meta'], newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
	for row in spamreader:
		samples.append(row[0])

def run_sga(sample): # Run SGA
	print(sample)
	result = run("%s somatic-variant-filters --tumor-bam %s/%s/analysis/provisional.bam --normal-bam %s/%s/analysis/provisional.bam --reference %s %s/%s:%s%s\">\"%s.vcf -t 10 --annotate-only" % \
(sga,infolder,sample,infolder,sample,ref,out_folder,chr,start,refl,altl),  stdout=PIPE, shell=True)
	DP=int(float(result.stdout.decode('utf-8').split('\t')[7].split('TumorTotalDepth=')[1].split(';')[0]))
	VAF=float(result.stdout.decode('utf-8').split('\t')[7].split('TumorVAF=')[1].split(';')[0])
	print(f'DP={DP}, AF={VAF}')
#	if DP >= 200:
	with open("%s/%s:%s%s>%s" % (out_folder,chr,start,refl,altl), "a") as file:
		file.write('\n'+':'.join([str(VAF),sample]))

if not os.path.exists("%s/%s:%s%s>%s" % (out_folder,chr,start,refl,altl)):
	pool = mp.Pool(10)
	pool.map(run_sga, [str(sample) for sample in samples])
	pool.close()
	run("grep -v %s %s/%s:%s%s\">\"%s | tr . , | sort -nr -t':' -k1 | tr , . > %s/%s:%s%s\">\"%s_sorted" % (target_sample, out_folder,chr,start,refl,altl, out_folder,chr,start,refl,altl), shell=True)
	run("sed -i 1i`grep %s %s/%s:%s%s\">\"%s` %s/%s:%s%s\">\"%s_sorted" % (target_sample, out_folder,chr,start,refl,altl, out_folder,chr,start,refl,altl), shell=True)
	run ("mv %s/%s:%s%s\">\"%s_sorted %s/%s:%s%s\">\"%s" % (out_folder,chr,start,refl,altl,out_folder,chr,start,refl,altl), shell=True)
	run("rm %s/%s:%s%s\">\"%s.vcf" % (out_folder,chr,start,refl,altl), shell=True)
	run("sed -i 1i\'MAF=%s\' %s\'%s\'" % (maf, out_folder, variant), shell = True)
else:
	print("%s/%s:%s%s>%s Is already exist" % (out_folder,chr,start,refl,altl))
