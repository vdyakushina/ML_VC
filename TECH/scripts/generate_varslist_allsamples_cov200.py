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

(variant, infolder, out_folder, maf) = sys.argv[1:5]

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

with open("%s/%s:%s%s>%s" % (out_folder,chr,start,refl,altl), "w") as file:
	file.write("VAF:sample")

for sample in samples:
	# Run SGA
	print(sample)
	result = run("%s somatic-variant-filters --tumor-bam %s/%s/analysis/provisional.bam --normal-bam %s/%s/analysis/provisional.bam --reference %s %s/%s:%s%s\">\"%s.vcf -t 10 --annotate-only" % \
(sga,infolder,sample,infolder,sample,ref,out_folder,chr,start,refl,altl),  stdout=PIPE, shell=True)
	DP=int(float(result.stdout.decode('utf-8').split('\t')[7].split('TumorTotalDepth=')[1].split(';')[0]))
	VAF=float(result.stdout.decode('utf-8').split('\t')[7].split('TumorVAF=')[1].split(';')[0])
	print(f'DP={DP}, AF={VAF}')
	if DP >= 200:
		with open("%s/%s:%s%s>%s" % (out_folder,chr,start,refl,altl), "a") as file:
			file.write('\n'+':'.join([str(VAF),sample]))

run("sed \'1d\' %s/%s:%s%s\">\"%s | tr . , | sort -nr -t':' -k1 | tr , . > %s/%s:%s%s\">\"%s_sorted" % (out_folder,chr,start,refl,altl,out_folder,chr,start,refl,altl), shell=True)
run ("mv %s/%s:%s%s\">\"%s_sorted %s/%s:%s%s\">\"%s" % (out_folder,chr,start,refl,altl,out_folder,chr,start,refl,altl), shell=True)
run("rm %s/%s:%s%s\">\"%s.vcf" % (out_folder,chr,start,refl,altl), shell=True)


run("sed -i 1i\'MAF=%s\' %s\'%s\'" % (maf, out_folder, variant), shell = True)