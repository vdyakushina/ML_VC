#!/usr/bin/env python3.8
from subprocess import run
from subprocess import PIPE
import re
import sys
import os
from func import read_resources
import multiprocessing as mp

(variant, bam, out_folder) = sys.argv[1:4]

if not os.path.exists(out_folder):
  run("mkdir %s" % out_folder, shell=True)

out_folder=out_folder
chr=variant.split(':')[0]
start=re.split(r"[A-Z]", variant)[0].split(':')[1]
(refl,altl)=re.findall(r"[A-Z]+", variant)

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
samtools = path['samtools']
sga = path['sga']
ref = path['ref']

# Generate bam with reads from variant region

print(f'get names')
run("%s view -@10 /home/gkhvorykh/samples/%s/analysis/provisional.bam %s:%s-%s | cut -f 1 | sort | uniq > %s/%s:%s%s\">\"%s.%s.names.txt" % (samtools, bam, chr, start, start, out_folder, chr,start,refl,altl,bam), shell=True)
print(f'get header')
run("%s view -H /home/gkhvorykh/samples/%s/analysis/provisional.bam > %s/%s:%s%s\">\"%s.%s.sam" % (samtools, bam, out_folder, chr,start,refl,altl,bam), shell=True)
print(f'get reads')
run("%s view -@10 /home/gkhvorykh/samples/%s/analysis/provisional.bam %s | grep -f %s/%s:%s%s\">\"%s.%s.names.txt >> %s/%s:%s%s\">\"%s.%s.sam" % (samtools, bam, chr, out_folder, chr,start,refl,altl,bam, out_folder, chr,start,refl,altl,bam), shell=True)
print(f'fastq')
run("%s fastq -1 %s/%s:%s%s\">\"%s.tmp_1_%s.fq -2 %s/%s:%s%s\">\"%s.tmp_2_%s.fq %s/%s:%s%s\">\"%s.%s.sam" % (samtools,out_folder, chr,start,refl,altl, bam,out_folder, chr,start,refl,altl, bam,out_folder, chr,start,refl,altl, bam), shell=True)
run("cat %s/%s:%s%s\">\"%s.tmp_1_%s.fq >> %s/%s:%s%s\">\"%s.negative.R1.fq" % (out_folder, chr,start,refl,altl, bam,out_folder,chr,start,refl,altl), shell=True)
run("cat %s/%s:%s%s\">\"%s.tmp_2_%s.fq >> %s/%s:%s%s\">\"%s.negative.R2.fq" % (out_folder, chr,start,refl,altl, bam,out_folder,chr,start,refl,altl), shell=True) 
run("rm %s/%s:%s%s\">\"%s.%s.sam %s/%s:%s%s\">\"%s.%s.names.txt %s/%s:%s%s\">\"%s.tmp_1_%s.fq %s/%s:%s%s\">\"%s.tmp_2_%s.fq" % (out_folder, chr,start,refl,altl, bam, out_folder, chr,start,refl,altl, bam, out_folder, chr,start,refl,altl, bam, out_folder, chr,start,refl,altl, bam), shell=True)
