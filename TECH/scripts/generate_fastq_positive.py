#!/usr/bin/env python3.8
from subprocess import run, PIPE
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

# Generate vcf with target variant for SGA
run("echo \"%s\t%s\t.\t%s\t%s\t.\t.\t.\t.\t.\" > %s/target.vcf" % (chr,start,refl,altl,out_folder), shell=True)

# Generate bam with reads from variant region & get Read Names for furthe analysis with SGA
run("%s view %s %s:%s-%s > %s/target.sam" % (samtools, bam, chr, start, start, out_folder), shell=True)
with open("%s/target.sam" % (out_folder), 'r') as file:
  name=[]
  for line in file.readlines():
    if line.split('\t')[0] not in name:
      name.append(line.split('\t')[0])

header=run("%s view -H %s" % (samtools, bam),  stdout=PIPE, shell=True).stdout.decode('utf-8')

def generate_fastq(item):
  global positive
  global negative

  # Get both sequenses from bam of two reads
  reads1=run("%s view -f 64 %s | grep %s" % (samtools, bam, item),  stdout=PIPE, shell=True).stdout.decode('utf-8')
  if len(reads1)>0:
    read1=reads1.split('\t')[9]
  else:
    read1=''
  reads2=run("%s view -f 128 %s | grep %s" % (samtools, bam, item),  stdout=PIPE, shell=True).stdout.decode('utf-8')
  if len(reads2)>0:
    read2=reads2.split('\t')[9]
  else:
    read2=''

  # Add reads to fastq if sequenses matches positive dict
  if read1 in positive and positive[read1]==read2:
    print(f'index: {name.index(item)}, pool:{mp.current_process()._identity}, positive len: {len(positive)}')
    with open("%s/target_%s.sam" % (out_folder, name.index(item)), "w") as file:
      file.write(header)
    with open("%s/target_%s.sam" % (out_folder, name.index(item)), "a") as file:
      file.write(reads2)
      file.write(reads1)
    run("%s view -b %s/target_%s.sam | %s sort > %s/target_%s.bam" % (samtools,out_folder, name.index(item), samtools,out_folder, name.index(item)), shell=True) 
    run("%s index %s/target_%s.bam" % (samtools,out_folder, name.index(item)), shell=True)
    run("%s fastq -1 %s/tmp_1_%s.fq -2 %s/tmp_2_%s.fq %s/target_%s.bam" % (samtools,out_folder, name.index(item),out_folder, name.index(item),out_folder, name.index(item)), shell=True)
    run("cat %s/tmp_1_%s.fq >> %s/%s:%s%s\">\"%s.positive.R1.fq" % (out_folder, name.index(item),out_folder,chr,start,refl,altl), shell=True)
    run("cat %s/tmp_2_%s.fq >> %s/%s:%s%s\">\"%s.positive.R2.fq" % (out_folder, name.index(item),out_folder,chr,start,refl,altl), shell=True) 
    run("rm %s/target_%s.sam %s/target_%s.bam %s/target_%s.bam.bai %s/tmp_1_%s.fq %s/tmp_2_%s.fq" % (out_folder, name.index(item), out_folder, name.index(item), out_folder, name.index(item), out_folder, name.index(item), out_folder, name.index(item)), shell=True)

  # Skip reads if sequenses matches negative dict
  elif (read1 in negative) and (negative[read1]==read2):
    print(f'index: {name.index(item)}, pool:{mp.current_process()._identity}, negative len:{len(negative)}')

  # Run SGA if reads don't match dict
  else:
    # Generate bam with pair of mates (by Read Name)
    with open("%s/target_%s.sam" % (out_folder, name.index(item)), "w") as file:
      file.write(header)
    with open("%s/target_%s.sam" % (out_folder, name.index(item)), "a") as file:
      file.write(reads2)
      file.write(reads1)
    run("%s view -b %s/target_%s.sam | %s sort > %s/target_%s.bam" % (samtools,out_folder, name.index(item), samtools,out_folder, name.index(item)), shell=True) 
    run("%s index %s/target_%s.bam" % (samtools,out_folder, name.index(item)), shell=True)
    result = run("%s somatic-variant-filters -t 2 --tumor-bam %s/target_%s.bam --normal-bam %s/target_%s.bam --reference %s %s/target.vcf --annotate-only" % (sga,out_folder,name.index(item),out_folder,name.index(item), ref,out_folder),  stdout=PIPE, shell=True)
    AF=int(float(result.stdout.decode('utf-8').split('\t')[7].split('TumorVAF=')[1].split(';')[0]))
    if AF == 1:
      positive[read1]=read2
      print(f'index: {name.index(item)}, pool:{mp.current_process()._identity}, positive len:{len(positive)}')
      run("%s fastq -1 %s/tmp_1_%s.fq -2 %s/tmp_2_%s.fq %s/target_%s.bam" % (samtools,out_folder, name.index(item),out_folder, name.index(item),out_folder, name.index(item)), shell=True)
      run("cat %s/tmp_1_%s.fq >> %s/%s:%s%s\">\"%s.positive.R1.fq" % (out_folder, name.index(item),out_folder,chr,start,refl,altl), shell=True)
      run("cat %s/tmp_2_%s.fq >> %s/%s:%s%s\">\"%s.positive.R2.fq" % (out_folder, name.index(item),out_folder,chr,start,refl,altl), shell=True) 
      run("rm %s/target_%s.sam %s/target_%s.bam %s/target_%s.bam.bai %s/tmp_1_%s.fq %s/tmp_2_%s.fq" % (out_folder, name.index(item), out_folder, name.index(item), out_folder, name.index(item), out_folder, name.index(item), out_folder, name.index(item)), shell=True)
    else:
      negative[read1]=read2
      print(f'index: {name.index(item)}, pool:{mp.current_process()._identity}, negative len:{len(negative)}')
      run("rm %s/target_%s.sam %s/target_%s.bam %s/target_%s.bam.bai" % (out_folder, name.index(item), out_folder, name.index(item), out_folder, name.index(item)), shell=True)

# For each name check if the mate pair is positive for variant

positive={}
negative={}

pool = mp.Pool(10)
pool.map(generate_fastq, [item for item in name])
pool.close()

