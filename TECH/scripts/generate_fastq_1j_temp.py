#!/usr/bin/env python3.8
from subprocess import run
from subprocess import PIPE
import re
import sys
import os
from func import read_resources
import multiprocessing as mp
import numpy as np

(variant, f, DP) = sys.argv[1:4]

res = os.getcwd()+"/resources.csv"
path = read_resources(res)

##Read sequences and quals from positive and negative fastq files
seqs_p_1={}
seqs_p_2={}
seqs_n_1={}
seqs_n_2={}

nn=0
seqn=1
qn=3
with open('/home/gkhvorykh/train_positive/positive_fastq/%s.positive.R1.fq' % variant, 'r') as inpf:
	for i, line in enumerate(inpf):
		if i==nn:
			name=line.strip()
			nn+=4
		elif i == seqn:
			seqs_p_1[name]=[line.strip()]
			seqn+=4
		elif i == qn:
			seqs_p_1[name].append(line.strip())
			qn+=4

nn=0
seqn=1
qn=3
with open('/home/gkhvorykh/train_positive/positive_fastq/%s.positive.R2.fq' % variant, 'r') as inpf:
	for i, line in enumerate(inpf):
		if i==nn:
			name=line.strip()
			nn+=4
		elif i == seqn:
			seqs_p_2[name]=[line.strip()]
			seqn+=4
		elif i == qn:
			seqs_p_2[name].append(line.strip())
			qn+=4


nn=0
seqn=1
qn=3
with open('/home/gkhvorykh/train_negative/negative_fastq/%s.negative.R1.fq' % variant, 'r') as inpf:
	for i, line in enumerate(inpf):
		if i==nn:
			name=line.strip()
			nn+=4
		elif i == seqn:
			seqs_n_1[name]=[line.strip()]
			seqn+=4
		elif i == qn:
			seqs_n_1[name].append(line.strip())
			qn+=4

nn=0
seqn=1
qn=3
with open('/home/gkhvorykh/train_negative/negative_fastq/%s.negative.R2.fq' % variant, 'r') as inpf:
	for i, line in enumerate(inpf):
		if i==nn:
			name=line.strip()
			nn+=4
		elif i == seqn:
			seqs_n_2[name]=[line.strip()]
			seqn+=4
		elif i == qn:
			seqs_n_2[name].append(line.strip())
			qn+=4


def run_func(DP):
	##Prepare directories where fastq, raw bams and analysis will be stored
	folder=variant.replace('>','_').replace(':','_')+'_'+str(f)+'_'+str(DP)
	if not os.path.exists('/home/gkhvorykh/train_true_samples/%s/raw/' % folder):
		run('mkdir -p \'/home/gkhvorykh/train_true_samples/%s/raw/\'' % folder, shell=True)
	if not os.path.exists('/home/gkhvorykh/train_true_samples/%s/analysis/' % folder):
		run('mkdir -p \'/home/gkhvorykh/train_true_samples/%s/analysis/\'' % folder, shell=True)

	##Define paths for fastq and bam
	out_fq1='/home/gkhvorykh/train_true_samples/%s/raw/%s_R1.fq' % (folder,folder)
	out_fq2='/home/gkhvorykh/train_true_samples/%s/raw/%s_R2.fq' % (folder,folder)
	out_bam='/home/gkhvorykh/train_true_samples/%s/raw/%s.bam' % (folder,folder)
	##Generate resulting fastq R1 & R2
	run('rm %s %s' % (out_fq1, out_fq2), shell=True)
	i=1
	while i <= int(DP):
		F=np.random.randint(1,10000)
		if F < 10000*float(f):
			k=0
			while k < 1:
				num=np.random.randint(1,len(seqs_p_1))
				name=list(seqs_p_1.keys())[num]
				name_rand=name+':'+str(np.random.randint(1000000,9999999))
				if name in seqs_p_2.keys():
					with open(out_fq1, 'a') as outf1, open(out_fq2, 'a') as outf2:
						outf1.write(name_rand+'\n'+'\n+\n'.join(seqs_p_1[name])+'\n')
						outf2.write(name_rand+'\n'+'\n+\n'.join(seqs_p_2[name])+'\n')
					k+=1
				else:
					continue
		else:
			c=0
			while c < 1:
				num=np.random.randint(1,len(seqs_n_1))
				name=list(seqs_n_1.keys())[num]
				name_rand=name+':'+str(np.random.randint(1000000,9999999))
				if name in seqs_n_2.keys():
					with open(out_fq1, 'a') as outf1, open(out_fq2, 'a') as outf2:
						outf1.write(name_rand+'\n'+'\n+\n'.join(seqs_n_1[name])+'\n')
						outf2.write(name_rand+'\n'+'\n+\n'.join(seqs_n_2[name])+'\n')
					c+=1
				else:
					continue
		i+=1
	run('%s mem -M %s %s %s | samtools sort -O BAM - -o %s' % \
	(path['bwa'], path['ref'], out_fq1, out_fq2, out_bam), shell = True)

	run('%s index %s' % (path['samtools'], out_bam) , shell = True)

	run('rm %s %s' % (out_fq1, out_fq2), shell=True)
	with open('/home/gkhvorykh/train_true_samples/meta.tsv', 'a') as fout:
		fout.write(folder+'\n')

run_func(DP)