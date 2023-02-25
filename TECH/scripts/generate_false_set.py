#!/usr/bin/env python3


from subprocess import run, PIPE
from pathlib import Path
from func import read_resources
import sys
import os.path

for i in ['SNV', 'INDEL']:
	samples={}
	with open('/home/gkhvorykh/samples/false_variants_%s.tsv' % i, 'r') as in_s:
		for line in in_s:
			line=line.strip()
			vart=line.split('\t')[0]
			sample=line.split('\t')[1]
			if sample in samples:
				samples[sample].append(vart)
			else:
				samples[sample]=[vart]
	for sample in samples.keys():
		for var in samples[sample]:
			result=run('grep \'%s\' /home/gkhvorykh/samples/%s/analysis/variant.raw.%s.tsv' %(var, sample, i), shell = True, stdout=PIPE).stdout.decode()
			if not os.path.exists('/home/gkhvorykh/samples/train_false_results_%s.tsv' % i):
				with open('/home/gkhvorykh/samples/train_false_results_%s.tsv' % i, 'w') as wf:
					clms=['sample', 'variant_name', 'Mutect2_TO_FILTER', 'Mutect2_TO_STRQ', 'Mutect2_TO_TLOD', 'Mutect2_TNP_FILTER', 'Strelka_TO_FILTER', 'Strelka_TO_QUAL', 'Strelka_TNP_FILTER' , 'SC_TO_FILTER', 'SC_TNP_FILTER', 'SINVICT', 'VD_FILTER', 'VD_Q', 'SGA_VAF', 'SGA_SB', 'SGA_RepeatRefCount', 'SGA_DP']
					if i == 'SNV':
						for x in ['SC_TO_FILTER', 'SC_TNP_FILTER']:
							clms.remove(x)
					wf.write('\t'.join(clms)+'\n')
			with open('/home/gkhvorykh/samples/train_false_results_%s.tsv' % i, 'a') as outf:
				outf.write(sample+'\t'+result)