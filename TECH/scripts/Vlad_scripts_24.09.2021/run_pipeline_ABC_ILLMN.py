#!/usr/bin/env python3
# Wrapper script to process the whole pipeline

from subprocess import run
from pathlib import Path
from func import read_resources
import sys
import os.path
import csv
import multiprocessing as mp

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
target = 'AODABCV1'

samples=['chr11:108121509A>G.1', 'chr11:108121509A>G.2', 'chr11:108121509A>G.3', 'chr11:108121509A>G.4', 'chr11:108121509A>G.5', 'chr11:108168115CTA>C.1', 'chr11:108168115CTA>C.2', 'chr11:108168115CTA>C.3', 'chr11:108168115CTA>C.4', 'chr11:108168115CTA>C.5', 'chr13:32912337CTG>C.1', 'chr13:32912337CTG>C.2', 'chr13:32912337CTG>C.3', 'chr13:32912337CTG>C.4', 'chr13:32912337CTG>C.5']
samples=['chr11:108121509A>G.1']
#with open(path['meta'], newline='') as csvfile:
#	spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
#	for row in spamreader:
#		samples.append(row[0])

def run_pipeline(sample):
	print('sample', sample)
#	bam="/pipeline/data/samples_xls/%s/analysis/provisional.bam" % sample
	bam="/pipeline/data/samples_xls/%s/raw/%s.sequences.sorted.bam" % (sample, sample)
	out="/pipeline/data/samples_xls/%s/analysis" % sample

	print('\necho has been started')
	print("echo \'sample: %s, bam: %s, out: %s\' >> /home/gkhvorykh/samples/samples_in_process.vlad_positive.txt" % (sample, bam, out))
	run("echo \'sample: %s, bam: %s, out: %s\' >> /home/gkhvorykh/samples/samples_in_process.vlad_positive.txt" % (sample, bam, out), shell = True)
	
	print('\npipeline has been started')
	print('%s/pipeline_ABC_ILLMN.py "%s" "%s" "%s"' % (path['scripts'], bam, target, out + "/"))
	run('%s/pipeline_ABC_ILLMN.py "%s" "%s" "%s"' % (path['scripts'], bam, target, out + "/"), shell = True)

pool = mp.Pool(10)
pool.map(run_pipeline, [sample for sample in samples])
pool.close()
