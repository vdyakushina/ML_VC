#!/usr/bin/env python3
# Wrapper script to process the whole pipeline

from subprocess import run
from pathlib import Path
from func import read_resources
import sys, csv, time
import os.path
import multiprocessing as mp

res = os.getcwd()+"/resources.csv"
path = read_resources(res)

runpath='/home/onco-admin/RnD/Yakushina/MIPT20220429/'
samples=[ name for name in os.listdir(runpath) if os.path.isdir('%s/%s' % (runpath,name))]
samples=[runpath+item for item in samples if item not in ['Stats', 'Undetermined_S0', 'Reports']]
target='AODHRD15'

def run_pipeline(sample):
	print('sample', sample)
	fastq_dir=sample+'/Fastq'
	out=sample+'/analysis'
	with open('%s/samples_in_process.txt' % runpath, 'a') as in_process_list:
		in_process_list.write(f'sample:{sample}, fastq_dir:{fastq_dir}, out:{out}\n')
	start_time=time.perf_counter()
	run('%s/pipeline_ABC_ILLMN_trimmed.py \'%s\' %s \'%s\'' % (path['scripts'], fastq_dir, target, out), shell = True)
	total_time=round((time.perf_counter()-start_time)/60, 2)
	with open('%s/samples_finished_time.txt' % runpath, 'a') as ready_list:
		ready_list.write(f'sample:{sample}, minutes:{total_time}\n')

pool = mp.Pool(10)
pool.map(run_pipeline, [str(sample) for sample in samples])
pool.close()
