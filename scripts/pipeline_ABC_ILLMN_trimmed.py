#!/usr/bin/env python3
# Wrapper script to process the whole pipeline

from subprocess import run, PIPE
from pathlib import Path
from func import read_resources
import sys, time
import os.path
import multiprocessing as mp
import pandas as pd


def help():
	print('''
	Script runs all steps of the pipeline.
	The output files are placed in the '/analysis' subfolder,
	which is created in the current folder.

	Usage: pipeline_ABC_ILLMN.py <bam> target
		bam: path/to/file.bam with aligned reads
		target: target id
	''')
	sys.exit(1)

### Initiate
samples={}
if '--single' in sys.argv:
	(fastq1, fastq2, dir, panel)=sys.argv[2:6]
	name=os.path.commonprefix([os.path.basename(fastq1), os.path.basename(fastq2)]).strip('_').replace('_R', '')
	samples[name]={'fastq1':fastq1, 'fastq2':fastq2, 'panel':panel 'out':dir}

else:
	(run_info) = sys.argv[1]
	with open(run_info, 'r') as in_f:
		for line in in_f:
			line=line.strip()
			samples[line.split('\t')[0]]={}
			(samples[line.split('\t')[0]]['fastq1'], samples[line.split('\t')[0]]['fastq2'], samples[line.split('\t')[0]]['panel'], samples[line.split('\t')[0]]['out'])=line.split('\t')[1:5]



res = os.getcwd()+"/resources.csv"
path = read_resources(res)

### Preprocess data
def preprocess(sample):
	### Create output folder if it doesn't exist
	Path(samples[sample]['out']).mkdir(exist_ok = True)
	start_time=time.perf_counter()
	run("%s/preprocess.py %s %s %s %s" % (path['scripts'], samples[sample]['fastq1'], samples[sample]['fastq2'], samples[sample]['out'], samples[sample]['panel']), shell = True)
	total_time=round((time.perf_counter()-start_time)/60, 2)
	with open('%s/process.txt' % samples[sample]['out'], 'a') as ready_list:
		ready_list.write(f'preprocessing\tsample:{sample}\tminutes:{total_time}\n')

pool = mp.Pool(5)
pool.map(preprocess, [str(sample) for sample in samples.keys()])
pool.close()

### QC
def qc(sample):
	panel=samples[sample]['panel']
	if panel=='unknown':
		panel=open('%s/panel_name.csv' % out).read().split('\t')[0]
	start_time=time.perf_counter()
	run("%s/qc.py %s %s %s" % (path['scripts'], samples[sample]['out']+'/'+'provisional.bam', 'provisional', panel), shell = True)
	total_time=round((time.perf_counter()-start_time)/60, 2)
	with open('%s/process.txt' % samples[sample]['out'], 'a') as ready_list:
		ready_list.write(f'qc\tsample:{sample}\tminutes:{total_time}\n')

pool = mp.Pool(5)
pool.map(qc, [str(sample) for sample in samples.keys()])
pool.close()

### Call variants
def variant_calling(sample):
	panel=samples[sample]['panel']
	if panel=='unknown':
		panel=open('%s/panel_name.csv' % samples[sample]['out']).read().split('\t')[0]
	start_time=time.perf_counter()
	run("%s/variant_calling.py %s %s" % (path['scripts'], samples[sample]['out'], panel), shell = True)
	total_time=round((time.perf_counter()-start_time)/60, 2)
	with open('%s/process.txt' % samples[sample]['out'], 'a') as ready_list:
		ready_list.write(f'variant_calling\tsample:{sample}\tminutes:{total_time}\n')

pool = mp.Pool(5)
pool.map(variant_calling, [str(sample) for sample in samples])
pool.close()

### ML
def ML(sample):
	run("%s/RandomForest.py %s %s" % (path['scripts'], samples[sample]['out']+'/variant.raw.SNV.tsv', 'SNV'), shell = True)
	run("%s/RandomForest.py %s %s" % (path['scripts'], samples[sample]['out']+'/variant.raw.INDEL.tsv', 'INDEL'), shell = True)

pool = mp.Pool(10)
pool.map(ML, [str(sample) for sample in samples])
pool.close()

### Add annotation
def annotation(sample):
	start_time=time.perf_counter()
	for vtype in ['SNV', 'INDEL']:
		run("%s/interpretation.py %s/variant.raw.%s.vcf.gz annotation_%s.tsv" % (path['scripts'], samples[sample]['out'], vtype, vtype), shell = True)
		calling=pd.read_csv('%s/variant.raw.%s.tsv' % (samples[sample]['out'], vtype), sep='\t')
		annotation=pd.read_csv('%s/annotation_%s.tsv' % (samples[sample]['out'], vtype), sep='\t')
		calling=calling.merge(right=annotation, left_on='variant_name', right_on='name', how='left')
		calling.drop("name", 1).to_csv('%s/variant.raw.%s.tsv' % (samples[sample]['out'], vtype), sep='\t', index=False)
	total_time=round((time.perf_counter()-start_time)/60, 2)
	with open('%s/samples_finished_time.txt' % runpath, 'a') as ready_list:
		ready_list.write(f'annotation\tsample:{sample}\tminutes:{total_time}\n')

pool = mp.Pool(10)
pool.map(annotation, [str(sample) for sample in samples])
pool.close()

if '--single' in sys.argv:
	print(f'"SINGLE" module is not sutible for further varint classification with UEBA')
	sys.exit()
'''
### Prepare gathered vcf
runpath=os.path.commonprefix([i['out'] for i in samples.values()])
run("bcftools concat -a -d all --no-version %s/*/analysis/variant.raw.*.vcf.gz -O v -o %s/gathered.vcf" % (runpath, runpath), shell = True)

### Run UEBA
if panel=='unknown':
	panel=open('%s/analysis/panel_name.csv' % samples[0]).read().split('\t')[0]
bed=path[panel]
UEBA_config=path[panel+'_UEBA_config']
#run('perl /home/onco-admin//RnD/UEBAcall/make_distributions.pl -l %s/list_bam -v %s/gathered.vcf -p %s -bdata %s/bdata -mode append -n 13' % 
#(UEBA_config, runpath, bed, UEBA_config), shell=True)

def ueba(sample):
	out=samples[sample]['out']
	panel=samples[sample]['panel']
	if panel=='unknown':
		panel=open('%s/panel_name.csv' % out).read().split('\t')[0]
	start_time=time.perf_counter()
	run("%s/UEBA.py %s %s" % (path['scripts'], out, panel), shell = True)
	run("%s/UEBA_process.py %s" % (path['scripts'], out), shell = True)
	total_time=round((time.perf_counter()-start_time)/60, 2)
	with open('%s/samples_finished_time.txt' % runpath, 'a') as ready_list:
		ready_list.write(f'UEBA\tsample:{sample}\tminutes:{total_time}\n')

for sample in samples:
	ueba(sample)
'''