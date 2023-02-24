#!/usr/bin/env python3
# Estimate the quality of sequensing data. 
# The output is written to 'qc/' folder which is created 
# in the current working folder. 

import argparse, sys, os, subprocess, json
from subprocess import run, PIPE
from pathlib import Path
from ast import literal_eval
from func import read_resources
from func import ts

# Initiate the variables

resources = os.getcwd()+"/resources.csv"

# Get command line arguments
parser = argparse.ArgumentParser(description = 'Estimates QC metrics for coverage.')
parser.add_argument('bam', help = 'path/to/sample.bam file')
parser.add_argument('sample', help = 'sample id')
parser.add_argument('target', help = 'panel id')

args = parser.parse_args()
bam = Path(args.bam).resolve()
sample = args.sample
target = args.target

pools = [target+'_'+str(i) for i in [1,2]]
genes = ["BRCA"]


# Get paths to resources
path = read_resources(resources)

# Check imput
if not os.path.isfile(bam): sys.exit(f'{bam} doesn\'t exist!')
if not os.path.isfile(path[target]): sys.exit(f'{path[target]} doesn\'t exist!')

# Determine target
if panel=='unknown':
	rbam=pref+'.bam'
	target=run('%s/panel_search.py %s "AODHRD15,AODCPV1,AODABCV1"'% (path['scripts'], rbam),  stdout=PIPE, shell=True).stdout.decode('utf-8').strip()
	open('%s/panel_name.csv' % out_f, 'w').write(target)
	panel = target.split('\t')[0]
bed = path.get(panel)

# Initiate the dictionary for QC metrics
qc = {'sample': sample}

# === functions ===

def cov(bam, target, sample):
	res = run("python3 %s/depthnuc.py %s %s %s" % (path['scripts'], bam, target, sample),stdout=PIPE, shell=True).stdout.decode('utf-8').strip()
	return res

# === main ===

# Get working directory
wd = os.getcwd()
# Create the qc folder
Path(os.path.dirname(bam)+'/qc').mkdir(parents = True, exist_ok = True)
# Step into qc folder
os.chdir(os.path.dirname(bam)+'/qc')

# Estimate the coverage and metrics
print(ts(), "Estimate coverage")
res = cov(bam, target, sample)
qc.update(literal_eval(res))

# Estimate on-target 
ontarget = run("python3 %s/ontarget.py %s %s" % (path['scripts'], bam, target), stdout=PIPE, shell=True).stdout.decode('utf-8').strip()
qc.update(literal_eval(ontarget))

# Split into pools and estimate the coverage and qc metrics for each pool
qc.update({'pool':[]})
k = 1
for pool in pools:
	pref="%s.pool%s" % (sample, k)
	subprocess.check_call("%s/subset_reads.sh %s %s %s" % (path['scripts'], bam, pool, pref),shell = True)
	# Estimate the coverage and metrics
	res = cov(pref+".bam", pool, pref)
	t = {'number' : k}
	t.update(literal_eval(res))
	qc['pool'].append(t)
	k += 1

# Estimate the sensitivity
print(ts(), "Estimate sensitivity")
subprocess.check_call("%s/run_ephagen.sh %s %s" % (path['scripts'], bam, sample), shell = True)

# Parser outputs from EphaGen
for gene in genes:
	with open("%s.%s.tsv" % (sample, gene), 'r') as f:
		l = f.readlines()[1]
		l = l.strip()
		d = l.split()
		qc.update({ gene: [] })
		qc[gene].append({
			'sensitivity': float(d[2])
		})

# Save quality metrics as json
with open(sample + ".json", 'w') as f:
	json.dump(qc, f, indent = 2)

# Come back into working directory
os.chdir(wd)
