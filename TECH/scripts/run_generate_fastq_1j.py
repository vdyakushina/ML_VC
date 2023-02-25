#!/usr/bin/env python3

from subprocess import run
from pathlib import Path
from func import read_resources
import sys
import os.path
import csv
import multiprocessing as mp
import numpy as np

(f, n) = sys.argv[1:3]

res = os.getcwd()+"/resources.csv"
path = read_resources(res)

vars=[]
with open('/home/gkhvorykh/train_positive/File1.tsv', 'r') as inf:
	for line in inf:
		line=line.strip()
		var=line.split('\t')[0]
		vars.append(var)

def run_func(variant):
	run("./generate_fastq_1j.py '%s' %s %s" % (variant, f, n), shell = True)

pool = mp.Pool(10)
pool.map(run_func, [variant for variant in vars])
pool.close()