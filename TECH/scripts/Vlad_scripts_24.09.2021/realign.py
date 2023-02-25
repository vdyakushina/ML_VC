#!/usr/bin/env python3
# Realign locally around the indels
# Source: https://qcb.ucla.edu/wp-content/uploads/sites/14/2016/03/GATKwr12-3-IndelRealignment.pdf

import sys
import os
from subprocess import run
from func import read_resources

# Initiate
(bam, pref) = sys.argv[1:3]
res = '/pipeline/scripts/resources.csv'

# Check input
if not os.path.isfile(bam):
	print(f'{bam} doesn\'t exist!')
	sys.exit(1)

path = read_resources(res)

# Create target intervals list 
out = pref + ".realigner.intervals"
out = "'" + out + "'"
bam = "'" + bam + "'"
print("\nTARGET INTERVALS\n%s -jar %s -T RealignerTargetCreator -R %s -I %s -o %s" % (path['java8'], path['gatk3.8'], path['ref'], bam, out))
run("%s -jar %s -T RealignerTargetCreator -R %s -I %s -o %s" %
(path['java8'], path['gatk3.8'], path['ref'], bam, out), shell = True)

# Realign reads
inp = out  # [1:-1]
arr_pref = '\>'.join(pref.split('>'))
print('arr_pref', arr_pref)
out = "'" + arr_pref + ".realigned.bam" + "'"
print("%s -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s --maxConsensuses 990000 --maxReadsForConsensuses 990000 --maxReadsForRealignment 990000 --maxReadsInMemory 990000" % (path['java8'], path['gatk3.8'], path['ref'], bam, inp, out))
run("%s -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s --maxConsensuses 990000 --maxReadsForConsensuses 990000 --maxReadsForRealignment 990000 --maxReadsInMemory 990000" %
(path['java8'], path['gatk3.8'], path['ref'], bam, inp, out), shell = True)

# Clean
try:
	os.remove(inp)
except:
	print("Error deleting ", inp)

