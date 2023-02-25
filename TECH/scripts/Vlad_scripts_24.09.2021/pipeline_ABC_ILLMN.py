#!/usr/bin/env python3
# Wrapper script to process the whole pipeline

from subprocess import run
from pathlib import Path
from func import read_resources
import sys
import os.path

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

### Check number of arguments
if len(sys.argv) != 4:
	print("Wrong number of arguments!")
	help()

### Initiate
(bam, target, out) = sys.argv[1:4]
# out = "analysis"
res = "/pipeline/scripts/resources.csv"
path = read_resources(res)

### Check input
bed = path.get(target)
if not bed:
	print(f'{target} is unknown target id')
	sys.exit(1)

for f in [bam, bed]:
	if not Path(f).is_file():
		print(f'{f} doesn\'t exist')
		sys.exit(1)


### Get sample id
sample = Path(bam).stem

### Show input
print(f'\nbam: {bam}\nTarget: {target}\nSample: {sample}\nbed: {bed}\nOutput: {out}\n')

### Create output folder if it doesn't exist
Path(out).mkdir(exist_ok = True)


### Preprocess data
print('\n### Preprocess data')
pref = "'" + out + "/" + 'provisional' + "'"
print('passed bam:', bam)
print('pref:', pref)
print("python3 %s/preprocess.py '%s' %s %s" % (path['scripts'], bam, bed, pref))
run("python3 %s/preprocess.py '%s' %s %s" % (path['scripts'], bam, bed, pref), shell = True)


# pref = out + "/" + 'provisional'

### Call variants
print('\n### Call variants')
callers=['vardict', 'scalpel', 'sinvict']
for caller in callers:
        run("sh %s/%s.sh %s %s %s %s" % (path['scripts'], caller, pref + ".bam", out + "/", path['ref'], target), shell = True)
        run("%s/%s_process.py %s " % (path['scripts'], caller, out + "/"), shell = True)

callers=['strelka_tnp', 'strelka_to', 'mutect2_tnp', 'mutect2_to']
for caller in callers:
        run("%s/%s.py %s %s %s %s" % (path['scripts'], caller, pref + ".bam", out + "/", path['ref'], target), shell = True)
        run("%s/%s_process.py %s " % (path['scripts'], caller, out + "/"), shell = True)


### Add VAF and SB from sga to variant.raw
print('\n### Add VAF and SB from sga to variant.raw')
run("%s/%s.py %s %s" % (path['scripts'], 'sga', pref + ".bam", out + "/"), shell = True)

'''
### Remove variants found in ExAC
print('\n### Remove variants found in ExAC')
run("sh %s/%s.sh %s" % (path['scripts'], "exac_subtraction_from_raw", out + "/" + "variant.raw.INDEL.vcf.gz"), shell=True)
run("sh %s/%s.sh %s" % (path['scripts'], "exac_subtraction_from_raw", out + "/" + "variant.raw.SNV.vcf.gz"), shell=True)
'''

### Convert variant.raw*.vcf.gz to variant.raw*.tsv
print('\n### Convert variant.raw*.vcf.gz to variant.raw*.tsv')
run("python %s/%s.py %s" % (path['scripts'], 'vcf_to_tsv', out + "/" + "variant.raw.INDEL.vcf.gz"), shell = True)
run("python %s/%s.py %s" % (path['scripts'], 'vcf_to_tsv', out + "/" + "variant.raw.SNV.vcf.gz"), shell = True)

