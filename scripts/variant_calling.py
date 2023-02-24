#!/usr/bin/env python3
# Wrapper script to process variant calling

from subprocess import run
from pathlib import Path
from func import read_resources
import sys, glob, os.path

### Initiate
(out, target) = sys.argv[1:4]
res = os.getcwd()+"/resources.csv"
path = read_resources(res)
pref = out + "/" + 'provisional'

if len(glob.glob('%s/*variant*vcf*' % out))>=1:
	run('rm %s/*variant*vcf*' % out, shell=True)
if len(glob.glob('%s/variant.raw.*.tsv' % out))>=1:
	run('rm %s/variant.raw.*.tsv' % out, shell=True)

### Call variants
print('\n### Call variants')

callers=['vardict', 'sinvict']
for caller in callers:
	run("sh %s/%s.sh %s %s %s %s" % (path['scripts'], caller, pref + ".bam", out + "/", path['ref'], target), shell = True)
	run("%s/%s_process.py %s " % (path['scripts'], caller, out + "/"), shell = True)

callers=['scalpel_to', 'scalpel_tnp', 'strelka_tnp', 'strelka_to', 'mutect2_tnp', 'mutect2_to']
for caller in callers:
	run("%s/%s.py %s %s %s %s" % (path['scripts'], caller, pref + ".bam", out + "/", path['ref'], target), shell = True)
	run("%s/%s_process.py %s " % (path['scripts'], caller, out + "/"), shell = True)

### Add VAF, DP and SB from sga to variant.raw
print('\n### Add VAF and SB from sga to variant.raw')
run("%s/%s.py %s %s" % (path['scripts'], 'sga', pref + ".bam", out + "/"), shell = True)

### Convert variant.raw*.vcf.gz to variant.raw*.tsv
print('\n### Convert variant.raw*.vcf.gz to variant.raw*.tsv')
run("%s/%s.py %s" % (path['scripts'], 'vcf_to_tsv', out + "/" + "variant.raw.INDEL.vcf.gz"), shell = True)
run("%s/%s.py %s" % (path['scripts'], 'vcf_to_tsv', out + "/" + "variant.raw.SNV.vcf.gz"), shell = True)

### Run mpileup
run("samtools mpileup -d 0 --reference %s %s > %s" % (path['ref'], pref + ".bam", out + "/variant.mpileup.vcf"), shell = True)
run("perl /home/onco-admin/RnD/UEBAcall/mpileupToVcf.pl -input %s -output %s -limit 0" % (out + "/variant.mpileup.vcf", out + "/variant.mpileup_processed.vcf"), shell = True)
run("bcftools norm -m -any -f %s %s | bgzip > %s" % (path['ref'], out + "/variant.mpileup_processed.vcf", out + "/variant.mpileup.vcf.gz"), shell = True)
run("tabix -f %s" % out + "/variant.mpileup.vcf.gz", shell = True)