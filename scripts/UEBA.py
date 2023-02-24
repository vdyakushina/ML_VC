#!/usr/bin/env python3
# Run mutect2 in tumor only mode 
import sys
import os
from subprocess import run
from func import read_resources

### Initiate
(folder, target) = sys.argv[1:5]

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
bed=path[target]
UEBA_config=path[target+'_UEBA_config']

for type in ['SNV', 'INDEL']:
	run('gunzip -k -f %s/variant.raw.%s.vcf.gz' % (folder, type), shell=True)
	
	run('perl /home/onco-admin//RnD/UEBAcall/get_counts_for_sample.pl -s %s/provisional.bam -v %s/variant.raw.%s.vcf -p %s -n 6 -o %s/UEBA.cdata' %
	(folder, folder, type, bed, folder), shell=True)
	
	run('perl /home/onco-admin//RnD/UEBAcall/make_call.pl -v %s/variant.raw.%s.vcf -bdata %s/bdata -cdata %s/UEBA.cdata -n 6 > %s/UEBA.calls' %
	(folder, type, UEBA_config, folder, folder), shell=True)
	
	run('perl /home/onco-admin//RnD/UEBAcall/makeVCF.pl -input %s/UEBA.calls -output %s/variant.UEBA.%s.vcf -sample %s/provisional.bam' %
	(folder, folder, type, folder), shell=True)
	
	run('perl /home/onco-admin//RnD/UEBAcall/filterVCF.pl -input %s/variant.UEBA.%s.vcf' % (folder, type), shell=True) 

	run('bgzip -f %s/variant.UEBA.%s.vcf' % (folder, type), shell=True)
	run('tabix -f %s/variant.UEBA.%s.vcf.gz' % (folder, type), shell=True)
	
	run('rm %s/variant.raw.%s.vcf %s/UEBA.cdata %s/UEBA.calls' % (folder, type, folder, folder), shell=True)