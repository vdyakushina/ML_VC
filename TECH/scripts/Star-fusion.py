#!/usr/bin/env python3.8
from subprocess import run
from subprocess import PIPE
import sys
import os

with open('/media/aod/DATA/STAR-fusion/samples_fusions.csv', 'r') as in_f:
	for line in in_f:
		samples.append(line.split('\t')[1])

for sample in samples:
	## create directory for results
	try:
		os.mkdir('samples/22382-01-04')
	except OSError:
		print ("Creation of the directory samples/%s failed" % '22382-01-04')
	else:
		print ("Creation of the directory samples/%s SUCCESS" % '22382-01-04')
	## make fastq from bam
	try:
		bam=run("ls /media/aod/DATA/data/samples/%s/raw/*.bam" % sample, stdout=PIPE, shell=True).stdout.decode('utf-8').strip().split('DATA')[1]
		run("sudo docker run -v /media/aod/DATA/:/home/ trinityctat/starfusion java -jar picard.jar SamToFastq I=/home%s FASTQ=/home/STAR-fusion/samples/%s/%s.fastq" %(bam, sample, sample), shell=True)
		## run STAR-fusion
		run("sudo docker run -v /media/aod/DATA/:/home/ trinityctat/starfusion \
			STAR-Fusion --genome_lib_dir /home/STAR-fusion/STAR-fusion_db/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ \
			--left_fq /home/STAR-fusion/samples/%s/%s.fastq --output_dir /home/STAR-fusion/samples/%s/star_fusion_outdir \
			--no_remove_dups --FusionInspector validate --no_RT_artifact_filter --CPU 14" %(sample, sample, sample), shell=True)
	except:
		with open('/media/aod/DATA/STAR-fusion/samples_truncated_bams.csv', 'a') as trunc_f:
			trunc_f.write(sample+'\n')