#!/usr/bin/env python3
# Preprocess data

from subprocess import run
from func import read_resources
import sys
import os

# Initiate
(bam, bed, pref) = sys.argv[1:4]
res = "/pipeline/scripts/resources.csv"
path = read_resources(res)

# Trimm Illumina adapters
# java -jar /pipeline/tools/picard.jar SamToFastq I=/pipeline/data/samples_test/7IoPucH4PKOgl/raw/7IoPucH4PKOgl.bam FASTQ=prov_R1.fq SECOND_END_FASTQ=prov_R2.fq
# java -jar ~/VYakushina/20210610_155149/Trimmomatic-0.39/trimmomatic-0.39.jar PE ./prov_R1.fq ./prov_R2.fq ./prov_adaptertrimmed_R1.fq ./prov_adaptertrimmed_R1_unpaired.fq ./prov_adaptertrimmed_R2.fq ./prov_adaptertrimmed_R2_unpaired.fq ILLUMINACLIP:/home/immichail/VYakushina/Illumina.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
# bwa mem -M /pipeline/data/reference/hg19.fa prov_adaptertrimmed_R1.fq prov_adaptertrimmed_R2.fq | samtools sort -O BAM - -o adapter_trimmed.bam

# Subset reads by regions
out = pref + ".t1.bam"
run("%s view -b -L %s -o %s %s" % (path['samtools'], bed, out, bam), shell = True)

# Make group of reads
inp = out
out = pref + ".t2.bam"
run("java -jar %s AddOrReplaceReadGroups I=%s O=%s RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20" %
		(path['picard'], inp, out), shell = True)

# Sort and index
inp = out
out = pref + ".target.group.bam"
run("%s sort %s -o %s" % (path['samtools'], inp, out), shell = True)
run("%s index %s" % (path['samtools'], out), shell = True)

# Clean
for p in [".t1.bam", ".t2.bam"]:
	x = pref + p
	try:
		os.remove(x)
	except:
		print("Error while deleting ", x)

# Run GATK Indel Realignment
inp = out
run("python3 %s/realign.py %s %s" % (path['scripts'], inp, pref), shell = True)

# Run GATK Base quality recalibration
inp = pref + ".realigned.bam"
run("python3 %s/bqsr.py %s %s" % (path['scripts'], inp, pref), shell = True)

# Formulate final working .bam file

os.rename(pref + ".bqsr.bam", pref + ".bam");
os.rename(pref + ".bqsr.bai", pref + ".bam.bai");
