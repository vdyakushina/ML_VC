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

print('\npreprocess bam', bam)
print('preprocess pref', pref)
qbam = "'" + bam + "'"
qpref = "'" + pref + "'"
print('\nquoted preprocess bam', qbam)
print('quoted preprocess pref', qpref)


# Trimm Illumina adapters
# java -jar /pipeline/tools/picard.jar SamToFastq I=/pipeline/data/samples_test/7IoPucH4PKOgl/raw/7IoPucH4PKOgl.bam FASTQ=prov_R1.fq SECOND_END_FASTQ=prov_R2.fq
# java -jar ~/VYakushina/20210610_155149/Trimmomatic-0.39/trimmomatic-0.39.jar PE ./prov_R1.fq ./prov_R2.fq ./prov_adaptertrimmed_R1.fq ./prov_adaptertrimmed_R1_unpaired.fq ./prov_adaptertrimmed_R2.fq ./prov_adaptertrimmed_R2_unpaired.fq ILLUMINACLIP:/home/immichail/VYakushina/Illumina.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
# bwa mem -M /pipeline/data/reference/hg19.fa prov_adaptertrimmed_R1.fq prov_adaptertrimmed_R2.fq | samtools sort -O BAM - -o adapter_trimmed.bam

# Subset reads by regions
out = "'" + pref + ".t1.bam" + "'"
print("\nSUBSETTING\n%s view -b -L %s -o %s %s" % (path['samtools'], bed, out, qbam))
run("%s view -b -L %s -o %s %s" % (path['samtools'], bed, out, qbam), shell = True)

# Make group of reads
inp = out
out = "'" + pref + '.t2.bam' + "'"
print("\nREAD GROUPS\njava -jar %s AddOrReplaceReadGroups I=%s O=%s RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20" % (path['picard'], inp, out))
run("java -jar %s AddOrReplaceReadGroups I=%s O=%s RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20" %
		(path['picard'], inp, out), shell = True)

# Sort and index
inp = out
print('sort & index inp', inp)
out = "'" + pref + ".target.group.bam" + "'"
print('sort & index out', out)
print("%s sort %s -o %s" % (path['samtools'], inp, out))
run("%s sort %s -o %s" % (path['samtools'], inp, out), shell = True)
print("%s index %s" % (path['samtools'], out))
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
print("\nINDEL REALIGNMENT\npython3 %s/realign.py %s %s" % (path['scripts'], inp, qpref))
run("python3 %s/realign.py %s %s" % (path['scripts'], inp, qpref), shell = True)

# Run GATK Base quality recalibration
inp = "'" + pref + ".realigned.bam" + "'"
print("\nBQSR\npython3 %s/bqsr.py %s %s" % (path['scripts'], inp, qpref))
run("python3 %s/bqsr.py %s %s" % (path['scripts'], inp, qpref), shell = True)

# Formulate final working .bam file

os.rename(pref + ".bqsr.bam", pref + ".bam");
os.rename(pref + ".bqsr.bai", pref + ".bam.bai");
