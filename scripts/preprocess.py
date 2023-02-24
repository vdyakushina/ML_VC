#!/usr/bin/env python3
# Preprocess data

from subprocess import run, PIPE
from func import read_resources
import sys, os, re, glob

### Initiate
(fastq1, fastq2, out_folder, panel) = sys.argv[1:5]
res = os.getcwd()+"/resources.csv"
path = read_resources(res)
bwa='bwa'
pref = out_folder + "/" + 'provisional'
name=os.path.commonprefix([os.path.basename(fastq1), os.path.basename(fastq2)]).strip('_').replace('_R', '')

try:
	os.mkdir(out_folder)
except FileExistsError: 
	print(f'{out_folder} exists')

### Remove files resalting from previous preprocessing
if len(glob.glob('%s/*.bam*' % out_folder))>=1:
	run('rm %s/*.bam*' % out_folder, shell=True)
if os.path.isdir('%sbamreadcount' % out_folder):
	run('rm -rf %sbamreadcount' % out_folder, shell=True)

### Trimm Illumina adapters 
print('trimmomatic')
run('java -jar %s/trimmomatic-0.39.jar PE %s %s %s%s_adaptertrimmed_R1.fq %s%s_adaptertrimmed_R1_unpaired.fq %s%s_adaptertrimmed_R2.fq %s%s_adaptertrimmed_R2_unpaired.fq \
	ILLUMINACLIP:%s/adapters/AODABCV1_Illumina.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30' % \
	(path['trimmomatic'], fastq1, fastq2, out_folder, name, out_folder, name, out_folder, name, out_folder, name, path['trimmomatic']), shell = True)

### Align 
print('bwa')
run('%s mem -M %s %s%s_adaptertrimmed_R1.fq %s%s_adaptertrimmed_R2.fq | samtools sort -O BAM - -o %s/%s.bam' % \
	(bwa, path['ref'], out_folder, name, out_folder, name, out_folder, name), shell = True)

run('%s index %s/%s.bam' %('samtools', out_folder, name), shell = True)

### Remove short reads (M in SIGAR < 40%)
print('remove short reads')
run("%s view -H %s/%s.bam > %s.sam" % ('samtools', out_folder, name, pref), shell=True)
reads=run('samtools view %s/%s.bam'% (out_folder, name),  stdout=PIPE, shell=True).stdout.decode('utf-8').split('\n')
with open('%s.sam' % pref, 'a') as out_f:
	for line in reads:
		try:
			sigar = line.split('\t')[5]
			s={}
			for l,n in zip(re.split(r'[0-9]+', sigar)[1:], re.split(r'[A-Z]', sigar)[:-1]):
				if l in s:
					s[l] += int(n)
				else:
					s[l] = int(n)
			if s['M'] > 40:
				out_f.write(line+'\n')
		except:
			continue

run('samtools view -bS %s.sam -o %s.bam' % (pref, pref), shell=True) 
run('rm %s.sam' % pref, shell=True)

# Determine target
if panel=='unknown':
	rbam=pref+'.bam'
	target=run('%s/panel_search.py %s "AODHRD15,AODCPV1,AODABCV1"'% (path['scripts'], rbam),  stdout=PIPE, shell=True).stdout.decode('utf-8').strip()
	open('%s/panel_name.csv' % out_f, 'w').write(target)
	panel = target.split('\t')[0]
bed = path.get(panel)

# Subset reads by regions
inp = pref+'.bam'
out = pref + ".t1.bam"
run("%s view -b -L %s -o %s %s" % ('samtools', bed, out, inp), shell = True)

# Make group of reads
inp = out
out = pref + ".t2.bam"
run("java -jar %s AddOrReplaceReadGroups I=%s O=%s RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=tumor" %
		(path['picard'], inp, out), shell = True)

# Sort and index
inp = out
out = pref + ".target.group.bam"
run("%s sort %s -o %s" % ('samtools', inp, out), shell = True)
run("%s index %s" % ('samtools', out), shell = True)

# Remove reads with -F2, -f 2048, from different chr, from different loci (more then 200000bp)
inp = out
out = pref + "_filtered.bam"
run("%s view -H %s > %s" % ('samtools', inp, pref + ".target.group.sam"), shell=True)
run("%s view -f 2 -F 2048 %s | awk -F \'\t\' \'$7==\"=\"\' | awk -F \'\t\' \'sqrt(($4 - $8)^2) < 200000\' >> %s" % ('samtools', inp, pref + ".target.group.sam"), shell=True)
run("%s view -b %s | samtools sort - -o %s" % ('samtools', pref + ".target.group.sam", out), shell=True)
run("%s index %s" % ('samtools', out), shell=True)

# Clean
for p in [".t1.bam", ".t2.bam", ".target.group.bam", ".target.group.bam.bai", ".target.group.sam"]:
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
os.rename(pref + ".bqsr.bam", pref + ".bam")
os.rename(pref + ".bqsr.bai", pref + ".bam.bai")

# Clean
for p in ["_filtered.bam", "_filtered.bam.bai", ".realigned.bam", ".realigned.bai", ".recall.table"]:
	x = pref + p
	try:
		os.remove(x)
	except:
		print("Error while deleting ", x)

run('rm %s%s_adaptertrimmed_R1.fq %s%s_adaptertrimmed_R1_unpaired.fq %s%s_adaptertrimmed_R2.fq %s%s_adaptertrimmed_R2_unpaired.fq' % (out_folder, name, out_folder, name, out_folder, name, out_folder, name), shell=True)