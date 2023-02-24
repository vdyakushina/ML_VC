#!/usr/bin/env python3
# Run mutect2 in tumor only mode 
import sys
import os
from subprocess import run
from func import read_resources
from pathlib import Path

# Initiate
(bam, folder, ref, target) = sys.argv[1:5]

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
#bcftools=path['bcftools']
#samtools=path['samtools']
bcftools='bcftools'
samtools='samtools'
bed=path[target]
meta=path['meta']+f'{target}_testkit_ILLMN/meta_{target}.csv'

# Check input
if not os.path.isfile(bam):
  print(bam, "doesn't exist")
  sys.exit(1)

if not os.path.isfile(bed):
  print(bed, "doesn't exist")
  sys.exit(1)

#Run scalpel with 10 pseudo norma
for pseudo in range(1,11):

  ## Generate pseudo norma
  norma = folder + "norma.bam"
  run("%s/generate_validation_sample.py %s %s" % (path['scripts'], meta, norma), shell=True)
  run("%s index %s" % (samtools, norma), shell=True)

  ##Replase sample name to "norma"
  run("java -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGID=%s RGSM=\"norma\" RGLB=LIB RGPL=ILLUMINA RGPU=unit1" %
    (path['picard'], norma, norma+"_RG", pseudo), shell=True)
  run("mv %s %s" % (norma+"_RG", norma), shell=True)
  run("%s index %s" % (samtools, norma), shell=True)

  ##Create folders for scalpel
  Path(folder+'/scalpel/main/tumor/logs/').mkdir(exist_ok = True, parents=True)
  Path(folder+'/scalpel/main/normal/logs/').mkdir(exist_ok = True, parents=True)
  Path(folder+'/scalpel/logs/').mkdir(exist_ok = True, parents=True)

  ##Run Scalpel
  run("docker run --rm -v /:/media hanfang/scalpel:0.5.3 scalpel-discovery --somatic --normal /media/%s --tumor /media/%s --bed /media/%s --ref /media/%s --dir /media/%s --numprocs 10 --log" %
    (norma, bam, bed, ref, folder+"scalpel/"), shell = True)
  run("cp %s\"/scalpel/main/somatic.indel.vcf\" %s\"/scalpel/somatic.indel.vcf\"" % (folder, folder), shell = True)
  run("bgzip -f %s\"/scalpel/somatic.indel.vcf\"" % (folder), shell = True)
  run("tabix -f %s\"/scalpel/somatic.indel.vcf.gz\"" % (folder), shell = True)

  ## Normalise vcf and rename
  output = folder + str(pseudo)+"_variant.scalpel_tnp.vcf.gz"
  run("echo \"norma\ntumor\" > %s" % (folder+"header.txt"), shell = True)
  run("%s reheader -s %s %s\"/scalpel/somatic.indel.vcf.gz\" > %s\"/scalpel/variant.indel.vcf.gz\"" %(bcftools, folder+"header.txt", folder, folder), shell = True)
  run("tabix -f %s\"/scalpel/variant.indel.vcf.gz\" " % (folder), shell = True)
  run("%s view -s 'tumor' %s\"/scalpel/variant.indel.vcf.gz\" | %s norm -m -any -f %s - | bgzip -f > %s" %(bcftools, folder, bcftools, ref, output), shell = True)
  run("tabix -f %s " % output, shell = True)
  run("rm -rf %s/scalpel/ %s/norma.bam %s/norma.bam.bai %s" % (folder, folder, folder, folder+"header.txt"), shell = True)
