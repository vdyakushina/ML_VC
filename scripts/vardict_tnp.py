#!/usr/bin/env python3
# Run mutect2 in tumor only mode 
import sys
import os
from subprocess import run
from func import read_resources

# Initiate
(bam, folder, ref, target) = sys.argv[1:5]

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
bcftools=path['bcftools']
samtools=path['samtools']
bed=path[target]+'.gz'

# Check input
if not os.path.isfile(bam):
  print(bam, "doesn't exist")
  sys.exit(1)

if not os.path.isfile(bed):
  print(bed, "doesn't exist")
  sys.exit(1)

#run("mkdir %s" % folder, shell=True)

#Run strelka with 10 pseudo norma
for pseudo in range(1,11):

  ## Generate pseudo norma
  norma = folder + "norma.bam"
  run("%s/generate_validation_sample.py %s %s" % (path['scripts'], path['meta'], norma), shell=True)
  run("%s index %s" % (samtools, norma), shell=True)

  ##Replase sample name to "norma"
  run("java -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGID=%s RGSM=\"norma\" RGLB=LIB RGPL=ILLUMINA RGPU=unit1" %
    (path['picard'], norma, norma+"_RG", pseudo), shell=True)
  run("mv %s %s" % (norma+"_RG", norma), shell=True)
  run("%s index %s" % (samtools, norma), shell=True)
 
  output = folder + str(pseudo) + "_variant.strelka_tnp.vcf.gz"
  
  ##Run VarDict 
  
  run("vardict -G /pipeline/data/reference/hg19.fa -N tumor -b "/home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/norma.bam|/home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/provisional.bam" -c 1 -S 2 -E 3 -g 4 /pipeline/data/targets/AODABCV1.designed.bed | testsomatic.R | var2vcf_paired.pl -N "tumor|normal" | bgzip > /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/1_variant.vardict_tnp.vcf.gz" %
      (path['strelka'], norma, bam, ref, bed, folder+"strelka/"), shell = True)
  tabix /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/1_variant.vardict_tnp.vcf.gz
  
  run("rm -rf %s/norma.bam %s/norma.bam.bai" % (folder, folder, folder, folder), shell = True)