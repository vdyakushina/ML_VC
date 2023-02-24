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
#bcftools=path['bcftools']
#samtools=path['samtools']
bcftools='bcftools'
samtools='samtools'
bed=path[target]+'.gz'
meta=path['meta']+f'{target}_testkit_ILLMN/meta_{target}.csv'

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
  run("%s/generate_validation_sample.py %s %s" % (path['scripts'], meta, norma), shell=True)
  run("%s index %s" % (samtools, norma), shell=True)

  ##Replase sample name to "norma"
  run("java -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGID=%s RGSM=\"norma\" RGLB=LIB RGPL=ILLUMINA RGPU=unit1" %
    (path['picard'], norma, norma+"_RG", pseudo), shell=True)
  run("mv %s %s" % (norma+"_RG", norma), shell=True)
  run("%s index %s" % (samtools, norma), shell=True)
 
  output = folder + str(pseudo) + "_variant.strelka_tnp.vcf.gz"
  
  ##Run Strelka 
  
  run("%s/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam %s \
    --tumorBam %s \
    --referenceFasta %s \
    --callRegions %s \
    --targeted \
    --runDir %s" %
      (path['strelka'], norma, bam, ref, bed, folder+"strelka/"), shell = True)
  
  run("%s" % folder+"/strelka/runWorkflow.py -m local -j 20", shell = True)

  ## Concatinate SNV & INDEL, remove NORMAL from results, normalise, rename
  run("echo \"tumor\" > %s" % (folder+"header.txt"), shell = True)
  run("%s concat -D -a %s | %s norm -m -any -f %s - | %s view -s TUMOR - | %s reheader -s %s | bgzip > %s" %(bcftools, folder+"/strelka/results/variants/somatic.*.vcf.gz", bcftools, ref, bcftools, bcftools, folder+"header.txt", output), shell = True)
  run("tabix %s" % output, shell = True)
  run("rm -rf %s/strelka/ %s/norma.bam %s/norma.bam.bai %s" % (folder, folder, folder, folder+"header.txt"), shell = True)