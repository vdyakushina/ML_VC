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
  name=bam.split('/')[-3]
  run ("echo %s > %s" % (name, folder+"header.txt"), shell = True)
  run("%s concat -D -a %s | %s norm -m -any -f %s - | %s view -s TUMOR - | %s reheader -s %s | bgzip > %s" %(bcftools, folder+"/strelka/results/variants/somatic.*.vcf.gz", bcftools, ref, bcftools, bcftools, folder+"header.txt", output), shell = True)
  run("tabix %s" % output, shell = True)
  run("rm -rf %s/strelka/ %s/header.txt %s/norma.bam %s/norma.bam.bai" % (folder, folder, folder, folder), shell = True)




/pipeline/tools/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/norma.bam \
    --tumorBam /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/provisional.bam \
    --referenceFasta /pipeline/data/reference/hg19.fa \
    --callRegions /pipeline/data/targets/AODABCV1.designed.bed.gz \
    --targeted \
    --runDir /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/strelka/


/pipeline/tools/bcftools-1.12/bin/bcftools concat -D -a /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/strelka/results/variants/somatic.*.vcf.gz \
| /pipeline/tools/bcftools-1.12/bin/bcftools norm -m -any -f /pipeline/data/reference/hg19.fa - \
| /pipeline/tools/bcftools-1.12/bin/bcftools reheader -s <(echo -e "norma\ntumour") \
| bgzip > /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/smurf/sample-strelka.vcf.gz