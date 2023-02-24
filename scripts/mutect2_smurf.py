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
bcftools = path['bcftools']
samtools = path['samtools']
bed=path[target]

# Check input
if not os.path.isfile(bam):
  print(bam, "doesn't exist")
  sys.exit(1)


#Run mutect2 with 10 pseudo norma
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
  
  ##Run mutect2
  output = folder + str(pseudo) + "_variant.mutect2_tnp.vcf.gz"
  run("%s --java-options \"-Xmx8G\" Mutect2 \
    -R %s \
    -I %s \
    -I %s \
    -normal norma \
    -L %s \
    --max-reads-per-alignment-start 0 \
    --genotype-germline-sites true \
    -O %s" %
      (path['gatk'], ref, bam, norma, bed, output), shell = True)
  run("rm %s %s" % (norma, norma + ".bai"), shell=True)

  # Put filter tag
  input = folder + str(pseudo) + "_variant.mutect2_tnp.vcf.gz"
  output = folder + str(pseudo) + "_variant.mutect2_tnp.filter.vcf"
  run("%s --java-options \"-Xmx8G\" FilterMutectCalls \
   -R %s \
   -V %s \
   -O %s" %
   (path['gatk'], ref, input, output), shell = True)
  
  ## Normalise vcf and rename
  name=bam.split('/')[-3]
  run ("echo %s > %s" % (name, folder+"header.txt"), shell = True)
  run("%s norm -m -any -f %s %s | %s view -s ^norma - | %s reheader -s %s | bgzip > %s" % (bcftools, ref, output, bcftools, bcftools, folder+"header.txt", input), shell = True)
  run("tabix %s " % input, shell = True)
  run("rm %s" % folder+"header.txt", shell = True)

'''
/pipeline/tools/gatk-4.2.0.0/gatk --java-options "-Xmx8G" Mutect2 \
    -R /pipeline/data/reference/hg19.fa \
    -I /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/provisional.bam \
    -I /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/norma.bam \
    -normal norma \
    -L /pipeline/data/targets/AODABCV1.designed.bed \
    --max-reads-per-alignment-start 0 \
    -genotype-germline-sites true \
    -o /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/smurf/sample-mutect.vcf.gz
'''
 