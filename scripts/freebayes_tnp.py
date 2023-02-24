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

#

freebayes -f /pipeline/data/reference/hg19.fa --pooled-continuous --pooled-discrete -F 0.03 -C 2 /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/provisional.bam /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/norma.bam | /pipeline/tools/bcftools-1.12/bin/bcftools norm -m -any -f /pipeline/data/reference/hg19.fa - -Oz -o /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/variant.freebayes.vcf.gz