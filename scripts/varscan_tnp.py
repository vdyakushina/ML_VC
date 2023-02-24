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

/pipeline/tools/samtools-1.11/bin/samtools mpileup -f /pipeline/data/reference/hg19.fa /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/provisional.bam > /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/provisional.pileup
varscan somatic /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/norma.pileup /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/provisional.pileup /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/variant.varscan_tnp --output-vcf

bgzip /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/variant.varscan_tnp.snp.vcf
tabix /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/variant.varscan_tnp.snp.vcf.gz

bgzip /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/variant.varscan_tnp.indel.vcf
tabix /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/variant.varscan_tnp.indel.vcf.gz

/pipeline/tools/bcftools-1.12/bin/bcftools concat -D -a /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/variant.varscan_tnp.snp.vcf.gz /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/variant.varscan_tnp.indel.vcf.gz | /pipeline/tools/bcftools-1.12/bin/bcftools norm -m -any -f /pipeline/data/reference/hg19.fa - -Oz -o /home/gkhvorykh/samples/7edvtxB2cF6I5/analysis/variant.varscan_tnp.vcf.gz