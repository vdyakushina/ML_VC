#!/usr/bin/env python3
# Run scalpel in tumor only mode 
import sys
import os
from subprocess import run
from func import read_resources
from pathlib import Path

## Initiate
(bam, folder, ref, target) = sys.argv[1:5]

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
#bcftools=path['bcftools']
bcftools='bcftools'
bed=path[target]

## Check input
if not os.path.isfile(bam):
  print(bam, "doesn't exist")
  sys.exit(1)

if not os.path.isfile(bed):
  print(bed, "doesn't exist")
  sys.exit(1)

Path(folder+'/scalpel/main/').mkdir(exist_ok = True, parents=True)
Path(folder+'/scalpel/logs/').mkdir(exist_ok = True, parents=True)

## Run Scalpel

run("docker run --rm -v /:/media hanfang/scalpel:0.5.3 scalpel-discovery --single --bam /media/%s --bed /media/%s --ref /media/%s --dir /media/%s --numprocs 10 --log" %
   (bam, bed, ref, folder+"scalpel/"), shell = True)

run("bgzip -f %s > %s" % (folder+"scalpel/variants.indel.vcf", folder+"scalpel/variants.indel.vcf.gz"), shell = True)
run("tabix -f %s" % folder+"scalpel/variants.indel.vcf.gz", shell = True)

## Normalise
output = folder + "variant.scalpel_to.vcf.gz"
run("echo \"tumor\" > %s" % (folder+"header.txt"), shell = True)
run("%s norm -m -any -f %s %s | %s reheader -s %s | bgzip -f > %s" %(bcftools, ref, folder+"scalpel/variants.indel.vcf.gz", bcftools, folder+"header.txt", output), shell = True)
run("tabix -f %s" % output, shell = True)

run("rm -rf %s %s" % (folder+"scalpel/", folder+"header.txt"), shell = True)

