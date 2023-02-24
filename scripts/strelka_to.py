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
bed=path[target]+'.gz'

# Check input
if not os.path.isfile(bam):
  print(bam, "doesn't exist")
  sys.exit(1)

if not os.path.isfile(bed):
  print(bed, "doesn't exist")
  sys.exit(1)

#run("mkdir %s" % folder, shell=True)


output = folder + "variant.strelka_to.vcf.gz"
  
##Run Strelka 
  
run("%s/bin/configureStrelkaGermlineWorkflow.py \
  --bam %s \
  --referenceFasta %s \
  --callRegions %s \
  --targeted \
  --runDir %s" %
    (path['strelka'], bam, ref, bed, folder+"strelka/"), shell = True)
  
run("%s" % folder+"/strelka/runWorkflow.py -m local -j 10", shell = True)

## Concatinate SNV & INDEL, remove NORMAL from results, normalise, rename
run("echo \"tumor\" > %s" % (folder+"header.txt"), shell = True)
run("%s concat -D -a %s | %s norm -m -any -f %s - | %s reheader -s %s | bgzip > %s" %('bcftools', folder+"/strelka/results/variants/variants.vcf.gz", 'bcftools', ref, 'bcftools', folder+"header.txt", output), shell = True)
run("tabix %s" % output, shell = True)

run("rm -rf %s %s" % (folder+"strelka/", folder+"header.txt"), shell = True)
