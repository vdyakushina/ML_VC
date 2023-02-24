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
#bcftools = path['bcftools']
bed=path[target]

# Check input
if not os.path.isfile(bam):
  print(bam, "doesn't exist")
  sys.exit(1)

output = folder + "variant.mutect2_to.vcf.gz"
run("%s --java-options \"-Xmx8G\" Mutect2 \
  -R %s \
  -I %s \
  -L %s \
  --max-reads-per-alignment-start 0 \
  -O %s" %
    (path['gatk'], ref, bam, bed, output), shell = True)

# Put filter tag
input = folder + "variant.mutect2_to.vcf.gz"
output = folder + "variant.mutect2_to.filter.vcf"
run("%s --java-options \"-Xmx8G\" FilterMutectCalls \
 -R %s \
 -V %s \
 -O %s" %
 (path['gatk'], ref, input, output), shell = True)

# Normalize vcf & rename
run("echo \"tumor\" > %s" % (folder+"header.txt"), shell = True)
run("%s norm -m -any -f %s %s | %s reheader -s %s | bgzip > %s" % ('bcftools', ref, output, 'bcftools', folder+"header.txt", input), shell = True)
run("tabix %s " % input, shell = True)
run("rm %s/variant.mutect2_to.filter.vcf %s/variant.mutect2_to.filter.vcf.idx %s/variant.mutect2_to.vcf.gz.stats %s/variant.mutect2_to.filter.vcf.filteringStats.tsv %s" % (folder, folder, folder, folder, folder+"header.txt"), shell = True)

