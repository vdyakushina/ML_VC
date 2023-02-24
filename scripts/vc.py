#!/usr/bin/env python3
# A wrapper script to call for variants 


from subprocess import run
from pathlib import Path
from func import read_resources
import sys
import os.path

res = "/pipeline/scripts/resources.csv"
path = read_resources(res)


### Sinvict
sinvict.py analysis/k1.bqsr.bam AODABCV1 analysis/vcf k1
sinvict_to_vcf.py analysis/vcf/calls_level1.sinvict analysis/vcf/k1

