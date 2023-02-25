import glob
import pandas as pd
import numpy as np
import sys
import os
import time
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools
from subprocess import run
import pprint
import random

def input_path():
    from func import read_resources
    res = "/pipeline/scripts/resources.csv"
    path = read_resources(res)
    return path

if __name__ == '__main__':
    print('USAGE:', 'python rename_reads.py <R1.fastq> <R2.fastq>')
    
    R1_new_name = sys.argv[1] + '.renamed';
    
    with open(sys.argv[1]) as R1, \
         open(R1_new_name, 'w') as R1_new:
        for line in R1:
            match = re.match('^@(?P<prefix>.+)(?P<tile>\d+):(?P<x_pos>\d+):(?P<y_pos>\d+)$', line)
            if (match):
                #new_name = match.groupdict[0] + match.groupdict[1] + ':' + match.groupdict[2] + ':' + random.randint(1000, 999999)
                print(match.groupdict())
                #R1_new.write(line)

    sys.exit();

