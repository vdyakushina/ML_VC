# Load a bed file with the average coverage of the target, 
# subset nearest target by threshold value, and estimate the MAPD metric 

import pandas as pd
import io
import sys
import os.path
import subprocess
from statistics import median
from subprocess import run, check_output, PIPE

# Initiate
bed = sys.argv[1]
# TODO: uncomment for debug
#bed = "/home/gennady/proj/atlas-co/data/k1.cov.bed"

maxgap = 100000
resources = '/home/onco-admin/ATLAS_software/aod-pipe/Pipe/AODTKAbatch/pipeline/scripts/resources.csv'

def read_resources(csv):
  df = pd.read_csv(csv, header = None, usecols = [0, 1], names = ["id", "path"])
  return df.set_index('id').T.to_dict('records')[0]

# Get paths to resources
path = read_resources(resources)

# Check input
if not os.path.isfile(bed): sys.exit("%s doesn't exist!" % bed)

# Get gaps between targets
st = run("%s spacing -i %s" % (path['bedtools'], bed), stdout=PIPE, shell=True).stdout.decode('utf-8')

# Convert string into dataframe
obj = io.StringIO(st)
df1 = pd.read_csv(obj, header = None, sep = "\t", na_values=".")

# Subset the nearest targets by threshold
df2 = df1[df1[4] < maxgap]

# Subset column with mean coverage of the target
v = df2[3]

# Build consequtive index
v.index = range(0, len(v))

# Construct the vector of consecutive pairwise 
# differences scaled by the mean
d = [(v[i] - v[i+1])/v.mean() for i in range(0, len(v)-1)]
# Compute the median absolute deviation of d[i]
dev = [abs(value - median(d)) for value in d]
# Compute MAPD
mapd = median(dev)
# Output
print(round(mapd, 2), end = " ")
