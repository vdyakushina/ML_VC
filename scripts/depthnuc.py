# Estimates the everage depth from 'samtools depth' output
# and counts some metrics

import sys, os, io
import pandas as pd
from subprocess import run, check_output, PIPE

# Initiate
bam, target, pref = sys.argv[1:4]
resources = '/home/onco-admin/ATLAS_software/aod-pipe/Pipe/AODTKAbatch/pipeline/scripts/resources.csv'

def read_resources(csv):
  df = pd.read_csv(csv, header = None, usecols = [0, 1], names = ["id", "path"])
  return df.set_index('id').T.to_dict('records')[0]

# Load resources
path = read_resources(resources)
bed = path[target]
#samtools=path['samtools']
samtools='samtools'


# Check input
for fn in [bam , bed]:
	if not os.path.isfile(fn): sys.exit("%s doesn't exist!" % fn)

## Estimate the coverage per nucleotide
# Load the results of 'samtools depth'
res = run(
	"""%s depth -b %s %s |
	awk 'BEGIN {OFS=\"\t\"} {print $1, $2-1, $2, $3}' |
	%s intersect -a stdin -b %s -wa -wb |
	awk 'BEGIN {OFS = \"\t\"} {print $5, $6, $7, $4, $8}'""" % (
	samtools, bed, bam, 'bedtools', bed),
	stdout=PIPE, shell=True).stdout.decode('utf-8')

# Convert string into dataframe
obj = io.StringIO(res)
df1 = pd.read_csv(obj, header = None, sep = "\t",
			names = ['chr', 'start', 'end', 'cov', 'target'])

# Estimate the mean coverage
avgcov = round(df1['cov'].mean())#.astype("int64")

# Estimate the coverage for each target region
df2 = df1.groupby(['chr', 'start', 'end'])['cov'].mean().reset_index()

# Round the mean coverage 
df2['cov'] = df2['cov'].round().astype("int64")

# Save the target coverage
output = pref+".cov.bed"
df2.to_csv(output, sep = "\t", header = None, index = None)

# Estimate uniformity
n20 = sum(df2['cov'] >= avgcov * 0.2)
total = df2.shape[0]
unf = round(n20/total, 2)

# Estimate eveness score
escore = check_output(
  "Rscript %s/eveness_score.R %s" % (path['scripts'], output),
  shell = True)

# Estimate mapd
mapd = check_output(
  "python3 %s/mapd.py %s " % (path['scripts'], output),
	shell = True)

# Output the results
print({
	'coverage': avgcov,
	'uniformity': unf,
	'escore': float(escore),
	'mapd': float(mapd)
	})
