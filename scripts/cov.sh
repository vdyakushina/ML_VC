#!/bin/bash
# Estimates the target coverage

# Initiate
bam=$1
target=$2
pref=$3
res=/pipeline/scripts/resources.csv
cpu=8

# Get resources
declare -A path=() 
while IFS="," read -r key value; do
	path[${key}]+=${value}
done < <(cut -d"," -f1,2 ${res})

# Estimate the coverage
${path[mosdepth]} -t ${cpu} -b ${path[${target}]} -n -Q 1 ${pref} ${bam} 

# Count the average coverage
avgcov=`${path[scripts]}/get_average_coverage.py ${pref}.regions.bed.gz`

# Get breadth
breadth=`${path[scripts]}/get_breadth.py ${pref}.mosdepth.region.dist.txt`

# Get uniformity
unf=`${path[scripts]}/get_uniformity.py ${pref}.regions.bed.gz ${avgcov}`

# Plot the comulative coverage
Rscript ${path[scripts]}/plot_comulative_coverage.R ${pref}.mosdepth.region.dist.txt ${avgcov} ${pref}

echo "${avgcov} ${breadth} ${unf}"

exit 0
