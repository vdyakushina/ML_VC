#!/bin/bash

# Initiate
fq1=$1
fq2=$2
res=/pipeline/scripts/resources.csv
cpu=8
wd=`pwd`

# Get resources
declare -A path=() 
while IFS="," read -r key value; do
	path[${key}]+=${value}
done < <(cut -d"," -f1,2 ${res})

# Run fastQC
${path[fastqc]} --quiet --threads ${cpu} --outdir ${wd} ${fq1} ${fq2}  
