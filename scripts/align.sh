#!/bin/bash

# Initiate
fq1=$1
fq2=$2
pref=$3
res=/pipeline/scripts/resources.csv
cpu=8

# Get resources
declare -A path=() 
while IFS="," read -r key value; do
	path[${key}]+=${value}
done < <(cut -d"," -f1,2 ${res})

# Align, sort and index
${path[bwa]} mem -t ${cpu} -R "@RG\tID:EXOME\tSM:${pref}" ${path[ref]} ${fq1} ${fq2} \
| ${path[sambamba]} view -S -f bam -t ${cpu} /dev/stdin \
| ${path[sambamba]} sort -t ${cpu} -o ${pref}.bam /dev/stdin

[ -f ${pref}.bam ] || { echo "error: ${pref}.bam doesn't exist"; exit 0; }
${path[sambamba]} index -t ${cpu} ${pref}.bam
