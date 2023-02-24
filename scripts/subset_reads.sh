#!/bin/bash

# Initiate
bam=$1
target=$2
pref=$3
res='/home/onco-admin/ATLAS_software/aod-pipe/Pipe/AODTKAbatch/pipeline/scripts/resources.csv'
cpu=8

# Get resources
declare -A path=() 
while IFS="," read -r key value; do
	path[${key}]+=${value}
done < <(cut -d"," -f1,2 ${res})

# Subset reads by target
samtools view -L ${path[${target}]} -b ${bam} | \
samtools sort -@ ${cpu} -O bam -o ${pref}.bam
samtools index ${pref}.bam 
