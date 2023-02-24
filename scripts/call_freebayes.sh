#!/bin/bash

# Initiate
bam=$1
pref=$2
res=/pipeline/scripts/resources.csv
cpu=8

# Get resources
declare -A path=() 
while IFS="," read -r key value; do
	path[${key}]+=${value}
done < <(cut -d"," -f1,2 ${res})

# Add paths to the tools required
export PATH="${path[freebayes-scripts]}:${PATH}"

# Check input
[ -f ${bam} ] || { echo "${bam} doesn't exists"; exit 1; }

# Call variants with freebayes
bash freebayes-parallel <(fasta_generate_regions.py ${path[ref]}.fai 10000000) ${cpu} -n 2 --haplotype-length 0 -f ${path[ref]} ${bam} > ${pref}.freebayes.raw.vcf 
