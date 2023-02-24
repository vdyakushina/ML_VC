#!/bin/bash
# Estimate the sensetivity with EphaGen 

# Initiate
bam=$1
pref=$2
res='/home/onco-admin/ATLAS_software/aod-pipe/Pipe/AODTKAbatch/pipeline/scripts/resources.csv'

# Get resources
declare -A path=() 
while IFS="," read -r key value; do
	path[${key}]+=${value}
done < <(cut -d"," -f1,2 ${res})

# run ephagen
parallel perl ${path[ephagen]}/src/ephagen.pl \
	--bam $bam \
	--ref ${path[ref]} \
	--vcf ${path[ephagen]}/reference/{}.hg19.vcf \
	--out $pref.{}.tsv \
	--out_vcf $pref.{}.vcf ::: BRCA 
