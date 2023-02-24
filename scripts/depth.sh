#!/bin/bash
# Estimate the coverage 

# Initiate
bam=$1
target=$2
pref=$3
res=/pipeline/scripts/resources.csv

# Get resources
declare -A path=() 
while IFS="," read -r key value; do
	path[${key}]+=${value}
done < <(cut -d"," -f1,2 ${res})

# Remove output file if exists
out=${pref}.cov.txt
[ -f ${out} ] && rm ${out}

# Loop via all targets and estimate the coverage
while read line1; do 
	read -r -a arr1 <<< ${line1}
	region="${arr1[0]}:${arr1[1]}-${arr1[2]}"
	${path[samtools]} coverage --region ${region} --no-header --min-MQ 1 ${bam} >> ${out}
done < ${path[${target}]}

# Get average coverage and uniformity
res=`${path[scripts]}/depth_stat.py ${out} `

# Estimate percentage of reads mapped on target 
all=`${path[samtools]} view -q 1 -c ${bam}`
on=`${path[samtools]} view -q 1 -c -L ${path[${target}]} ${bam}`

ontarget=`python3 -c "f=$on/$all; print('%.2f' % f)"`

echo "${res} ${ontarget}"
