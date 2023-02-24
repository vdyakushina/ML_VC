#!/bin/bash
# A wrapper script to run VarDictJava 

# Initiate and check input
help=$(cat << EOF

Usage: `basename $0` <bam> target id
	bam: path/to/sample.bam
	target: target id
	id: sample id (default: basename of bam file )

EOF
)

af=0.01 # minimum allele fraction
res=/pipeline/scripts/resources.csv
vdict=/pipeline/tools/VarDictJava

declare -A path=() 
while IFS="," read -r key value; do
	path[${key}]+=${value}
done < <(cut -d"," -f1,2 ${res})

[ -z $1 ] && { echo "Missing arguments... $help"; exit 1; }
[ -z $2 ] && { echo "Missing arguments... $help"; exit 1; }

bam=$1
target=${path[$2]}
id=$3

[ -f ${bam} ] || { echo "${bam} doesn't exist!"; exit 1; }
[ -f ${target} ] || { echo "${target} doesn't exist!"; exit 1; }
[ -z ${id} ] && id=$(basename "${bam}" .bam)

out=${id}.vd

# Show input
cat <<EOF

Input: ${bam}
Target: ${target}
Sample: ${id}
Minimum allele fraction: ${af}
Output prefix: ${out}

EOF

# Call variants
echo "Calling variants..."
$vdict/build/install/VarDict/bin/VarDict -G ${path['ref']} \
	-f $af -N $id -b $bam -z 0 -c 1 -S 2 -E 3 -g 4 ${target} |\
	$vdict/VarDict/teststrandbias.R |\
	$vdict/VarDict/var2vcf_valid.pl -N $id -E -f $af > ${out}.vcf

