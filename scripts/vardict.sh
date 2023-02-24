#!/bin/bash
echo "USAGE: sh vardict.sh <PATH_TO_BAM> <PATH_TO_OUTPUT_DIR> <PATH_TO_REF_SEQ> <PATH_TO_TARGETS>"

# Constants
resources=resources.csv
vardict="$(echo "$(grep -m1 "vardict" "$resources")" | cut -d ',' -f2)"
echo "Running VarDict: $vardict"
#bcftools="$(echo "$(grep -m1 "bcftools" "$resources")" | cut -d ',' -f2)"
echo "Running bcftools: $bcftools"

# Variables
echo "BAM: $1";
echo "OUT: $2";
echo "REF: $3";
panel="$(echo "$(grep -m1 "$4" "$resources")" | cut -d ',' -f2)"
echo "TARGET: $4; BED: $panel";
AF_THR="0.01"
echo "Min AF: $AF_THR";

# Extracting sample name from the path
dir_name="$(dirname -- "$1")"
sample_name="$(basename -- "$(dirname -- "$(dirname -- "$1")")")"
echo "Sample name: $sample_name";

# Running VarDict
echo
echo "Command:"
echo "vardict -G $3 -f $AF_THR -N $sample_name -b $1 -c 1 -S 2 -E 3 -g 4 $panel | teststrandbias.R | var2vcf_valid.pl -E -f $AF_THR > $2/variant.vardict.vcf"
$vardict/vardict -G $3 -f $AF_THR -N $sample_name -b $1 -c 1 -S 2 -E 3 -g 4 $panel | teststrandbias.R | var2vcf_valid.pl -E -f $AF_THR | \
	bcftools norm -m -any -f $3 -O z -o $2/variant.vardict.no_filters.vcf.gz

tabix -f $2/variant.vardict.no_filters.vcf.gz

# Testing whether output VCF exists
if [ -f "$2/variant.vardict.no_filters.vcf.gz" ]; then
    echo "VCF: $2/variant.vardict.no_filters.vcf.gz"
fi

# Filtering
bcftools filter $2/variant.vardict.no_filters.vcf.gz -e '((INFO/AF * INFO/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 2.0) || (INFO/DP < 10) || (INFO/QUAL < 45)))' \
	-O z -o	$2/variant.vardict.vcf.gz

tabix -f $2/variant.vardict.vcf.gz