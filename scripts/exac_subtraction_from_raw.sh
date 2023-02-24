#!/bin/bash
echo "USAGE: sh exac_subtraction_from_raw.sh <INPUT_VCF_FILE>"
echo

resources=resources.csv
exac="$(echo "$(grep -m1 "exac_abc_filter" "$resources")" | cut -d ',' -f2)"
#bcftools="$(echo "$(grep -m1 "bcftools" "$resources")" | cut -d ',' -f2)"

out_name=$1
echo $out_name
#out_name="${out_name//.vcf.gz}"
#out_name="${out_name//.vcf}"
out_name=$out_name.notcommon.vcf.gz
echo $out_name

echo "$bcftools isec -C $1 $exac -O z -o $out_name"
bcftools isec -C "$1" "$exac" -w 1 -O z -o "$out_name"

tabix -f "$1".notcommon.vcf.gz

#echo "mv $1.notcommon.vcf.gz $1"
#mv $1.notcommon.vcf.gz $1
#echo "tabix -f $1"
#tabix -f $1
