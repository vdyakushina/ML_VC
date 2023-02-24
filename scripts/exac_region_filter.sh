#!/bin/bash
echo "USAGE: sh exac_region_filter.sh <PANEL_NAME> <OUTPUT_DIR>"

resources=resources.csv
exac="$(echo "$(grep -m1 "exac_chr" "$resources")" | cut -d ',' -f2)"
ref="$(echo "$(grep -m1 "ref" "$resources")" | cut -d ',' -f2)"
bcftools="$(echo "$(grep -m1 "bcftools" "$resources")" | cut -d ',' -f2)"

panel="$(echo "$(grep -m1 "$1" "$resources")" | cut -d ',' -f2)"
echo "Panel name: $1";
echo "Panel path: $panel";
echo "Output dir: $2"
echo

echo "$bcftools view $exac --threads 6 -R $panel -O z -o $2/ExAC.regions.$1.vcf.gz"
$bcftools view $exac --threads 6 -R $panel -O z -o $2/ExAC.regions.$1.vcf.gz
tabix -f $2/ExAC.regions.$1.vcf.gz
echo

echo "$bcftools norm -m -any -f [REF] [IN.vcf] > [OUT.vcf]"
$bcftools norm -m -any -f $ref $2/ExAC.regions.$1.vcf.gz -O z -o $2/ExAC.regions.$1.norm.vcf.gz
mv $2/ExAC.regions.$1.norm.vcf.gz $2/ExAC.regions.$1.vcf.gz
tabix -f $2/ExAC.regions.$1.vcf.gz
echo

echo "$bcftools filter $2/ExAC.regions.$1.vcf.gz -i 'INFO/AF > 0.01' -O z -o $2/ExAC.regions.$1.filter.vcf.gz"
$bcftools filter $2/ExAC.regions.$1.vcf.gz -i 'INFO/AF > 0.01' -O z -o $2/ExAC.regions.$1.filter.vcf.gz
tabix -f $2/ExAC.regions.$1.filter.vcf.gz
echo

if [ -f "$2/ExAC.regions.$1.vcf.gz" ]; then
    echo "Panel VCF: $2/ExAC.regions.$1.vcf.gz"
fi
if [ -f "$2/ExAC.regions.$1.filter.vcf.gz" ]; then
    echo "Filtered (AF > 0.01) Panel VCF: $2/ExAC.regions.$1.filter.vcf.gz"
fi
