#!/bin/bash

resources=resources.csv
ref="$(echo "$(grep -m1 "ref" "$resources")" | cut -d ',' -f2)"
sga="$(echo "$(grep -m1 "sga" "$resources")" | cut -d ',' -f2)"

vcf_indel_name=variant.raw.INDEL.JOIN.vcf.gz
sga_indel_name=variant.raw.INDEL.JOIN.SGA.vcf
vcf_snv_name=variant.raw.SNV.JOIN.vcf.gz
sga_snv_name=variant.raw.SNV.JOIN.SGA.vcf


for bam_file in /pipeline/data/samples_test/join/analysis/provisional.bam; do

dir_name="$(dirname -- "$bam_file")"
echo $dir_name

echo "$sga somatic-variant-filters --tumor-bam $bam_file --normal-bam $bam_file --reference $ref $dir_name/$vcf_indel_name --annotate-only > $dir_name/$sga_indel_name"
$sga somatic-variant-filters --tumor-bam $bam_file --normal-bam $bam_file --reference $ref $dir_name/$vcf_indel_name --annotate-only > $dir_name/$sga_indel_name
echo "bgzip -f $dir_name/$sga_indel_name"
bgzip -f $dir_name/$sga_indel_name
echo "tabix -f $dir_name/$sga_indel_name.gz"
tabix -f $dir_name/$sga_indel_name.gz

echo "mv $dir_name/$sga_indel_name.gz $dir_name/$vcf_indel_name"
echo "mv $dir_name/$sga_indel_name.gz.tbi $dir_name/$vcf_indel_name.tbi"
mv $dir_name/$sga_indel_name.gz $dir_name/$vcf_indel_name
mv $dir_name/$sga_indel_name.gz.tbi $dir_name/$vcf_indel_name.tbi


echo "$sga somatic-variant-filters --tumor-bam $bam_file --normal-bam $bam_file --reference $ref $dir_name/$vcf_snv_name --annotate-only > $dir_name/$sga_snv_name"
$sga somatic-variant-filters --tumor-bam $bam_file --normal-bam $bam_file --reference $ref $dir_name/$vcf_snv_name --annotate-only > $dir_name/$sga_snv_name
echo "bgzip -f $dir_name/$sga_snv_name"
bgzip -f $dir_name/$sga_snv_name
echo "tabix -f $dir_name/$sga_snv_name.gz"
tabix -f $dir_name/$sga_snv_name.gz

echo "mv $dir_name/$sga_snv_name.gz $dir_name/$vcf_snv_name"
echo "mv $dir_name/$sga_snv_name.gz.tbi $dir_name/$vcf_snv_name.tbi"
mv $dir_name/$sga_snv_name.gz $dir_name/$vcf_snv_name
mv $dir_name/$sga_snv_name.gz.tbi $dir_name/$vcf_snv_name.tbi


if [ -f "$dir_name/$vcf_indel_name" ]; then
    echo "INDEL annotations written to $dir_name/$vcf_indel_name"
fi

if [ -f "$dir_name/$vcf_snv_name" ]; then
    echo "SNV annotations written to $dir_name/$vcf_snv_name"
fi


echo
done

