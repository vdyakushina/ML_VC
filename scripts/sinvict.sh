#!/bin/bash
echo "USAGE: sh sinvict.sh <PATH_TO_BAM> <PATH_TO_OUTPUT_DIR> <PATH_TO_REF_SEQ> <PATH_TO_TARGETS>"

# Constants
resources=resources.csv
bamreadcount="$(echo "$(grep -m1 "bamreadcount" "$resources")" | cut -d ',' -f2)"
echo "Running Sinvict: $bamreadcount"
sinvict="$(echo "$(grep -m1 "sinvict" "$resources")" | cut -d ',' -f2)"
echo "Running Sinvict: $sinvict"
bcftools="$(echo "$(grep -m1 "bcftools" "$resources")" | cut -d ',' -f2)"
echo "Running bcftools: $bcftools"

# Variables
echo "Running bam-readcount --> SiNVICT"
echo "BAM: $1";
echo "OUT: $2";
echo "REF: $3";
panel="$(echo "$(grep -m1 "$4" "$resources")" | cut -d ',' -f2)"
echo "TARGET: $4; BED: $panel";

## Extracting sample name from the filename
#file_name="$(basename -- "$1")"
#sample_name="${file_name%%.*}"

# Extracting sample name from the path
dir_name="$(dirname -- "$1")"
sample_name="$(basename -- "$(dirname -- "$(dirname -- "$1")")")"
echo "Sample name: $sample_name";

# Running bam-readcount
echo
echo "Running bam-readcount"
mkdir -p $2/bamreadcount
echo "$bamreadcount -f $3 -l $panel $1 > $2/bamreadcount/$sample_name.bamreadcount"
$bamreadcount -f $3 $1 > $2/bamreadcount/$sample_name.bamreadcount

# Running SiNVICT
echo
echo "Running SiNVICT"
for sinvict_file in $2/bamreadcount/*.sinvict; do
if [ -f "$sinvict_file" ]; then
    echo "Found prior *.sinvict file; removing to avoid segfault!"
    rm $sinvict_file
fi
done
echo "$sinvict -t $2/bamreadcount -o $2/bamreadcount"
$sinvict -t $2/bamreadcount -o $2/bamreadcount

#echo "perl /pipeline/scripts/sinvict_callsToVcf.pl $2/bamreadcount/calls_level1.sinvict > $2/variant.sinvict_perlscript.vcf"
#perl /pipeline/scripts/sinvict_callsToVcf.pl $2/bamreadcount/calls_level1.sinvict > $2/variant.sinvict_perlscript.vcf
#echo "DOES NOT RUN sort -k1,1 -k2,2 -t$'\t' $2/variant.sinvict_perlscript.vcf > $2/variant.sinvict_perlscript.sort.vcf"
##sort -k1,1 -k2,2 -t $'\t' $2/variant.sinvict_perlscript.vcf > $2/variant.sinvict_perlscript.sort.vcf
echo "python /pipeline/scripts/sinvict_to_vcf.VL.py $2/bamreadcount/calls_level1.sinvict $2/variant.sinvict.vcf"
python3 ./sinvict_to_vcf.VL.py $2/bamreadcount/calls_level1.sinvict $2/variant.sinvict_vcf.vcf

bgzip -f $2/variant.sinvict_vcf.vcf
tabix -f $2/variant.sinvict_vcf.vcf.gz

# Normalize vcf & rename
bcftools norm -m -any -f ${3} ${2}/variant.sinvict_vcf.vcf.gz | bcftools reheader -s `echo "tumor"` | bcftools sort | bgzip > ${2}/variant.sinvict.vcf.gz
tabix ${2}/variant.sinvict.vcf.gz
rm $2/variant.sinvict_vcf.vcf.gz* 