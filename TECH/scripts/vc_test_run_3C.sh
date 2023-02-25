for bam_file in /pipeline/data/samples_test/7*/analysis/provisional.bam; do
#echo $bam_file

dir_name="$(dirname -- "$bam_file")"
echo $dir_name

#sh /pipeline/scripts/vardict.sh $bam_file $dir_name /pipeline/data/reference/hg19.fa /pipeline/data/targets/AODABCV1.designed.bed
python /pipeline/scripts/vc_process.ALL_types.vardict.py $dir_name

#sh /pipeline/scripts/sinvict.sh $bam_file $dir_name /pipeline/data/reference/hg19.fa /pipeline/data/targets/AODABCV1.designed.bed
python /pipeline/scripts/vc_process.ALL_types.sinvict.py $dir_name

#sh /pipeline/scripts/scalpel.sh $bam_file $dir_name/scalpel /pipeline/data/reference/hg19.fa /pipeline/data/targets/AODABCV1.designed.bed
python /pipeline/scripts/vc_process.INDEL.scalpel.py $dir_name

echo
done

