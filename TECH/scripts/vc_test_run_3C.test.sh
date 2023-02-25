for bam_file in /pipeline/data/samples_test/bai_d*/analysis/provisional.bam; do
#echo $bam_file

dir_name="$(dirname -- "$bam_file")"
echo $dir_name

#sh /pipeline/scripts/vardict.sh $bam_file $dir_name /pipeline/data/reference/hg19.fa /pipeline/data/targets/AODABCV1.designed.bed
#python /pipeline/scripts/vardict_process.py $dir_name

sh /pipeline/scripts/sinvict.sh $bam_file $dir_name /pipeline/data/reference/hg19.fa /pipeline/data/targets/AODABCV1.designed.bed
python /pipeline/scripts/sinvict_process.py $dir_name

#sh /pipeline/scripts/scalpel.sh $bam_file $dir_name /pipeline/data/reference/hg19.fa /pipeline/data/targets/AODABCV1.designed.bed
#python /pipeline/scripts/scalpel_process.py $dir_name

#bam_meta=/pipeline/data/samples_test/meta/meta.bam
#sh /pipeline/scripts/scalpel.somatic.sh $bam_file $bam_meta $dir_name/scalpel-somatic /pipeline/data/reference/hg19.fa /pipeline/data/targets/AODABCV1.designed.bed


#python3 vc_process.py $dir_name *INDEL.vcf
#python3 vc_process.INDEL.py $dir_name
#python3 vc_process.SNV.py $dir_name

#python3 vc_process.old.INDEL.py $dir_name
#python3 vc_process.old.SNV.py $dir_name

echo
done

