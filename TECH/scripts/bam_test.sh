for bam_file in /pipeline/data/samples_tsv/*/raw/*bam; do
echo $bam_file

dir_name="$(dirname -- "$bam_file")"
echo $dir_name

sample_name="$(basename -- "$(dirname -- "$(dirname -- "$bam_file")")")"
echo $sample_name

echo
done
