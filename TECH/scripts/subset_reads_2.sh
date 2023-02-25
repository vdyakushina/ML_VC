samtools view -H /pipeline/data/samples_tsv/chr11_108098525G_A/raw_d0cNZaUxL5d7L/d0cNZaUxL5d7L.bam > /pipeline/data/samples_tsv/chr11_108098525G_A/raw_d0cNZaUxL5d7L_s/d0cNZaUxL5d7L_s.sam

while read read_name; do

echo $read_name
samtools view /pipeline/data/samples_tsv/chr11_108098525G_A/raw_d0cNZaUxL5d7L/d0cNZaUxL5d7L.bam | grep ${read_name} >> /pipeline/data/samples_tsv/chr11_108098525G_A/raw_d0cNZaUxL5d7L_s/d0cNZaUxL5d7L_s.sam

done < /pipeline/data/samples_tsv/chr11_108098525G_A/fqbam/only_read_names.txt
