#!/usr/bin/env python3
import sys
import pandas as pd
import os
from subprocess import run
from func import read_resources

(bam, folder) = sys.argv[1:3]

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
#path_bcftools = path['bcftools']
path_bcftools = 'bcftools'
#samtools = path['samtools']
samtools = 'samtools'
#sga = path['sga']
sga='/home/onco-admin/bin/sga/src/SGA/sga'
ref = path['ref']


#Split bam to forvard and revers
print('#Split bam to forvard and revers')
run("%s view -f80 -h %s > %s/provisional.revers.sam" % (samtools, bam, folder), shell=True)
run("%s view -f128 -F16 %s >> %s/provisional.revers.sam" % (samtools, bam, folder), shell=True)
run("%s view -b %s/provisional.revers.sam | samtools sort - -o %s/provisional.revers.bam" % (samtools, folder, folder), shell=True)
run("%s index %s/provisional.revers.bam" % (samtools, folder), shell=True)

run("%s view -f64 -F16 -h %s > %s/provisional.forvard.sam" % (samtools, bam, folder), shell=True)
run("%s view -f128 -f16 %s >> %s/provisional.forvard.sam" % (samtools, bam, folder), shell=True)
run("%s view -b %s/provisional.forvard.sam | samtools sort - -o %s/provisional.forvard.bam" % (samtools, folder, folder), shell=True)
run("%s index %s/provisional.forvard.bam" % (samtools, folder), shell=True)

run("rm %s/provisional.revers.sam %s/provisional.forvard.sam" % (folder, folder), shell=True)

#Count reads with sga
print('#Count reads with sga')
run("%s somatic-variant-filters --tumor-bam %s/provisional.forvard.bam --normal-bam %s/provisional.forvard.bam --reference %s %s/variant.raw.SNV.vcf.gz -t 10 --annotate-only > %s/variant.forvard.SNV.sga" % (sga, folder, folder, ref, folder, folder), shell=True)
run("%s somatic-variant-filters --tumor-bam %s/provisional.forvard.bam --normal-bam %s/provisional.forvard.bam --reference %s %s/variant.raw.INDEL.vcf.gz -t 10 --annotate-only > %s/variant.forvard.INDEL.sga" % (sga, folder, folder, ref, folder, folder), shell=True)
run("%s somatic-variant-filters --tumor-bam %s/provisional.revers.bam --normal-bam %s/provisional.revers.bam --reference %s %s/variant.raw.SNV.vcf.gz -t 10 --annotate-only > %s/variant.revers.SNV.sga" % (sga, folder, folder, ref, folder, folder), shell=True)
run("%s somatic-variant-filters --tumor-bam %s/provisional.revers.bam --normal-bam %s/provisional.revers.bam --reference %s %s/variant.raw.INDEL.vcf.gz -t 10 --annotate-only > %s/variant.revers.INDEL.sga" % (sga, folder, folder, ref, folder, folder), shell=True)

run("%s somatic-variant-filters --tumor-bam %s/provisional.bam --normal-bam %s/provisional.bam --reference %s %s/variant.raw.INDEL.vcf.gz -t 10 --annotate-only > %s/variant.INDEL.sga" % (sga, folder, folder, ref, folder, folder), shell=True)
run("%s somatic-variant-filters --tumor-bam %s/provisional.bam --normal-bam %s/provisional.bam --reference %s %s/variant.raw.SNV.vcf.gz -t 10 --annotate-only > %s/variant.SNV.sga" % (sga, folder, folder, ref, folder, folder), shell=True)

#Prepare annotation for info field
print('#Prepare annotation for info field')

try:
  forvard_SNV=pd.read_csv("%s/variant.forvard.SNV.sga" % folder, sep='\t', comment='#', header=None)[range(0,10)].rename({0:"#CHROM", 1:"POS", 2:"ID", 3:"REF", 4:"ALT", 5:"QUAL", 6:"FILTER", 7:"INFO", 8:"FMT", 9:"SMPL"}, axis=1)
except:
  forvard_SNV=pd.DataFrame({"#CHROM":[None], "POS":0, "ID":[None], "REF":[None], "ALT":[None], "QUAL":[None], "FILTER":'', "INFO":[None], "FMT":[None], "SMPL":[None]})

try:
  revers_SNV=pd.read_csv("%s/variant.revers.SNV.sga" % folder, sep='\t', comment='#', header=None)[range(0,10)].rename({0:"#CHROM", 1:"POS", 2:"ID", 3:"REF", 4:"ALT", 5:"QUAL", 6:"FILTER", 7:"INFO", 8:"FMT", 9:"SMPL"}, axis=1)
except:
  revers_SNV=pd.DataFrame({"#CHROM":[None], "POS":0, "ID":[None], "REF":[None], "ALT":[None], "QUAL":[None], "FILTER":'', "INFO":[None], "FMT":[None], "SMPL":[None]})

try:
  forvard_INDEL=pd.read_csv("%s/variant.forvard.INDEL.sga" % folder, sep='\t', comment='#', header=None)[range(0,10)].rename({0:"#CHROM", 1:"POS", 2:"ID", 3:"REF", 4:"ALT", 5:"QUAL", 6:"FILTER", 7:"INFO", 8:"FMT", 9:"SMPL"}, axis=1)
except:
  forvard_INDEL=pd.DataFrame({"#CHROM":[None], "POS":0, "ID":[None], "REF":[None], "ALT":[None], "QUAL":[None], "FILTER":'', "INFO":[None], "FMT":[None], "SMPL":[None]})

try:
  revers_INDEL=pd.read_csv("%s/variant.revers.INDEL.sga" % folder, sep='\t', comment='#', header=None)[range(0,10)].rename({0:"#CHROM", 1:"POS", 2:"ID", 3:"REF", 4:"ALT", 5:"QUAL", 6:"FILTER", 7:"INFO", 8:"FMT", 9:"SMPL"}, axis=1)
except:
  revers_INDEL=pd.DataFrame({"#CHROM":[None], "POS":0, "ID":[None], "REF":[None], "ALT":[None], "QUAL":[None], "FILTER":'', "INFO":[None], "FMT":[None], "SMPL":[None]})

try:
  SNV=pd.read_csv("%s/variant.SNV.sga" % folder, sep='\t', comment='#', header=None)[range(0,10)].rename({0:"#CHROM", 1:"POS", 2:"ID", 3:"REF", 4:"ALT", 5:"QUAL", 6:"FILTER", 7:"INFO", 8:"FMT", 9:"SMPL"}, axis=1)
except:
  SNV=pd.DataFrame({"#CHROM":[None], "POS":0, "ID":[None], "REF":[None], "ALT":[None], "QUAL":[None], "FILTER":'', "INFO":[None], "FMT":[None], "SMPL":[None]})

try:
  INDEL=pd.read_csv("%s/variant.INDEL.sga" % folder, sep='\t', comment='#', header=None)[range(0,10)].rename({0:"#CHROM", 1:"POS", 2:"ID", 3:"REF", 4:"ALT", 5:"QUAL", 6:"FILTER", 7:"INFO", 8:"FMT", 9:"SMPL"}, axis=1)
except:
  INDEL=pd.DataFrame({"#CHROM":[None], "POS":0, "ID":[None], "REF":[None], "ALT":[None], "QUAL":[None], "FILTER":'', "INFO":[None], "FMT":[None], "SMPL":[None]})

files=[forvard_SNV, revers_SNV, forvard_INDEL, revers_INDEL, SNV, INDEL]
for file in files:
  if len(file[~file['#CHROM'].isnull()]):
    for index, line in file.iterrows():
      VAF = float(line['INFO'].split('TumorVAF=')[1].split(';')[0])
      RepeatUnit = line['INFO'].split('RepeatUnit=')[1].split(';')[0]
      RepeatRefCount = int(float(line['INFO'].split('RepeatRefCount=')[1].split(';')[0].replace(',', '.')))
      DP = int(float(line['INFO'].split('TumorTotalDepth=')[1].split(';')[0].replace(',', '.')))
      file.loc[index, 'VAF']=VAF
      file.loc[index, 'RepeatUnit']=RepeatUnit
      file.loc[index, 'RepeatRefCount']=RepeatRefCount
      file.loc[index, 'DP']=DP
  else:
    for index, line in file.iterrows():
      file.loc[index, 'VAF']=0
      file.loc[index, 'RepeatUnit']=None
      file.loc[index, 'RepeatRefCount']=0
      file.loc[index, 'DP']=0

SNV['RepeatRefCount']=SNV['RepeatRefCount'].astype(int)
INDEL['RepeatRefCount']=INDEL['RepeatRefCount'].astype(int)
SNV['DP']=SNV['DP'].astype(int)
INDEL['DP']=INDEL['DP'].astype(int)

annot_SNV=pd.merge(forvard_SNV[["#CHROM","POS","REF","ALT", "VAF"]], revers_SNV[["#CHROM","POS","REF","ALT", "VAF"]], on=["#CHROM","POS","REF","ALT"], suffixes=('_forvard', '_revers'))
annot_INDEL=pd.merge(forvard_INDEL[["#CHROM","POS","REF","ALT", "VAF"]], revers_INDEL[["#CHROM","POS","REF","ALT", "VAF"]], on=["#CHROM","POS","REF","ALT"], suffixes=('_forvard', '_revers'))

for index, line in annot_SNV.iterrows():
 annot_SNV.loc[index,"SB"]=(min(line['VAF_forvard'], line['VAF_revers'])+0.00001)/(max(line['VAF_forvard'], line['VAF_revers'])+0.00001)

for index, line in annot_INDEL.iterrows():
 annot_INDEL.loc[index,"SB"]=(min(line['VAF_forvard'], line['VAF_revers'])+0.00001)/(max(line['VAF_forvard'], line['VAF_revers'])+0.00001)

pd.merge(annot_SNV[["#CHROM","POS","REF","ALT", "SB"]], SNV[["#CHROM","POS","REF","ALT", "VAF", "RepeatRefCount", "DP", "RepeatUnit"]], on=["#CHROM","POS","REF","ALT"]).to_csv("%s/annotation_SNV.txt" % folder, sep='\t', index=False, float_format='%.5f')
pd.merge(annot_INDEL[["#CHROM","POS","REF","ALT", "SB"]], INDEL[["#CHROM","POS","REF","ALT", "VAF", "RepeatRefCount", "DP", "RepeatUnit"]], on=["#CHROM","POS","REF","ALT"]).to_csv("%s/annotation_INDEL.txt" % folder, sep='\t', index=False, float_format='%.5f')

print('#bgzip SNV')
run("bgzip -f %s" % folder+"annotation_SNV.txt", shell=True)
run("tabix -f -s1 -b2 -e2 %s" % folder+"annotation_SNV.txt.gz", shell=True)

print('#bgzip INDEL')
run("bgzip -f %s" % folder+"annotation_INDEL.txt", shell=True)
run("tabix -f -s1 -b2 -e2 %s" % folder+"annotation_INDEL.txt.gz", shell=True)

run("echo \'##INFO=<ID=SGA_SB,Number=1,Type=Float,Description=\"Strend bias; Ratio of min VAF in revers or forward reads to max VAF in reverse or forward reads; VAF from sga\">\' > %s" % folder+"header.info.txt", shell=True)
run("echo \'##INFO=<ID=SGA_VAF,Number=1,Type=Float,Description=\"Variant allele frequency from sga\">\' >> %s" % folder+"header.info.txt", shell=True)
run("echo \'##INFO=<ID=SGA_RepeatRefCount,Number=1,Type=Integer,Description=\"RepeatRefCount from sga\">\' >> %s" % folder+"header.info.txt", shell=True)
run("echo \'##INFO=<ID=SGA_DP,Number=1,Type=Integer,Description=\"DP from sga\">\' >> %s" % folder+"header.info.txt", shell=True)

#Add VAF & SB to variant.raw.SNV.vcf.gz
print('#Add VAF & SB to variant.raw.SNV.vcf.gz')
run("%s annotate -a %s/annotation_SNV.txt.gz -c CHROM,POS,REF,ALT,SGA_SB:=\"SB\",SGA_VAF:=\"VAF\",SGA_RepeatRefCount:=\"RepeatRefCount\",SGA_DP:=\"SGA_DP\" -h %s/header.info.txt %s/variant.raw.SNV.vcf.gz > %s/variant.raw.SNV.vcf" % (path_bcftools, folder, folder, folder, folder), shell=True)
run("%s annotate -a %s/annotation_INDEL.txt.gz -c CHROM,POS,REF,ALT,SGA_SB:=\"SB\",SGA_VAF:=\"VAF\",SGA_RepeatRefCount:=\"RepeatRefCount\",SGA_DP:=\"SGA_DP\" -h %s/header.info.txt %s/variant.raw.INDEL.vcf.gz > %s/variant.raw.INDEL.vcf" % (path_bcftools, folder, folder, folder, folder), shell=True)
run("bgzip -f %s/variant.raw.SNV.vcf" % folder, shell=True)
run("tabix -f %s/variant.raw.SNV.vcf.gz" % folder, shell=True)
run("bgzip -f %s/variant.raw.INDEL.vcf" % folder, shell=True)
run("tabix -f %s/variant.raw.INDEL.vcf.gz" % folder, shell=True)

run("rm %s/provisional.revers.bam %s/provisional.revers.bam.bai %s/provisional.forvard.bam %s/provisional.forvard.bam.bai %s/variant.forvard.SNV.sga %s/variant.forvard.INDEL.sga %s/variant.revers.SNV.sga %s/variant.revers.INDEL.sga %s/variant.SNV.sga %s/variant.INDEL.sga %s/header.info.txt %s/annotation_SNV.txt.gz %s/annotation_SNV.txt.gz.tbi %s/annotation_INDEL.txt.gz %s/annotation_INDEL.txt.gz.tbi" % (folder, folder, folder, folder, folder, folder, folder, folder, folder, folder, folder, folder, folder, folder, folder), shell=True)
