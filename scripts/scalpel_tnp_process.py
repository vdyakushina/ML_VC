#!/usr/bin/env python3
import sys
import pandas as pd
import os
from subprocess import run
from func import read_resources

folder=sys.argv[1]
name='tumor'

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
#path_bcftools = path['bcftools']
path_bcftools = 'bcftools'

print('folder', folder)

#Prepare annotation for info field
try:
  vcf1=pd.read_csv("%s/1_variant.scalpel_tnp.vcf.gz" % folder, sep='\t', comment='#', header=None, compression='gzip')[range(0,7)].rename({0:"#CHROM", 1:"POS", 2:"ID", 3:"REF", 4:"ALT", 5:"QUAL", 6:"FILTER_1"}, axis=1)
except:
  vcf1=pd.DataFrame({"#CHROM":[None], "POS":0, "ID":[None], "REF":[None], "ALT":[None], "QUAL":[None], "FILTER_1":''})
  run("mv %s/1_variant.scalpel_tnp.vcf.gz %s/1_variant.scalpel_tnp_empty.vcf.gz" % (folder, folder), shell=True, check=True)
  run("mv %s/1_variant.scalpel_tnp.vcf.gz.tbi %s/1_variant.scalpel_tnp_empty.vcf.gz.tbi" % (folder, folder), shell=True, check=True)
try:
  vcf2=pd.read_csv("%s/2_variant.scalpel_tnp.vcf.gz" % folder, sep='\t', comment='#', header=None, compression='gzip')[range(0,7)].rename({0:"#CHROM", 1:"POS", 2:"ID", 3:"REF", 4:"ALT", 5:"QUAL", 6:"FILTER_2"}, axis=1)
except:
  vcf2=pd.DataFrame({"#CHROM":[None], "POS":0, "ID":[None], "REF":[None], "ALT":[None], "QUAL":[None], "FILTER_2":''})
  run("mv %s/2_variant.scalpel_tnp.vcf.gz %s/2_variant.scalpel_tnp_empty.vcf.gz" % (folder, folder), shell=True, check=True)
  run("mv %s/2_variant.scalpel_tnp.vcf.gz.tbi %s/2_variant.scalpel_tnp_empty.vcf.gz.tbi" % (folder, folder), shell=True, check=True)

merged=vcf1.merge(vcf2.drop(["ID","QUAL"], axis=1), on=["#CHROM","POS","REF","ALT"], how="outer")

for i in range(3, 11):
  try:
    vcf=pd.read_csv("%s/%s_variant.scalpel_tnp.vcf.gz" % (folder, i), sep='\t', comment='#', header=None, compression='gzip')[range(0,7)].rename({0:"#CHROM", 1:"POS", 2:"ID", 3:"REF", 4:"ALT", 5:"QUAL", 6:"FILTER_%s" %i}, axis=1)
  except:
    vcf=pd.DataFrame({"#CHROM":[None], "POS":0, "ID":[None], "REF":[None], "ALT":[None], "QUAL":[None], "FILTER_%s" %i:''})
    run("mv %s/%s_variant.scalpel_tnp.vcf.gz %s/%s_variant.scalpel_tnp_empty.vcf.gz" % (folder, i, folder, i), shell=True, check=True)
    run("mv %s/%s_variant.scalpel_tnp.vcf.gz.tbi %s/%s_variant.scalpel_tnp_empty.vcf.gz.tbi" % (folder, i, folder, i), shell=True, check=True)  
  merged=merged.merge(vcf.drop(["ID","QUAL"], axis=1), on=["#CHROM","POS","REF","ALT"], how="outer")

merged=merged[~merged["#CHROM"].isnull()]

score={'PASS':0, 'HighChi2score':0, 'LowVaf':0, 'LowAltCnt':0, 'HighAltCnt':0, 'LowCovNormal':0, 'HighCovNormal':0, 'LowCovTumor':2, 'HighCovTumor':0, 'LowVafTumor':0, 'LowAltCntTumor':0, 'LowFisherScore':0, 'HighAltCntNormal':2, 'MS':4, 'HighVafNormal':4}

for index, line in merged.iterrows():
  sum=0
  for col in list(merged)[6:16]:
    if type(line[col])!=float:
      for filter in line[col].split(";"):
        sum+=score[filter]
    else:
       sum+=16   
    merged.loc[index, "RES"]=sum

merged=merged.sort_values(["#CHROM", "POS"])

if len(merged)==0:
	merged['RES']=None

merged[["#CHROM","POS","REF","ALT", "RES"]].to_csv("%s/annotation.txt" % folder, sep='\t', index=False)

run("bgzip -f %s" % folder+"annotation.txt", shell=True, check=True)
run("tabix -f -s1 -b2 -e2 %s" % folder+"annotation.txt.gz", shell=True, check=True)
run("echo \'##INFO=<ID=Scalpel_TNP_FILTER,Number=1,Type=String,Description=\"Filter results in tumor normal paires\">\' > %s" % folder+"header.info.txt", shell=True)

#Add INDELs to variant.raw.INDEL.vcf.gz
run("echo \'##fileformat=VCFv4.3\' > %s" % folder+"process.vcf", shell=True)
run("echo \'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\' >> %s" % (name, folder+"process.vcf"), shell=True)

if len(merged)!=0: 
	run("%s merge -m none --force-sample %s/*_variant.scalpel_tnp.vcf.gz -Ou | %s view -V snps,ref,bnd,other - -Ou | %s annotate -x QUAL,FILTER,INFO,FORMAT,FORMAT/GT - | %s view -s %s -H - >> %s/process.vcf" % (path_bcftools, folder, path_bcftools, path_bcftools, path_bcftools, name, folder), shell=True, check=True)

run("bgzip -f %s" % folder+"process.vcf", shell=True, check=True)
run("tabix -f %s" % folder+"process.vcf.gz", shell=True, check=True)

if not os.path.isfile(folder+"variant.raw.INDEL.vcf.gz"):
  run("%s annotate -a %s/annotation.txt.gz -c CHROM,POS,REF,ALT,INFO/Scalpel_TNP_FILTER:=\"RES\" -h %s/header.info.txt %s/process.vcf.gz -O z -o %s/variant.raw.INDEL.vcf.gz" % (path_bcftools, folder, folder, folder, folder), shell=True, check=True)
  run("tabix -f %s/variant.raw.INDEL.vcf.gz" % folder, shell=True, check=True)
  run("rm %s/process.vcf.gz %s/process.vcf.gz.tbi" % (folder, folder), shell=True, check=True)
else:
  run("%s concat -D -a %s/variant.raw.INDEL.vcf.gz %s/process.vcf.gz | %s annotate -a %s/annotation.txt.gz -c CHROM,POS,REF,ALT,INFO/Scalpel_TNP_FILTER:=\"RES\" -h %s/header.info.txt - -O z -o %s/variant.concat.INDEL.vcf.gz" % (path_bcftools, folder, folder, path_bcftools, folder, folder, folder), shell=True, check=True)
  run("mv %s/variant.concat.INDEL.vcf.gz %s/variant.raw.INDEL.vcf.gz" % (folder, folder), shell=True, check=True)
  run("tabix -f %s/variant.raw.INDEL.vcf.gz" % folder, shell=True, check=True)
  run("rm %s/process.vcf.gz %s/process.vcf.gz.tbi" % (folder, folder), shell=True, check=True)

run("rm %s/annotation.txt.gz %s/annotation.txt.gz.tbi %s/header.info.txt" % (folder, folder, folder), shell=True, check=True)

for i in range(1, 11):
  if os.path.exists('%s/%s_variant.scalpel_tnp_empty.vcf.gz' % (folder, i)):
    run("mv %s/%s_variant.scalpel_tnp_empty.vcf.gz %s/%s_variant.scalpel_tnp.vcf.gz" % (folder, i, folder, i), shell=True, check=True)
  if os.path.exists('%s/%s_variant.scalpel_tnp_empty.vcf.gz.tbi' % (folder, i)):
    run("mv %s/%s_variant.scalpel_tnp_empty.vcf.gz.tbi %s/%s_variant.scalpel_tnp.vcf.gz.tbi" % (folder, i, folder, i), shell=True, check=True)


