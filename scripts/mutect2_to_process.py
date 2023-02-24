#!/usr/bin/env python3
import sys
import pandas as pd
import os
from subprocess import run
from func import read_resources

folder=sys.argv[1]

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
#path_bcftools = path['bcftools']
path_bcftools = 'bcftools'

#Prepare annotation for info field
try:
  vcf=pd.read_csv("%s/variant.mutect2_to.vcf.gz" % folder, sep='\t', comment='#', header=None, compression='gzip')[range(0,10)].rename({0:"#CHROM", 1:"POS", 2:"ID", 3:"REF", 4:"ALT", 5:"QUAL", 6:"FILTER", 7:"INFO", 8:"FMT", 9:"SMPL"}, axis=1)
except:
  vcf=pd.DataFrame({"#CHROM":[None], "POS":0, "ID":[None], "REF":[None], "ALT":[None], "QUAL":[None], "FILTER":'', "INFO":[None], "FMT":[None], "SMPL":[None]})

vcf=vcf[~vcf["#CHROM"].isnull()]

infos=['STRQ', 'TLOD']
score={'PASS':0, 'germline':0, 'normal_artifact':0, 'orientation':2, 'haplotype':2, 'position':2, 'low_allele_frac':2, 'n_ratio':2, 'strict_strand':4, 'weak_evidence':4, 'slippage':4, 'base_qual':4, 'duplicate':4, 'fragment':4, 'strand_bias':8, 'multiallelic':8, 'contamination':8, 'map_qual':8, 'FAIL':8, 'possible_numt':8, 'clustered_events':64}

for index, line in vcf.iterrows():
  sum=0
  for filter in line['FILTER'].split(";"):
    sum+=score[filter]
  vcf.loc[index, "FILTER"]=sum

for item in infos:
  for index, line in vcf.iterrows():
    if item in line['INFO']:
      value = line['INFO'].split(item+"=")[1].split(';')[0]
    else:
      value = int(0)   
    vcf.loc[index, item] = value

vcf[["#CHROM","POS","REF","ALT", "FILTER", "STRQ", "TLOD"]].to_csv("%s/annotation.txt" % folder, sep='\t', index=False)

run("bgzip -f %s" % folder+"annotation.txt", shell=True)
run("tabix -f -s1 -b2 -e2 %s" % folder+"annotation.txt.gz", shell=True)

run("echo \'##INFO=<ID=Mutect2_TO_FILTER,Number=1,Type=String,Description=\"Results of filtration\">\' > %s" % folder+"header.info.txt", shell=True)
run("echo \'##INFO=<ID=Mutect2_TO_STRQ,Number=1,Type=Float,Description=\"Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors\">\' >> %s" % folder+"header.info.txt", shell=True)
run("echo \'##INFO=<ID=Mutect2_TO_TLOD,Number=1,Type=Float,Description=\"Log 10 likelihood ratio score of variant existing versus not existing\">\' >> %s" % folder+"header.info.txt", shell=True)

name='tumor'

#Add SNV to variant.raw.SNV.vcf.gz
run("echo \'##fileformat=VCFv4.3\' > %s" % folder+"process.vcf", shell=True)
run("echo \'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\' >> %s" % (name, folder+"process.vcf"), shell=True)
run("%s view -V indels,mnps,ref,bnd,other %s/variant.mutect2_to.vcf.gz -Ou | %s annotate -x FILTER,INFO,FORMAT,FORMAT/GT - | %s view -H >> %s/process.vcf" % (path_bcftools, folder, path_bcftools, path_bcftools, folder), shell=True)
run("bgzip -f %s" % folder+"process.vcf", shell=True)
run("tabix -f %s" % folder+"process.vcf.gz", shell=True)
if not os.path.isfile(folder+"variant.raw.SNV.vcf.gz"):
  run("%s annotate -a %s/annotation.txt.gz -c CHROM,POS,REF,ALT,Mutect2_TO_FILTER:=\"FILTER\",Mutect2_TO_STRQ:=\"STRQ\",Mutect2_TO_TLOD:=\"TLOD\" -h %s/header.info.txt %s/process.vcf.gz -O z -o %s/variant.raw.SNV.vcf.gz" % (path_bcftools, folder, folder, folder, folder), shell=True)
  run("tabix -f %s/variant.raw.SNV.vcf.gz" % folder, shell=True)
  run("rm %s/process.vcf.gz %s/process.vcf.gz.tbi" % (folder, folder), shell=True)
else:
  run("%s concat -D -a %s/variant.raw.SNV.vcf.gz %s/process.vcf.gz | %s annotate -a %s/annotation.txt.gz -c CHROM,POS,REF,ALT,Mutect2_TO_FILTER:=\"FILTER\",Mutect2_TO_STRQ:=\"STRQ\",Mutect2_TO_TLOD:=\"TLOD\" -h %s/header.info.txt - -O z -o %s/variant.concat.SNV.vcf.gz" % (path_bcftools, folder, folder, path_bcftools, folder, folder, folder), shell=True)
  run("mv %s/variant.concat.SNV.vcf.gz %s/variant.raw.SNV.vcf.gz" % (folder, folder), shell=True)
  run("tabix -f %s/variant.raw.SNV.vcf.gz" % folder, shell=True)
  run("rm %s/process.vcf.gz %s/process.vcf.gz.tbi" % (folder, folder), shell=True)

#Add INDELs to variant.raw.INDEL.vcf.gz
run("echo \'##fileformat=VCFv4.3\' > %s" % folder+"process.vcf", shell=True)
run("echo \'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\' >> %s" % (name, folder+"process.vcf"), shell=True)
run("%s view -V snps,ref,bnd,other %s/variant.mutect2_to.vcf.gz -Ou | %s annotate -x FILTER,INFO,FORMAT,FORMAT/GT - | %s view -H >> %s/process.vcf" % (path_bcftools, folder, path_bcftools, path_bcftools, folder), shell=True)
run("bgzip -f %s" % folder+"process.vcf", shell=True)
run("tabix -f %s" % folder+"process.vcf.gz", shell=True)
if not os.path.isfile(folder+"variant.raw.INDEL.vcf.gz"):
  run("%s annotate -a %s/annotation.txt.gz -c CHROM,POS,REF,ALT,Mutect2_TO_FILTER:=\"FILTER\",Mutect2_TO_STRQ:=\"STRQ\",Mutect2_TO_TLOD:=\"TLOD\" -h %s/header.info.txt %s/process.vcf.gz -O z -o %s/variant.raw.INDEL.vcf.gz" % (path_bcftools, folder, folder, folder, folder), shell=True)
  run("tabix -f %s/variant.raw.INDEL.vcf.gz" % folder, shell=True)
  run("rm %s/process.vcf.gz %s/process.vcf.gz.tbi" % (folder, folder), shell=True)
else:
  run("%s concat -D -a %s/variant.raw.INDEL.vcf.gz %s/process.vcf.gz | %s annotate -a %s/annotation.txt.gz -c CHROM,POS,REF,ALT,Mutect2_TO_FILTER:=\"FILTER\",Mutect2_TO_STRQ:=\"STRQ\",Mutect2_TO_TLOD:=\"TLOD\" -h %s/header.info.txt - -O z -o %s/variant.concat.INDEL.vcf.gz" % (path_bcftools, folder, folder, path_bcftools, folder, folder, folder), shell=True)
  run("mv %s/variant.concat.INDEL.vcf.gz %s/variant.raw.INDEL.vcf.gz" % (folder, folder), shell=True)
  run("tabix -f %s/variant.raw.INDEL.vcf.gz" % folder, shell=True)
  run("rm %s/process.vcf.gz %s/process.vcf.gz.tbi" % (folder, folder), shell=True)

run("rm %s/annotation.txt.gz %s/annotation.txt.gz.tbi %s/header.info.txt" % (folder, folder, folder), shell=True)