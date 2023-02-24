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

#Prepare annotation for info field
vcf=pd.read_csv("%s/variant.sinvict.vcf.gz" % folder, sep='\t', comment='#', header=None, compression='gzip')[range(0,8)].rename({0:"#CHROM", 1:"POS", 2:"ID", 3:"REF", 4:"ALT", 5:"QUAL", 6:"FILTER", 7:"INFO"}, axis=1)
vcf.loc[vcf['QUAL']=='.', 'QUAL']=0

vcf["SINVICT"]=1

vcf[["#CHROM","POS","REF","ALT", "SINVICT"]].to_csv("%s/annotation.txt" % folder, sep='\t', index=False)

run("bgzip -f %s" % folder+"annotation.txt", shell=True)
run("tabix -f -s1 -b2 -e2 %s" % folder+"annotation.txt.gz", shell=True)

run("echo \'##INFO=<ID=SINVICT,Number=1,Type=Integer,Description=\"Sinvict found the variant\">\' > %s" % folder+"header.info.txt", shell=True)

#Add SNV to variant.raw.SNV.vcf.gz

run("echo \'##fileformat=VCFv4.3\' > %s" % folder+"process.vcf", shell=True)
run("echo \'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\' >> %s" % (name, folder+"process.vcf"), shell=True)
run("%s view -V indels,mnps,ref,bnd,other %s/variant.sinvict.vcf.gz -Ou | %s annotate -x FILTER,INFO,FORMAT,FORMAT/GT - | %s view -H >> %s/process.vcf" % (path_bcftools, folder, path_bcftools, path_bcftools, folder), shell=True)
run("bgzip -f %s" % folder+"process.vcf", shell=True)
run("tabix -f %s" % folder+"process.vcf.gz", shell=True)
if not os.path.isfile(folder+"variant.raw.SNV.vcf.gz"):
  run("%s annotate -a %s/annotation.txt.gz -c CHROM,POS,REF,ALT,SINVICT:=\"SINVICT\" -h %s/header.info.txt %s/process.vcf.gz -O z -o %s/variant.raw.SNV.vcf.gz" % (path_bcftools, folder, folder, folder, folder), shell=True)
  run("tabix -f %s/variant.raw.SNV.vcf.gz" % folder, shell=True)
  run("rm %s/process.vcf.gz %s/process.vcf.gz.tbi" % (folder, folder), shell=True)
else:
  run("%s concat -D -a %s/variant.raw.SNV.vcf.gz %s/process.vcf.gz | %s annotate -a %s/annotation.txt.gz -c CHROM,POS,REF,ALT,SINVICT:=\"SINVICT\" -h %s/header.info.txt - -O z -o %s/variant.concat.SNV.vcf.gz" % (path_bcftools, folder, folder, path_bcftools, folder, folder, folder), shell=True)
  run("mv %s/variant.concat.SNV.vcf.gz %s/variant.raw.SNV.vcf.gz" % (folder, folder), shell=True)
  run("tabix -f %s/variant.raw.SNV.vcf.gz" % folder, shell=True)
  run("rm %s/process.vcf.gz %s/process.vcf.gz.tbi" % (folder, folder), shell=True)

#Add INDELs to variant.raw.INDEL.vcf.gz
run("echo \'##fileformat=VCFv4.3\' > %s" % folder+"process.vcf", shell=True)
run("echo \'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\' >> %s" % (name, folder+"process.vcf"), shell=True)
run("%s view -V snps,ref,bnd,other %s/variant.sinvict.vcf.gz -Ou | %s annotate -x FILTER,INFO,FORMAT,FORMAT/GT - | %s view -H >> %s/process.vcf" % (path_bcftools, folder, path_bcftools, path_bcftools, folder), shell=True)
run("bgzip -f %s" % folder+"process.vcf", shell=True)
run("tabix -f %s" % folder+"process.vcf.gz", shell=True)
if not os.path.isfile(folder+"variant.raw.INDEL.vcf.gz"):
  print("variant.raw.INDEL.vcf doesn't exist")
  run("%s annotate -a %s/annotation.txt.gz -c CHROM,POS,REF,ALT,SINVICT:=\"SINVICT\" -h %s/header.info.txt %s/process.vcf.gz -O z -o %s/variant.raw.INDEL.vcf.gz" % (path_bcftools, folder, folder, folder, folder), shell=True)
  run("tabix -f %s/variant.raw.INDEL.vcf.gz" % folder, shell=True)
  run("rm %s/process.vcf.gz %s/process.vcf.gz.tbi" % (folder, folder), shell=True)
else:
  print("variant.raw.INDEL.vcf exists")
  run("%s concat -D -a %s/variant.raw.INDEL.vcf.gz %s/process.vcf.gz | %s annotate -a %s/annotation.txt.gz -c CHROM,POS,REF,ALT,SINVICT:=\"SINVICT\" -h %s/header.info.txt - -O z -o %s/variant.concat.INDEL.vcf.gz" % (path_bcftools, folder, folder, path_bcftools, folder, folder, folder), shell=True)
  run("mv %s/variant.concat.INDEL.vcf.gz %s/variant.raw.INDEL.vcf.gz" % (folder, folder), shell=True)
  run("tabix -f %s/variant.raw.INDEL.vcf.gz" % folder, shell=True)
  run("rm %s/process.vcf.gz %s/process.vcf.gz.tbi" % (folder, folder), shell=True)

run("rm %s/annotation.txt.gz %s/annotation.txt.gz.tbi %s/header.info.txt" % (folder, folder, folder), shell=True)