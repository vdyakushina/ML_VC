#!/usr/bin/env python3
import sys
import pandas as pd
import os
from subprocess import run
from func import read_resources

folder=sys.argv[1]

res = os.getcwd()+"/resources.csv"
path = read_resources(res)
path_bcftools = 'bcftools'

#Prepare annotation for info field
for type in ['SNV', 'INDEL']:
	try:
		vcf=pd.read_csv("%s/variant.UEBA.%s.vcf.gz" % (folder, type), sep='\t', comment='#', header=None, compression='gzip')[range(0,10)].rename({0:"#CHROM", 1:"POS", 2:"ID", 3:"REF", 4:"ALT", 5:"QUAL", 6:"FILTER", 7:"INFO", 8:"FMT", 9:"SMPL"}, axis=1)
	except:
		vcf=pd.DataFrame({"#CHROM":[None], "POS":0, "ID":[None], "REF":[None], "ALT":[None], "QUAL":[None], "FILTER":'', "INFO":[None], "FMT":[None], "SMPL":[None]})

	vcf=vcf[~vcf["#CHROM"].isnull()]

	vcf['variant_name']=vcf['#CHROM'].astype(str)+':'+vcf['POS'].astype(str)+vcf['REF'].astype(str)+'>'+vcf['ALT'].astype(str)
	tsv=pd.read_csv("%s/variant.raw.%s.tsv" % (folder, type), sep='\t')
	tsv=tsv.merge(vcf[['variant_name','QUAL','FILTER','INFO']].rename(columns={'QUAL':'UEBA_QUAL', 'FILTER':'UEBA_FILTER', 'INFO':'UEBA_INFO'}), how='left', on='variant_name')
	tsv.to_csv("%s/variant.raw.%s.tsv" % (folder, type), index=False, sep='\t')
	
	vcf[['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER']].to_csv("%s/annotation.txt" % folder, sep='\t', index=False)
	run("bgzip -f %s" % folder+"/annotation.txt", shell=True)
	run("tabix -f -s1 -b2 -e2 %s" % folder+"/annotation.txt.gz", shell=True)
	
	run("echo \'##INFO=<ID=UEBA_QUAL,Number=1,Type=Float,Description=\"QUAL from UEBA\">\' > %s" % folder+"/header.info.txt", shell=True)
	run("echo \'##INFO=<ID=UEBA_FILTER,Number=1,Type=String,Description=\"FILTER from UEBA\">\' >> %s" % folder+"/header.info.txt", shell=True)
	name='tumor'
	
	run("%s annotate -a %s/annotation.txt.gz -c CHROM,POS,REF,ALT,UEBA_QUAL:=\"QUAL\",UEBA_FILTER:=\"FILTER\" -h %s/header.info.txt %s/variant.raw.%s.vcf.gz -Oz -o %s/variant.concat.%s.vcf.gz" % (path_bcftools, folder, folder, folder, type, folder, type), shell=True)
	run("tabix %s/variant.concat.%s.vcf.gz" %(folder, type), shell=True)
	run("%s annotate -a %s/variant.UEBA.%s.vcf.gz -c INFO %s/variant.concat.%s.vcf.gz -Oz -o %s/variant.raw.%s.vcf.gz" % (path_bcftools, folder, type, folder, type, folder, type), shell=True)
	
	run("tabix -f %s/variant.raw.%s.vcf.gz" % (folder, type), shell=True)
	
	run("rm %s/annotation.txt.gz %s/annotation.txt.gz.tbi %s/header.info.txt %s/variant.concat.%s.vcf.gz" % (folder, folder, folder, folder, type), shell=True)
