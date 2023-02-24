#!/usr/bin/env python3
import sys, os, glob, gzip, sqlite3, re, argparse
import pandas as pd
from func import read_resources
from subprocess import run, call, PIPE
from pathlib import Path
from datetime import datetime
from collections import Counter

sys.path.insert(1, str(Path(__file__).parent.resolve()))
sys.path.insert(1, os.environ['CLAUDIA'])
import AODDB
import libI

##Initiate

parser = argparse.ArgumentParser(description="Performes interpretation. Includes Annovar, SpliceAI, oc and algorithm for interpretation.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-sample", help="path to vcf file. File must be .vcf.gz.")
parser.add_argument("-outputf", help="Name of output file, not path, will bw used as prefix.")
args = vars(parser.parse_args())

sample=os.path.realpath(args["sample"])
output_f=args["outputf"]

res = os.path.dirname(os.path.realpath(__file__)) + '/resources.csv'
path = read_resources(res)

tmp_output_f=output_f.replace('.tsv','')
folder=os.path.dirname(sample)
tmp_folder=folder+'/tmp_'+tmp_output_f+'/'
fname=os.path.basename(sample)

if os.path.isdir(tmp_folder):
	run('rm -rf %s' % tmp_folder, shell=True)

if len(glob.glob('%s/%s.tsv' % (folder, tmp_output_f)))>=1:
	run('rm %s/%s.tsv' % (folder, tmp_output_f), shell=True)

db=pd.read_csv(path["db_ExonicFunc.refGene"], sep='\t', header=None)
variant_summary=path["variant_summary"]
var_sum=pd.read_csv(variant_summary, sep='\t', compression='gzip', dtype=str)
submission_summary=path["submission_summary"]
header=path["header.info"]
transcripts_dictiannary=path["transcripts_dictiannary"]
transcripts_dictiannary=pd.read_csv(transcripts_dictiannary, sep='\t').set_index("ensembl_transcript_id")["refseq_mrna"].to_dict()

def interpretation(annotated_input, splice_input, chasm_input, annotated_output):
	external = ['Color','GeneDx','Invitae','Ambry','Athena','Counsyl','Quest']
	annotated=pd.read_csv(annotated_input, sep='\t')
	splice=pd.read_csv(splice_input, sep='\t')
	splice=pd.concat([splice[["# CHROM", "POS", "REF", "ALT"]], splice["SpliceAI"].str.split("|", expand=True)[[2, 3, 4, 5]]], axis=1).rename({2:"DS_AG", 3:"DS_AL", 4:"DS_DG", 5:"DS_DL"}, axis=1)
	oc_table=pd.read_csv(oc_input, sep='\t')
	oc_table=oc_table.fillna('.')
	annotated=annotated.merge(splice, on=["# CHROM", "POS", "REF", "ALT"])
	annotated=annotated.merge(oc_table, on=["# CHROM", "POS", "REF", "ALT"], how='left')
	for item in ["MutPred_score", "AF_exome", "AF_genome", "CADD_phred"]:
		annotated.loc[annotated[item]== ".",item]=0
		annotated.loc[annotated[item]== "-",item]=0
	for item in ["DS_AG", "DS_AL", "DS_DG", "DS_DL"]:
		annotated.loc[annotated[item].isnull(),item]=0
	annotated.loc[annotated['SIFT_pred']=='Damaging', 'SIFT_pred']='D'
	annotated.loc[annotated['SIFT_pred']=='Tolerated', 'SIFT_pred']='T'
	for index, line in annotated.iterrows():
		print(index)
		gene=line["Gene.refGene"].split("\\")[0]
		chrm=line["# CHROM"]
		chrm_n=chrm.removeprefix("chr")
		pos=line["POS"]
		refn=line["REF"]
		altn=line["ALT"]
		var=var_sum.query(' Chromosome == "%s"  & PositionVCF == "%s" & ReferenceAlleleVCF == "%s" & AlternateAlleleVCF == "%s" ' % (str(chrm_n), str(pos), refn, altn))
		try:
			annotated.loc[index, "ClinVar_ID"]=var['VariationID'].values[0]
			annotated.loc[index, "CLNSIG_clinvar"]=var['ClinicalSignificance'].values[0]
			if bool(re.search(r'conflict', var['ClinicalSignificance'].values[0].lower())):
				var_ID = var["VariationID"].values[0] 
				all_submissions=run('zcat %s | awk -F \'\t\' -v OFS=\'\t\' \'$1==\'%s\'\'' % (submission_summary, var_ID), shell=True, stdout=PIPE).stdout.decode().strip().split('\n')
				sig_counts=Counter([item.split('\t')[1] for item in all_submissions])
				sig={}
				for submitter in ['Color', 'GeneDx', 'Invitae', 'Ambry', 'Athena', 'Quest', 'Counsyl']:
					sig[submitter]={}
					submittions=list(filter(lambda x:submitter in x, all_submissions))
					for i in submittions:
						try:
							sig[submitter][datetime.strptime(i.split('\t')[2], '%b %d, %Y')]=i.split('\t')[1]
						except ValueError:
							continue
					try:
						sig[submitter]=sig[submitter][max(sig[submitter].keys())]
					except ValueError:
						sig[submitter]='No'
						continue
					annotated.loc[index, submitter]=sig[submitter]
				try:
					sig_counts.pop('not provided')
					cnt=sum(sig_counts.values())
				except KeyError:
					cnt=sum(sig_counts.values())
				pgc=sum([sig_counts[i] for i in list(filter(lambda x:re.findall(r'pathogenic', x.lower()), sig_counts.keys()))])/cnt
				bn=sum([sig_counts[i] for i in list(filter(lambda x:re.findall(r'benign', x.lower()), sig_counts.keys()))])/cnt
				annotated.loc[index, "Benign_probability"]=bn
				annotated.loc[index, "Pathogenic_probability"]=pgc
				annotated.loc[index, "Number_evidence"]=cnt
				annotated.loc[index, "Submissions_counts"]=','.join([i+":"+str(j) for i,j in zip(sig_counts.keys(), sig_counts.values())])
		except:
			annotated.loc[index, "CLNSIG_clinvar"]='NaN'
		mutation_name = chrm+':'+ str(pos)+refn+'>'+ altn
		exfunc=line["ExonicFunc.refGene"]
		func=line["Func.refGene"]
		s = len(list(filter(re.compile("splice").match, list(map(lambda x: x.info['variantConsequence'], AODDB.Mutation(mutation_name).PrimaryAnnotation.VariantConsequences)))))
		print(f'splice {s}')
		splice_score=max(line[["DS_AG","DS_AL","DS_DG","DS_DL"]])
		if float(line["CADD_phred"])>=15:
			line["CADD_phred"]="D"
		else:
			line["CADD_phred"]="T"
		if float(line["MutPred_score"])>=0.8:
			line["MutPred_score"]="D"
		else:
			line["MutPred_score"]="T"
		predictors=list(line[["SIFT_pred", "fathmm_MKL_coding_pred", "MetaLR_pred",  "MutationTaster_pred", "CADD_phred", "MutPred_score"]])
		if gene in set(db[7].str.split('(', expand=True)[0].str.split('=', expand=True)[1]):
			if gene == "ATM" or gene == "BRCA2":
				value=max(db[(db[7].str.contains(gene)) & ((db[9]=="frameshift deletion") | (db[9]=="frameshift insertion") | (db[9]=="frameshift substitution") | (db[9]=="stopgain"))][1])
				strand="positive"
			else:
				value=min(db[(db[7].str.contains(gene)) & ((db[9]=="frameshift deletion") | (db[9]=="frameshift insertion") | (db[9]=="frameshift substitution") | (db[9]=="stopgain"))][1])
				strand="negative"
		if (len(db[(db[0]==chrm) & (db[1]==int(pos)) & (db[3]==refn) & (db[4]==altn)])>=1):
			annotated.loc[index, "annotation_somatic"]="Pathogenic"
			annotated.loc[index, "annotation_germline"]="Pathogenic"
			annotated.loc[index, "level"]="n1"
			continue
		else:
			try:
				if (bool(re.search(r'conflict', var['ClinicalSignificance'].values[0].lower()))) & (cnt >=10) & (pgc >= 0.8):
					annotated.loc[index, "annotation_somatic"]="Pathogenic"
					annotated.loc[index, "annotation_germline"]="Pathogenic"
					annotated.loc[index, "level"]="n2"
					continue
				else:
					pass
			except:
				pass 
		if float(line["AF_exome"]) >=0.001 or float(line["AF_genome"]) >=0.001:
			annotated.loc[index, "annotation_somatic"]="Neutral"
			annotated.loc[index, "annotation_germline"]="Neutral"
			annotated.loc[index, "level"]="n3"
			continue 
		elif exfunc== "stoploss":
			annotated.loc[index, "annotation_somatic"]="Pathogenic"
			annotated.loc[index, "annotation_germline"]="Pathogenic"
			annotated.loc[index, "level"]="n4"
			continue
		else: 
			try:
				if (bool(re.search(r'conflict', var['ClinicalSignificance'].values[0].lower()))) and (bn >= 0.8) and (cnt >=10): 
					annotated.loc[index, "annotation_somatic"]="Neutral"
					annotated.loc[index, "annotation_germline"]="Neutral"
					annotated.loc[index, "level"]="n10"
					continue
				else:
					pass
			except:
				pass
		if exfunc== "stopgain" or exfunc=="frameshift_deletion" or exfunc=="frameshift_insertion":
			if (gene in set(db[7].str.split('(', expand=True)[0].str.split('=', expand=True)[1])) and (strand=="positive" and int(pos) <= value): 
				annotated.loc[index, "annotation_somatic"]="Pathogenic"
				annotated.loc[index, "annotation_germline"]="Pathogenic"
				annotated.loc[index, "level"]="n5"
				continue
			elif (gene in set(db[7].str.split('(', expand=True)[0].str.split('=', expand=True)[1])) and (strand=="negative" and int(pos) >= value): 
				annotated.loc[index, "annotation_somatic"]="Pathogenic"
				annotated.loc[index, "annotation_germline"]="Pathogenic"
				annotated.loc[index, "level"]="n5"
				continue
			elif gene in set(db[7].str.split('(', expand=True)[0].str.split('=', expand=True)[1]):
				annotated.loc[index, "annotation_somatic"]="Neutral"
				annotated.loc[index, "annotation_germline"]="Neutral"
				annotated.loc[index, "level"]="n12"
				continue
			else:
				annotated.loc[index, "annotation_somatic"]="Likely pathogenic"
				annotated.loc[index, "annotation_germline"]="Likely pathogenic"
				annotated.loc[index, "level"]="n5"
				continue
		else:
			try:
				if var['ClinicalSignificance'].values[0].split(',')[0].lower() in ["benign", "benign/likely benign"]:
					annotated.loc[index, "annotation_somatic"]="Neutral"
					annotated.loc[index, "annotation_germline"]="Neutral"
					annotated.loc[index, "level"]="n6"
					continue
				else:
					pass
			except:
				pass
		try:
			if (var['ClinicalSignificance'].values[0].lower().startswith('pathogenic')) or (var['ClinicalSignificance'].values[0].split(',')[0].lower() == "likely pathogenic"):
				annotated.loc[index, "annotation_somatic"]="Likely pathogenic"
				annotated.loc[index, "annotation_germline"]="Likely pathogenic"
				annotated.loc[index, "level"]="n7"
				continue
			else:
				pass
		except:
			pass
		try:
			if (s > 0) and (float(splice_score)>=0.8):
				func="splicing"
				annotated.loc[index, "Func.refGene"]="splicing" 
				annotated.loc[index, "annotation_somatic"]="Pathogenic"
				annotated.loc[index, "annotation_germline"]="Pathogenic"
				annotated.loc[index, "level"]="n8"
				continue
		except:
			pass
		if func != "exonic": 
			annotated.loc[index, "annotation_somatic"]="Neutral"
			annotated.loc[index, "annotation_germline"]="Neutral"
			annotated.loc[index, "level"]="n9"
		elif exfunc == "nonframeshift_insertion" or exfunc == "nonframeshift_deletion":
			annotated.loc[index, "annotation_somatic"]="VUS"
			annotated.loc[index, "annotation_germline"]="VUS"
			annotated.loc[index, "level"]="n11"
		elif exfunc!= "nonsynonymous_SNV":
			annotated.loc[index, "annotation_somatic"]="Neutral"
			annotated.loc[index, "annotation_germline"]="Neutral"
			annotated.loc[index, "level"]="n12"
		elif predictors.count("D") in range(5,7):
			annotated.loc[index, "annotation_somatic"]="Suspected deleterious"
			annotated.loc[index, "annotation_germline"]="VUS"
			annotated.loc[index, "level"]="n13"
		elif  predictors.count("D") in range (3,5):
			annotated.loc[index, "annotation_somatic"]="VUS"
			annotated.loc[index, "annotation_germline"]="VUS"
			annotated.loc[index, "level"]="n15"
		elif predictors.count("D") in range (0,3):
			annotated.loc[index, "annotation_somatic"]="VUS"
			annotated.loc[index, "annotation_germline"]="VUS"
			annotated.loc[index, "level"]="n14"
	annotated['name'] = annotated['# CHROM'].astype('str') + ':'  + annotated['POS'].astype('str') + annotated['REF'].astype('str') + '>' + annotated['ALT'].astype('str')
	annotated=annotated.drop(['# CHROM', 'POS', 'REF', 'ALT'], 1)
	annotated=pd.concat([annotated[['name', 'Gene.refGene', 'Func.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'AF_exome', 'AF_genome', 'ClinVar_ID', 'CLNSIG_clinvar', 'Number_evidence', 'Submissions_counts', 'Benign_probability', 'Pathogenic_probability']],annotated.drop(['name', 'Gene.refGene', 'Func.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'AF_exome', 'AF_genome', 'ClinVar_ID', 'CLNSIG_clinvar', 'Number_evidence', 'Submissions_counts', 'Benign_probability', 'Pathogenic_probability'], 1)],1)
	for ex in external:
		if ex not in list(annotated):
			annotated[ex]=None
	for i, row in annotated.iterrows():
		if (libI._unify_interpretation(row['annotation_germline']) == 'vus'):
			external_neutral_found = False
			external_damaging_found = False
			for ex in external:
				if (libI._unify_interpretation(row[ex]) == 'neutral'):
					external_neutral_found = True
				if (libI._unify_interpretation(row[ex]) == 'damaging'):
					external_damaging_found = True
			if ((row['Pathogenic_probability'] == 0)or
				row['Benign_probability'] > row['Pathogenic_probability']):
				if ((external_neutral_found == True)and(external_damaging_found == False)):
					annotated.loc[i, 'annotation_somatic'] = 'neutral'
					annotated.loc[i, 'annotation_germline'] = 'neutral'
					annotated.loc[i, 'level'] = 'n10.1'
					continue
			if ((row['Benign_probability'] == 0)or
				(row['Pathogenic_probability'] > row['Benign_probability'])):
				if ((external_damaging_found == True)and(external_neutral_found == False)):
					annotated.loc[i, 'annotation_somatic'] = 'damaging'
					annotated.loc[i, 'annotation_germline'] = 'damaging'
					annotated.loc[i, 'level'] = 'n10.2'
					continue
	annotated.to_csv(annotated_output, sep='\t', index=False)

##Check and prepare input
if os.path.isfile(sample) and '.vcf' in sample:
	print('file exists')
	if '.gz' in sample and not os.path.isfile(sample+'tbi'): #Check tabix file exists
		print("tbi doesn't exist, performing tabix")
		run("tabix -f %s" % (sample), shell = True)
elif os.path.isfile(sample) and '.tsv' in sample[-4:]:
	print('tsv file exists input must be vcf file')
	exit()
else:
	print('smth wrong with input')
	exit()

if not os.path.isdir(tmp_folder):
	Path(tmp_folder).mkdir(exist_ok = True)

## Run annotation with ANNOVAR
with gzip.open('%s' % sample, 'rb') as inf, open('%s/variant.raw.%s.vcf' % (tmp_folder, tmp_output_f), 'w') as outf:
	for line in inf.read().decode().split('\n'):
		if '#' in line:
			outf.write(line+'\n')
		elif 'chr' in line:
			line='\t'.join(line.split('\t')[:8])+'\tGT\t0/1\n'
			outf.write(line)

run("%stable_annovar.pl %s/variant.raw.%s.vcf %s --thread 6 -buildver hg19 -out %s/%s -remove -protocol dbnsfp42a,refGene,gnomad211_exome,gnomad211_genome -operation f,g,f,f -arg ,,, -nastring . --vcfinput" % \
	(path["annovar"], tmp_folder, tmp_output_f, path["annovar_db"], tmp_folder, tmp_output_f), shell = True)
run("sed -i 's/fathmm-MKL_coding_pred/fathmm_MKL_coding_pred/g' %s/%s.hg19_multianno.vcf" % (tmp_folder,tmp_output_f), shell = True) 

## Run annotation with SPLICEAI
run("spliceai -I %s -O %s/%s_splice.vcf -R %s -A grch37" % (sample, tmp_folder, tmp_output_f, path['ref']), shell = True)

## Run annotation with Open Cravat
run("oc run -n %s_oc -a vest polyphen2 sift cscape_coding provean cadd_exome chasmplus -l hg19 -d %s %s" %(tmp_output_f, tmp_folder, sample), shell = True)
oc_sq = sqlite3.connect("%s/%s_oc.sqlite" % (tmp_folder, tmp_output_f))
oc_a = pd.read_sql_query("SELECT original_input__chrom, extra_vcf_info__pos, extra_vcf_info__ref, extra_vcf_info__alt, cadd_exome__score, cadd_exome__phred, chasmplus__pval, chasmplus__score, cscape_coding__score, cscape_coding__rankscore, provean__score, provean__rankscore, provean__prediction, sift__prediction, sift__confidence, sift__score, sift__rankscore, vest__score, vest__pval from variant", oc_sq)
oc_a=oc_a.rename(columns={'original_input__chrom':'# CHROM', 'extra_vcf_info__pos':'POS', 'extra_vcf_info__ref':'REF', 'extra_vcf_info__alt':'ALT', 'sift__prediction':'SIFT_pred', 'provean__prediction':"PROVEAN_pred",'cadd_exome__phred':'CADD_phred'})
oc_a.to_csv("%s/%s_oc_final.tsv" % (tmp_folder, tmp_output_f), sep='\t', index=False)
oc_sq.close()

## Select required columnes
vcf=pd.read_csv("%s/%s.hg19_multianno.vcf" % (tmp_folder, tmp_output_f), sep='\t', comment='#', header=None)
items=["Gene.refGene", "Func.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene", "AF_exome", "AF_genome", "fathmm_MKL_coding_pred", "MetaLR_pred", "MutationTaster_pred",  "MutPred_score"]
for index,line in vcf.iterrows():
	for item in items: 
		if item in line[7]:
			vcf.loc[index, item]=line[7].split(item+"=")[1].split(';')[0]
		else:
			vcf.loc[index, item]='.'
vcf=pd.concat([vcf[[0,1, 3,4]], vcf[items]], axis=1).rename({0:"# CHROM", 1:"POS", 3:"REF", 4:"ALT"}, axis=1)
vcf.to_csv("%s/%s_final.vcf" % (tmp_folder, tmp_output_f), sep='\t', index=False)

splice=pd.read_csv("%s/%s_splice.vcf" % (tmp_folder, tmp_output_f), sep='\t', comment='#', header=None, engine='python')
for index,line in splice.iterrows():
	if "SpliceAI=" in line[7]:
		splice.loc[index, "SpliceAI"]=line[7].split("SpliceAI=")[1]
	else:
		splice.loc[index, "SpliceAI"]='.'
splice=pd.concat([splice[[0,1,3,4]], splice["SpliceAI"]], axis=1).rename({0:"# CHROM", 1:"POS", 3:"REF", 4:"ALT"}, axis=1)
splice.to_csv("%s/%s_splice_final.vcf" % (tmp_folder, tmp_output_f), sep='\t', index=False)

## Run Interpretation function
annotated_input = "%s/%s_final.vcf" % (tmp_folder, tmp_output_f)
splice_input = "%s/%s_splice_final.vcf" % (tmp_folder, tmp_output_f)
oc_input = "%s/%s_oc_final.tsv" % (tmp_folder, tmp_output_f)
annotated_output = "%s/%s" % (folder, output_f)
interpretation(annotated_input, splice_input, oc_input, annotated_output)

## Remove some tmp files
run('rm -rf %s' % (tmp_folder), shell = True)



