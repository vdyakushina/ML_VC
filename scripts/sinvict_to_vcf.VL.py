#!/usr/bin/env python3

import sys
from re import match
from func import read_resources
from Bio import SeqIO
import pandas as pd
import os


def get_nuc(c, p):
	# n = int(c.replace('chr', '')) - 1
	s = records[c].seq
	nuc = s[p - 1]
	return nuc


def correct_record_pos(chromosome, position, ref, alt):
	# Correct records with deletion
	if match("[-]", alt):
		alt = alt.replace('-', '')
		position = position - 1
		nuc = get_nuc(chromosome, position)
		ref = nuc + alt
	return position


def correct_record_ref(chromosome, position, ref, alt):
	# Correct records with deletion
	if match("[-]", alt):
		alt = alt.replace('-', '')
		position = position - 1
		nuc = get_nuc(chromosome, position)
		ref = nuc + alt
	return ref


def correct_record_alt(chromosome, position, ref, alt):
	# Correct records with deletion
	if match("[-]", alt):
		alt = alt.replace('-', '')
		position = position - 1
		nuc = get_nuc(chromosome, position)
		ref = nuc + alt
		alt = nuc

	# Correct records with insertion
	if match("[+]", alt):
		alt = ref + alt.replace('+', '')

	return alt


def correct_record(chromosome, position, ref, alt):
	# Correct records with deletion
	if match("[-]", alt):
		alt = alt.replace('-', '')
		position = position - 1
		nuc = get_nuc(chromosome, position)
		ref = nuc + alt
		alt = nuc

	# Correct records with insertion
	if match("[+]", alt):
		alt = ref + alt.replace('+', '')

	pos_ref_alt = str(position) + ';' + ref + ';' + alt
	return pos_ref_alt  # position, ref, alt


def get_info(s_dp, s_adf, s_adr):
	dp = int(''.join(_ for _ in str(s_dp) if _.isdigit()))
	adf = int(''.join(_ for _ in str(s_adf) if _.isdigit()))
	adr = int(''.join(_ for _ in str(s_adr) if _.isdigit()))
	af = (adf + adr) / dp
	info = f'DP=' + str(dp) + ';' + 'AD=' + str(adf) + ',' + str(adr) + ';' + 'AF=' + str(af)
	return info


if __name__ == '__main__':
	# Initiate
	if len(sys.argv) != 3:
		print("Usage:\t" + sys.argv[0] + "\t<input_sinvict_file_path>\t<output_file_path>")
		exit(0)

	# Create index for records in fasta file
	res = os.getcwd()+"/resources.csv"
	path = read_resources(res)
	records = SeqIO.index(path['ref'], 'fasta')
	# records = SeqIO.index('/home/aurora/Atlas/Files_input/ref/hg19/server/hg19.fa', 'fasta')

	# Parse the sinvict file
	variants_df = pd.read_table(sys.argv[1], header=None)
	variants_df['INFO'] = variants_df.apply(lambda x: get_info(x[4], x[8], x[9]), axis=1)

	variants_df.rename(columns={0: '#CHROM'}, inplace=True)
	variants_df['POS_REF_ALT'] = variants_df.apply(lambda x: correct_record(x['#CHROM'], x[1], x[3], x[5]), axis=1)
	variants_df[['POS', 'REF', 'ALT']] = variants_df['POS_REF_ALT'].str.split(';', expand=True)

	variants_df.drop(columns=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 'POS_REF_ALT'], inplace=True)
	variants_df['ID'] = '.'
	variants_df['QUAL'] = '.'
	variants_df['FILTER'] = '.'
	variants_df['FORMAT'] = '.'
	variants_df['tumor'] = '.'
	variants_df = variants_df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'tumor']]

	variants_df.sort_values(by=['#CHROM', 'POS'], inplace=True)

	# Save the results in vcf file
	comment = '##fileformat=VCFv4.3\n'
	comment += """##source=sinvict_to_vcf.VL.py
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depdth">
##INFO=<ID=AD,Number=1,Type=String,Description="Allelic depths for the ref and alt alleles in the order listed">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele fractions of alternate alleles">\n"""
	contig = ''
	for c in sorted(set(variants_df['#CHROM'].to_list())):
		contig += '##contig=<ID=' + c + '>\n'
	comment += contig

	out_vcf_fn = sys.argv[2]
	if os.path.isfile(out_vcf_fn):
		os.remove(out_vcf_fn)
	with open(out_vcf_fn, 'a') as f:
		f.write(comment)
		variants_df.to_csv(f, sep='\t', header=True, index=False)
		print('Sinvict calls written to VCF:', out_vcf_fn)
