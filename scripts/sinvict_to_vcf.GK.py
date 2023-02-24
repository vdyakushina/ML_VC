#!/usr/bin/env python3

import sys
from re import match
from func import read_resources
from Bio import SeqIO

# Initiate
res = '/pipeline/scripts/resources.csv'
path = read_resources(res)

if len(sys.argv) != 3:
	print("Usage:\t" + sys.argv[0] + "\t<input_sinvict_file_path>\t<output_prefix>")
	exit(0)

def get_nuc(c, p):
	#n = int(c.replace('chr', '')) - 1
	s = records[c].seq
	nuc = s[p - 1]
	return nuc

# Create index for recors in fasta file
records = SeqIO.index(path['ref'], 'fasta')

# Parse the sinvict file
variants = {}
tp = ""
with open(sys.argv[1]) as infile:
	for line in infile:
		line = line.rstrip()
		tokens = line.split()
		chromosome = tokens[0]
		position = int(tokens[1])
		ref = tokens[3]
		alt = tokens[5]
		info = f'sinvict:{tokens[4]},{tokens[6]},{tokens[7]},{tokens[8]},{tokens[9]},{tokens[10]},{tokens[11]}'
		tp = "SNV"

		# Correct records with deletion 
		if match("[-]", alt):
			alt = alt.replace('-', '')
			position = position - 1
			nuc = get_nuc(chromosome, position)
			ref = nuc + alt
			alt = nuc
			tp = "indel"

		# Correct records with insertion
		if match("[+]", alt):
			alt = ref + alt.replace('+', '')
			tp = "indel"

		key = f'{chromosome}:{position}'
		if key not in variants:
			variants[key] = [ref, alt, info, tp]

# Save the results in vcf file
pref = sys.argv[2]
f1 = pref + ".sinvict.snv.vcf"
f2 = pref + ".sinvict.indel.vcf"

outfile1 = open(f1, "w")
outfile2 = open(f2, "w")

outfile1.write("##fileformat=VCFv4.3\n")
outfile1.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
outfile2.write("##fileformat=VCFv4.3\n")
outfile2.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

for key, value in variants.items():
	key_tokens = key.split(':')
	if value[3] == "SNV":
		outfile1.write(key_tokens[0] + "\t" + key_tokens[1] + "\t" + "." + "\t" + value[0] + "\t" + value[1] + "\t" + "." + "\t" + "PASS" + "\t" + value[2] + "\n")

	if value[3] == "indel":
		outfile2.write(key_tokens[0] + "\t" + key_tokens[1] + "\t" + "." + "\t" + value[0] + "\t" + value[1] + "\t" + "." + "\t" + "PASS" + "\t" + value[2] + "\n")
