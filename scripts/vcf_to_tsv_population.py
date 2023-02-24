#!/usr/bin/env python3

import glob

import numpy as np
from cyvcf2 import VCF, Writer
import pandas as pd
import sys
import os
from subprocess import run


def input_dialog():
    if len(sys.argv) == 1:
        in_vcf_fn = input('Specify the input VCF file: ')
        while not os.path.isfile(in_vcf_fn):
            in_vcf_fn = input('Specify the input VCF file in the correct format: ')

    else:
        in_vcf_fn = sys.argv[1]
        while not os.path.isfile(in_vcf_fn):
            in_vcf_fn = input('Specify the input VCF file in the correct format: ')

    return in_vcf_fn


def process_vcf(vcf_fn):
    in_vcf = VCF(vcf_fn)
    count_all = 0
    variant_table = []
    lens = []
    for v in in_vcf:
        variant_name = v.CHROM + ':' + str(v.start + 1) + v.REF + '>' + v.ALT[0]

        # variant_info_all = ';'.join(['='.join([str(_) for _ in x]) for x in list(v.INFO)])
        variant_info_all = str(v).split('\t')[7].split(';')
        info_Mutect2_TO_FILTER = [x.split('=')[-1] for x in variant_info_all if 'Mutect2_TO_FILTER' in x]
        info_Mutect2_TO_STRQ = [x.split('=')[-1] for x in variant_info_all if 'Mutect2_TO_STRQ' in x]
        info_Mutect2_TO_TLOD = [x.split('=')[-1] for x in variant_info_all if 'Mutect2_TO_TLOD' in x]
        info_Mutect2_TO_RPA = [x.split('=')[-1] for x in variant_info_all if 'Mutect2_TO_RPA' in x]
        info_Mutect2_TO_RPA_1 = [x.split('=')[-1] for x in variant_info_all if 'Mutect2_TO_RPA_1' in x]
        info_Mutect2_TO_RPA_2 = [x.split('=')[-1] for x in variant_info_all if 'Mutect2_TO_RPA_2' in x]
        info_Mutect2_TO_SB = [x.split('=')[-1] for x in variant_info_all if 'Mutect2_TO_SB' in x]
        info_Mutect2_TNP_FILTER = [x.split('=')[-1] for x in variant_info_all if 'Mutect2_TNP_FILTER' in x]
        info_Strelka_TO_FILTER = [x.split('=')[-1] for x in variant_info_all if 'Strelka_TO_FILTER' in x]
        info_Strelka_TO_SB = [x.split('=')[-1] for x in variant_info_all if 'Strelka_TO_SB' in x]
        info_Strelka_TO_SNVHPOL = [x.split('=')[-1] for x in variant_info_all if 'Strelka_TO_SNVHPOL' in x]
        info_Strelka_TO_REFREP = [x.split('=')[-1] for x in variant_info_all if 'Strelka_TO_REFREP' in x]
        info_Strelka_TO_IDREP = [x.split('=')[-1] for x in variant_info_all if 'Strelka_TO_IDREP' in x]
        info_Strelka_TO_QUAL = [x.split('=')[-1] for x in variant_info_all if 'Strelka_TO_QUAL' in x]
        info_Strelka_TNP_FILTER = [x.split('=')[-1] for x in variant_info_all if 'Strelka_TNP_FILTER' in x]
        info_SC_FILTER = [x.split('=')[-1] for x in variant_info_all if 'SC_FILTER' in x]

        info_SC_PENALTY = [x.split('=')[-1] for x in variant_info_all if 'SC_PENALTY' in x]
        info_SC_PENALTY = info_SC_PENALTY[0] if len(info_SC_PENALTY) == 1 else 'N/A'
        info_SC_PENALTY = info_SC_PENALTY if str(info_SC_PENALTY[0]).isdigit() else 'N/A'

        info_SINVICT = [x.split('=')[-1] for x in variant_info_all if 'SINVICT' in x]
        info_VD_FILTER = [x.split('=')[-1] for x in variant_info_all if 'VD_FILTER' in x]
        info_VD_Q = [x.split('=')[-1] for x in variant_info_all if 'VD_Q' in x]

        info_SGA_VAF = [x.split('=')[-1] for x in variant_info_all if 'SGA_VAF' in x]
        info_SGA_SB = [x.split('=')[-1] for x in variant_info_all if 'SGA_SB' in x]
        info_SGA_RepeatUnit = [x.split('=')[-1] for x in variant_info_all if 'SGA_RepeatUnit' in x]
        info_SGA_RepeatRefCount = [x.split('=')[-1] for x in variant_info_all if 'SGA_RepeatRefCount' in x]
        info_SGA_DP = [x.split('=')[-1] for x in variant_info_all if 'SGA_DP' in x]

        # print(variant_name, variant_info)
        variant_table.append(
            {
                'variant_name': variant_name,
                'Mutect2_TO_FILTER': info_Mutect2_TO_FILTER[0] if len(info_Mutect2_TO_FILTER) == 1 else 'N/A',
                'Mutect2_TO_STRQ': info_Mutect2_TO_STRQ[0] if len(info_Mutect2_TO_STRQ) == 1 else 'N/A',
                'Mutect2_TO_TLOD': info_Mutect2_TO_TLOD[0] if len(info_Mutect2_TO_TLOD) == 1 else 'N/A',
                'Mutect2_TNP_FILTER': info_Mutect2_TNP_FILTER[0] if len(info_Mutect2_TNP_FILTER) == 1 else 'N/A',
                'Strelka_TO_FILTER': info_Strelka_TO_FILTER[0] if len(info_Strelka_TO_FILTER) == 1 else 'N/A',
                'Strelka_TO_QUAL': info_Strelka_TO_QUAL[0] if len(info_Strelka_TO_QUAL) == 1 else 'N/A',
                'Strelka_TNP_FILTER': info_Strelka_TNP_FILTER[0] if len(info_Strelka_TNP_FILTER) == 1 else 'N/A',
                'SC_FILTER': info_SC_FILTER[0] if len(info_SC_FILTER) == 1 else 'N/A',
                'SC_PENALTY': info_SC_PENALTY,
                'SINVICT': info_SINVICT[0] if len(info_SINVICT) == 1 else 'N/A',
                'VD_FILTER': info_VD_FILTER[0] if len(info_VD_FILTER) == 1 else 'N/A',
                'VD_Q': info_VD_Q[0] if len(info_VD_Q) == 1 else 'N/A',
                'SGA_VAF': info_SGA_VAF[0] if len(info_SGA_VAF) == 1 else 'N/A',
                'SGA_SB': info_SGA_SB[0] if len(info_SGA_SB) == 1 else 'N/A',
                'SGA_RepeatUnit': info_SGA_RepeatUnit[0] if len(info_SGA_RepeatUnit) == 1 else 'N/A',
                'SGA_RepeatRefCount': info_SGA_RepeatRefCount[0] if len(info_SGA_RepeatRefCount) == 1 else 'N/A',
                'SGA_DP': info_SGA_DP[0] if len(info_SGA_DP) == 1 else 'N/A',
            }
        )

        # lens.append([
        #     len(info_Mutect2_TO_FILTER), len(info_Mutect2_TO_STRQ), len(info_Mutect2_TO_TLOD), len(info_Mutect2_TNP_FILTER),
        #     len(info_Strelka_TO_FILTER), len(info_Strelka_TO_QUAL), len(info_Strelka_TNP_FILTER),
        #     len(info_SC_FILTER), len(info_SC_PENALTY), len(info_SINVICT), len(info_VD_FILTER), len(info_VD_Q)
        # ])

        count_all += 1

    # print('lens', lens)
    # print('lens mean', [np.mean(x) for x in lens])
    print('Parsed', count_all, 'variants in', vcf_fn)
    variant_table = pd.DataFrame(variant_table)
    in_vcf.close()
    return variant_table


if __name__ == '__main__':
    print('USAGE:', 'python vcf_to_tsv.py path/to/vcf')
    # in_vcf_fn = '/home/aurora/Atlas/Illumina_pipeline_development/final_raw_files/variant.raw.INDEL.JOIN.SGA.vcf.gz'
    # in_vcf_fn = '/home/aurora/Atlas/Illumina_pipeline_development/final_raw_files/variant.raw.SNV.JOIN.SGA.vcf.gz'
    in_vcf_fn = input_dialog()

    # Parse VCF files in a given in_vcf_fn
    print('VCF:', in_vcf_fn)
    table = process_vcf(in_vcf_fn)
    print(table)

    from func import read_resources
    res = os.getcwd()+"/resources.csv"
    path = read_resources(res)
    gnomad_genome=path['gnomad_genome_abc']
    gnomad_exome=path['gnomad_exome_abc']

    gnomad={}
    with open('%s' % gnomad_exome, 'r') as in_g:
        next(in_g)
        for line in in_g:
            line=line.strip()
            gnomad[''.join([':'.join(['chr'+str(line.split('\t')[0]), line.split('\t')[1]]), '>'.join([line.split('\t')[3], line.split('\t')[4]])])]=line.split('\t')[5]


    with open('%s' % gnomad_genome, 'r') as in_g:
        next(in_g)
        for line in in_g:
            line=line.strip()
            varg=str(''.join([':'.join(['chr'+str(line.split('\t')[0]), line.split('\t')[1]]), '>'.join([line.split('\t')[3], line.split('\t')[4]])]))
            if varg in gnomad:
                varg_old=gnomad[varg]
                gnomad[varg]=max(float(line.split('\t')[5]), float(varg_old))
            else:
                gnomad[varg]=float(line.split('\t')[5])

    for index, line in table.iterrows():
        if line['variant_name'] not in gnomad:
            table=table.drop(index, 0)
        elif float(gnomad[line['variant_name']])<0.001:
            table=table.drop(index, 0)
        else:
            continue
            
    if 'vcf.gz' in in_vcf_fn[-7:]:
        out_tsv_fn = in_vcf_fn.split('.notcommon')[0][:-7] + '_population' + '.tsv'
    else:
        out_tsv_fn = in_vcf_fn.split('.notcommon')[0][:-4] + '_population' + '.tsv'
    table.to_csv(out_tsv_fn, sep='\t', header=True, index=False)


