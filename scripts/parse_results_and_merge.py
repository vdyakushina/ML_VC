import glob
import os
import sys
import numpy as np
import pandas as pd
from cyvcf2 import VCF
from datetime import datetime


def input_path():
    from func import read_resources
    res = "/pipeline/scripts/resources.csv"
    path = read_resources(res)
    return path


def process_tsv(tsv_fn, var_hgvs, is_snv, common_tsv):
    is_snv_table = False if 'indel' in tsv_fn.split('/')[-1].lower() else True
    if not is_snv_table:
        curr_columns = columns_tsv_indel
        curr_columns_disp = columns_tsv_disp_indel
    else:
        curr_columns = columns_tsv
        curr_columns_disp = columns_tsv_disp
    print('Found', tsv_fn.split('/')[-1], ':',
          '\ttable supposedly contains', 'SNV' if is_snv_table is True else 'INDEL',
          '\tSample variant is', 'SNV' if is_snv is True else 'INDEL')
    if is_snv_table == is_snv:
        try:
            new_tsv = pd.read_table(tsv_fn)
            print('These are the variants you are looking for')
        except pd.errors.EmptyDataError:
            new_tsv = pd.DataFrame(columns=curr_columns)
            print('This TSV is empty:', tsv_fn.split('/')[-1])

        # print('Variants in this sample (selected columns)', new_tsv.shape, ':\n', new_tsv[columns_tsv_disp[:-1]])

        new_tsv = new_tsv.loc[new_tsv['variant_name'] == var_hgvs]
        new_tsv_shape = new_tsv.shape

        # TODO rework this in accordance with VCF function, so that the missing variants are actually added
        if new_tsv.empty:
            # new_tsv.shape != new_tsv_shape:
            print('\nThe variant of this sample', '(' + sample_variant + ')', 'is not found in TSV', tsv_fn.split('/')[-1])
            # adding empty row - but only for suitable variants (SNV/INDEL)
            new_tsv = pd.DataFrame([[var_hgvs] + [''] * len(curr_columns[1:])], columns=curr_columns)
        else:
            new_tsv['samples'] = sample_variant
            print('\nFound the variant of this sample', '(' + sample_variant + ')', ':\n', new_tsv[curr_columns_disp])

        old_shape = common_tsv.shape
        common_tsv = pd.concat([common_tsv, new_tsv], ignore_index=True)
        new_shape = common_tsv.shape
        if common_tsv.empty:
            print('\nCommon TSV is still empty', new_shape, '- not updated with', tsv_fn.split('/')[-1], new_tsv.shape)
        else:
            if old_shape == new_shape:
                print('\nCommon TSV', old_shape, 'is not updated with', tsv_fn.split('/')[-1], new_tsv.shape)
            else:
                if new_tsv_shape[0] == 0:
                    print('\nCommon TSV', old_shape, 'is updated with NOT_FOUND row from', tsv_fn.split('/')[-1], new_tsv.shape, '\n', common_tsv[curr_columns_disp])
                else:
                    print('\nCommon TSV', old_shape, 'is updated with', tsv_fn.split('/')[-1], new_tsv.shape, '\n', common_tsv[curr_columns_disp])

        return common_tsv
    else:
        print('These are not the variants you are looking for')
        return common_tsv


def process_vcf(vcf_fn, sample_variant, var_chr, var_pos, var_RA, is_snv, common_vcf):
    is_snv_table = False if 'indel' in vcf_fn.split('/')[-1].lower() else True
    print('Found', vcf_fn.split('/')[-1], ':',
          '\ttable supposedly contains', 'SNV' if is_snv_table is True else 'INDEL',
          '\tSample variant is', 'SNV' if is_snv is True else 'INDEL')
    if is_snv_table == is_snv:
        try:
            new_vcf = pd.read_table(vcf_fn, comment='#', header=None)
            new_vcf.columns = columns_vcf
            print('These are the variants you are looking for')
        except pd.errors.EmptyDataError:
            new_vcf = pd.DataFrame(columns=columns_vcf)
            print('This VCF is empty:', vcf_fn.split('/')[-1])

        new_vcf = new_vcf.astype(str)
        # if not new_vcf.empty:
        #     print('Variants in this sample (selected columns)', new_vcf.shape, ':\n', new_vcf[['#CHROM', 'POS', 'REF', 'ALT']])

        # new_vcf_copy = new_vcf.copy()
        new_vcf = new_vcf.loc[(new_vcf['#CHROM'] == var_chr[0]) &
                              (new_vcf['POS'] == var_pos[0]) &
                              (new_vcf['REF'] == var_RA[0]) &
                              (new_vcf['ALT'] == var_RA[1])]
        new_vcf_shape = new_vcf.shape

        # TODO expand the table with new sample names (discuss if this is needed)
        if new_vcf.empty:
            # new_vcf.shape != new_vcf_shape:
            print('\nThe variant of this sample', '(' + sample_variant + ')', 'is not found in vcf', vcf_fn.split('/')[-1])
            new_vcf = pd.DataFrame([[var_chr[0], var_pos[0], '.', var_RA[0], var_RA[1],
                                    '.', '.', 'sample_variant='+sample_variant+'.NOT_FOUND', '.', sample_variant]], columns=columns_vcf)
        else:
            print('\nFound the variant of this sample', '(' + sample_variant + ')', ':\n',
                  new_vcf[['#CHROM', 'POS', 'REF', 'ALT']])
            new_vcf['INFO'] = new_vcf['INFO'] + ';sample_variant='+sample_variant+'.FOUND'
            new_vcf['samples'] = sample_variant

        old_shape = common_vcf.shape
        common_vcf = pd.concat([common_vcf, new_vcf], ignore_index=True)
        if common_vcf.empty:
            print('common_vcf', old_shape, '+', vcf_fn.split('/')[-1], new_vcf.shape, ': common table is empty')
        else:
            new_shape = common_vcf.shape
            if old_shape == new_shape:
                print('\nCommon VCF', old_shape, 'is not updated with', vcf_fn.split('/')[-1], new_vcf.shape)
            else:
                if new_vcf_shape[0] == 0:
                    print('\nCommon VCF', old_shape, 'is updated with NOT_FOUND row from', vcf_fn.split('/')[-1], new_vcf.shape, '\n', common_vcf)
                else:
                    print('\nCommon VCF', old_shape, 'is updated with', vcf_fn.split('/')[-1], new_vcf.shape, '\n', common_vcf)

        return common_vcf
    else:
        print('These are not the variants you are looking for')
        return common_vcf


def write_vcf_as_table(common_table, out_vcf):
    curr_date = datetime.now().strftime("%A, %d. %B %Y %I:%M%p")
    comment = """##fileformat=VCFv4.2\n##fileDate=""" + curr_date
    comment += """\n##source=parse_results_and_merge.py
##INFO=<ID=VD_Q,Number=1,Type=Float,Description="VarDict quality">
##INFO=<ID=VD_FILTER,Number=1,Type=Integer,Description="VarDict FILTER (PASS/excuse)">
##INFO=<ID=SINVICT,Number=1,Type=Integer,Description="Sinvict discovered the variant or not">
##INFO=<ID=SC_FILTER,Number=1,Type=String,Description="Scalpel FILTER (PASS/excuse)">
##INFO=<ID=SC_PENALTY,Number=1,Type=Float,Description="Scalpel penalty for having filters MS, LowCov, or HighCov (4 for each)">
##INFO=<ID=Mutect2_TO_FILTER,Number=1,Type=String,Description="Results of filtration">
##INFO=<ID=Mutect2_TO_STRQ,Number=1,Type=Float,Description="Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors">
##INFO=<ID=Mutect2_TO_TLOD,Number=1,Type=Float,Description="Log 10 likelihood ratio score of variant existing versus not existing">
##INFO=<ID=Mutect2_TO_RPA,Number=1,Type=Float,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
##INFO=<ID=Mutect2_TO_RPA_1,Number=1,Type=Float,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
##INFO=<ID=Mutect2_TO_RPA_2,Number=1,Type=Float,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
##INFO=<ID=Mutect2_TO_SB,Number=1,Type=Float,Description="Strand Bias calculated as ratio  allele proportion  in strand with minimum proportion to one with maximum proportion">
##INFO=<ID=Strelka_TO_FILTER,Number=1,Type=String,Description="Filter results in tumor normal paires for strelka to">
##INFO=<ID=Strelka_TO_SB,Number=1,Type=Float,Description="Strand Bias calculated as ratio allele proportion  in strand with minimum proportion to one with maximum proportion">
##INFO=<ID=Strelka_TO_SNVHPOL,Number=.,Type=Integer,Description="SNV contextual homopolymer length">
##INFO=<ID=Strelka_TO_REFREP,Number=.,Type=Integer,Description="Number of times RU is repeated in reference">
##INFO=<ID=Strelka_TO_IDREP,Number=.,Type=Integer,Description="Number of times RU is repeated in indel allele">
##INFO=<ID=Strelka_TO_QUAL,Number=.,Type=Integer,Description="Quality">
##INFO=<ID=Strelka_TNP_FILTER,Number=1,Type=String,Description="Filter results in tumor normal pairs for strelka tnp">
##INFO=<ID=Mutect2_TNP_FILTER,Number=1,Type=String,Description="Filter results in tumor normal pairs">
##INFO=<ID=sample_variant,Number=1,Type=String,Description="Indicator that the sought variant was not found in a given sample">\n"""
    contig = ''
    for c in sorted(set(common_table['#CHROM'].to_list())):
        contig += '##contig=<ID=' + c + '>\n'
    comment += contig

    with open(out_vcf, 'a') as f:
        f.write(comment)
        common_table.to_csv(f, sep='\t', header=True, index=False)


if __name__ == '__main__':
    # TODO Из RAW .vcf и .tsv файлов выдергиваем строку (по одной с каждого файла), соответствующую исследуемому варианту
    #  Если такой строки нет, то на обсуждение.
    #  Добавляем эти строки в общие .vcf и .tsv файлы - это и есть тестовая выборка с истинно-положительными вариантами

    # pd.set_option('display.max_columns', 100)
    pd.set_option('display.width', None)

    # path = input_path()

    samples_dir = '/home/aurora/Atlas/Illumina_pipeline_development/samples_tsv'
    # samples_dir = '/home/v/Atlas/Illumina_pipeline_development/samples_tsv'

    # wildcard_tsv_indel = 'variant.raw.INDEL.tsv'
    # wildcard_tsv_snv = 'variant.raw.SNV.tsv'
    # wildcard_vcf_indel = 'variant.raw.INDEL.vcf'
    # wildcard_vcf_snv = 'variant.raw.SNV.vcf'
    # wildcard_gz_indel = wildcard_vcf_indel + '.gz'
    # wildcard_gz_snv = wildcard_vcf_snv + '.gz'

    # glob_handle_tsv_1 = glob.glob(samples_dir + '/*/analysis/' + wildcard_tsv_indel)
    # glob_handle_tsv_2 = glob.glob(samples_dir + '/*/analysis/' + wildcard_tsv_snv)
    # glob_handle_tsv = glob_handle_tsv_1 + glob_handle_tsv_2

    # glob_handle_vcf_1 = glob.glob(samples_dir + '/*/analysis/' + wildcard_vcf_indel)
    # glob_handle_vcf_2 = glob.glob(samples_dir + '/*/analysis/' + wildcard_vcf_snv)
    # glob_handle_vcf_3 = glob.glob(samples_dir + '/*/analysis/' + wildcard_gz_indel)
    # glob_handle_vcf_4 = glob.glob(samples_dir + '/*/analysis/' + wildcard_gz_snv)
    # glob_handle_vcf = glob_handle_vcf_1 + glob_handle_vcf_2 + glob_handle_vcf_3 + glob_handle_vcf_4

    columns_tsv = ['variant_name', 'Mutect2_TO_FILTER', 'Mutect2_TO_STRQ', 'Mutect2_TO_TLOD', 'Mutect2_TNP_FILTER', 'Strelka_TO_FILTER', 'Strelka_TO_QUAL', 'Strelka_TNP_FILTER', 'SC_FILTER', 'SC_PENALTY', 'SINVICT', 'VD_FILTER', 'VD_Q', 'SGA_VAF', 'SGA_SB', 'SGA_RepeatUnit', 'SGA_RepeatRefCount', 'SGA_DP', 'samples']
    columns_tsv_indel = [x for x in columns_tsv if 'SC_' not in x]
    columns_tsv_disp = ['variant_name', 'Mutect2_TO_FILTER', 'Mutect2_TNP_FILTER', 'Strelka_TO_FILTER', 'Strelka_TNP_FILTER', 'SC_FILTER',  'SINVICT', 'VD_FILTER', 'VD_Q', 'SGA_VAF', 'samples']
    columns_tsv_disp_indel = [x for x in columns_tsv_disp if 'SC_' not in x]
    columns_vcf = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'samples']

    line_l = '\n~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~\n'
    line_m = '\n~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~\n'
    line_s = '\n~~~~~ ~~~~~ ~~~~~ ~~~~~ ~~~~~\n'
    is_snv = True

    print("glob.glob(samples_dir + '/*/analysis')", glob.glob(samples_dir + '/*/analysis'))

    mutation_types = ['INDEL', 'SNV']
    for mt in mutation_types:
        print(line_l, 'Processing mutations of type:', mt)

        if mt == 'INDEL':
            common_tsv = pd.DataFrame(columns=columns_tsv_indel)
        else:
            common_tsv = pd.DataFrame(columns=columns_tsv)
        common_vcf = pd.DataFrame(columns=columns_vcf)

        curr_wildcard = 'variant.raw.' + mt
        wildcard_tsv = curr_wildcard + '.tsv'
        wildcard_vcf = curr_wildcard + '.vcf'
        wildcard_gz = wildcard_vcf + '.gz'

        for sample in glob.glob(samples_dir + '/*/analysis'):
            if os.path.isdir(sample):
                sample_variant = sample.split('/')[-2]
                var_chr = sample_variant.split('_')[:1]
                var_pos = [''.join(char for char in sample_variant.split('_')[1] if char.isnumeric())]
                var_RA = ''.join(char for char in sample_variant if not char.isnumeric()).split('_')[-2:]
                is_snv = all([len(x) == 1 for x in var_RA])
                var_hgvs = sample_variant.split('_')[0] + ':' + \
                           sample_variant.split('_')[1] + '>' + \
                           sample_variant.split('_')[2]
                print(line_m + 'sample_variant', sample_variant, '|', var_chr + var_pos + var_RA, '|', var_hgvs)

                # Creating common table with columns from the first non-empty file, adding found variants to the table
                # glob_handle_tsv = glob.glob(sample+'/'+wildcard_tsv_indel) + glob.glob(sample+'/'+wildcard_tsv_snv)
                glob_handle_tsv = glob.glob(sample+'/'+wildcard_tsv)
                if not glob_handle_tsv:
                    print('No', mt, '.tsv found')
                for tsv_fn in glob_handle_tsv:
                    common_tsv = process_tsv(tsv_fn, var_hgvs, is_snv, common_tsv)
                    print(line_s)

                # glob_handle_vcf = glob.glob(sample+'/'+wildcard_vcf_indel) + glob.glob(sample+'/'+wildcard_vcf_snv)\
                #                   + glob.glob(sample+'/'+wildcard_gz_indel) + glob.glob(sample+'/'+wildcard_gz_snv)
                glob_handle_vcf = glob.glob(sample+'/'+wildcard_vcf) + glob.glob(sample+'/'+wildcard_gz)
                if not glob_handle_vcf:
                    print('No', mt, '.vcf found')
                for vcf_fn in glob_handle_vcf:
                    common_vcf = process_vcf(vcf_fn, sample_variant, var_chr, var_pos, var_RA, is_snv, common_vcf)
                    print(line_s)

        if not common_tsv.empty:
            out_tsv = samples_dir + '/common_table.' + mt + '.tsv'
            common_tsv.to_csv(out_tsv, sep='\t', header=True, index=False)

        if not common_vcf.empty:
            out_vcf = samples_dir + '/common_table.' + mt + '.vcf'
            if os.path.isfile(out_vcf):
                os.remove(out_vcf)
            write_vcf_as_table(common_vcf, out_vcf)

    # for i in sorted(glob_handle_tsv):
    #     print(i)
    # process_tsv

    # for i in sorted(glob_handle_vcf):
    #     process_vcf(i)
