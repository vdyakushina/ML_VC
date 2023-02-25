import glob
from cyvcf2 import VCF, Writer
import pandas as pd
from datetime import datetime
import sys
import os
from subprocess import run


# INPUT ARGUMENT: analysis dir, search for other sample_name.*.SNV.vcf and sample_name.*.INDEL.vcf
# PARSE all these files and merge all variants, adding new INFO fields by variant callers
# CREATE tables / VCFs sample_name.raw.SNV.vcf and sample_name.raw.INDEL.vcf with the merged variants
# -- check if raw VCFs exist first

# План - сканировать папку на предмет наличия этого файла и файлов коллеров,
# отдельно вытаскивать все варианты из последних и добавлять в общую таблицу с пометкой коллера.
# Вывод таблицы в vcf.

# TODO
# 	☑ Check existence of the output file in the vc_process.py script
# 	☑ Adjust the script for each caller
# 	- Expand to process indels and SNV in one script without splitting them first


def input_dialog():
    if len(sys.argv) == 1:
        input_folder = input('Specify the input directory: ')
        while not os.path.isdir(input_folder):
            input_folder = input('Specify the input directory in the correct format: ')

    else:
        input_folder = sys.argv[1]
        while not os.path.isdir(input_folder):
            input_folder = input('Specify the input directory in the correct format: ')

    return input_folder


def process_vcf(vcf_fn):
    in_vcf = VCF(vcf_fn)
    # count_pass = 0
    count_all = 0
    variant_table = []
    for v in in_vcf:
        # if v.FILTER == None or v.FILTER == '.':
        # count_pass += 1
        variant_name = v.CHROM + ';' + str(v.start + 1) + ';' + v.REF + ';' + v.ALT[0]

        # variant_info_all = ';'.join(['='.join([str(_) for _ in x]) for x in list(v.INFO)])
        variant_info_all = str(v).split('\t')[7].strip()
        # print('variant_info_all', variant_info_all)
        variant_info = variant_info_all

        # print(variant_name, variant_info)
        variant_table.append(
            {
                'variant_name': variant_name,
                'variant_info': variant_info,
            }
        )

        count_all += 1

    variant_table = pd.DataFrame(variant_table)
    in_vcf.close()
    return variant_table


def modify_target(target_str, source_str):
    # print('target_str', target_str)
    # print('source_str', source_str)

    target = sorted(target_str.split(';'))
    source = sorted(source_str.split(';'))

    target = [x.split('=') for x in target if x != '']
    source = [x.split('=') for x in source if x != '']

    # print('target', target)
    # print('source', source)

    target_names = [x[0] for x in target]
    source_names = [x[0] for x in source]
    diff_source = set(source_names).difference(set(target_names))
    diff_target = set(target_names).difference(set(source_names))
    isec = set(source_names).intersection(set(target_names))

    final_str = ''
    if len(diff_source) == 0:
        final_str = source_str
    else:
        for dt in diff_target:
            for t in target:
                if dt == t[0]:
                    # print('target original', dt, t[1])
                    final_str += dt + '=' + t[1] + ';'
        for i in isec:
            for si in source:
                if i == si[0]:
                    # print('target replaced', i, si[1])
                    final_str += i + '=' + si[1] + ';'
        for ds in diff_source:
            for sd in source:
                if ds == sd[0]:
                    # print('target extended', ds, sd[1])
                    final_str += ds + '=' + sd[1] + ';'
    if final_str[-1] == ';':
        final_str = final_str[:-1]
    print('\nfinal_str', final_str)
    return final_str


def write_vcf_as_table(raw_table, out_vcf_fn):
    curr_date = datetime.now().strftime("%A, %d. %B %Y %I:%M%p")
    comment = """##fileformat=VCFv4.2\n##fileDate=""" + curr_date
    comment += """\n##source=vc_process.py
##INFO=<ID=VD_Q,Number=1,Type=Float,Description="VarDict quality">
##INFO=<ID=VD_FILTER,Number=1,Type=Integer,Description="VarDict FILTER (PASS/excuse)">
##INFO=<ID=SINVICT,Number=1,Type=Integer,Description="Sinvict discovered the variant or not">
##INFO=<ID=SC_FILTER,Number=1,Type=String,Description="Scalpel FILTER (PASS/excuse)">
##INFO=<ID=SC_CHI2,Number=1,Type=Float,Description="Scalpel - k-mer Chi-Square score">
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
##INFO=<ID=Strelka_TNP_FILTER,Number=1,Type=String,Description="Filter results in tumor normal paires for strelka tnp">
##INFO=<ID=Mutect2_TNP_FILTER,Number=1,Type=String,Description="Filter results in tumor normal paires">\n"""
    ##INFO=<ID=SC_MINCOV,Number=1,Type=Integer,Description="Scalpel - (minimum) k-mer coverage of non-reference allele">
    contig = ''
    for c in sorted(set(raw_table['#CHROM'].to_list())):
        contig += '##contig=<ID=' + c + '>\n'
    comment += contig

    # raw_table.to_csv(input_folder + '/' + 'variant.raw._INDEL_.vcf', sep='\t', header=True, index=False)
    with open(out_vcf_fn, 'a') as f:
        f.write(comment)
        raw_table.to_csv(f, sep='\t', header=True, index=False)


def find_tools(res='/pipeline/scripts/resources.csv'):
    from func import read_resources
    return read_resources(res)


def compress_vcf(out_vcf_fn):
    # bgzip = find_tools()['bgzip']
    # tabix = find_tools()['tabix']
    bgzip = 'bgzip'
    tabix = 'tabix'
    run("%s -f %s" % (bgzip, out_vcf_fn), shell=True)
    run("%s -f %s" % (tabix, out_vcf_fn + '.gz'), shell=True)

    if os.path.isfile(out_vcf_fn + '.gz'):
        print('Raw VCF:', out_vcf_fn + '.gz')


if __name__ == '__main__':
    # input_folder = '/home/aurora/Atlas/Illumina_pipeline_development/samples/join'
    input_folder = input_dialog()

    wildcard_fn = '*raw*SNV.vcf'
    wildcard_fn2 = wildcard_fn + '.gz'
    # glob_handle = glob.glob(input_folder + '/' + wildcard_fn)
    # glob_handle = glob.glob(input_folder + '/' + wildcard_fn2)
    glob_handle = glob.glob(input_folder + '/' + wildcard_fn) + glob.glob(input_folder + '/' + wildcard_fn2)

    out_vcf_fn = input_folder + '/' + 'variant.raw.SNV.JOIN.vcf'

    # Parse VCF files in a given input_folder
    tables = []
    count_main = 0
    for vcf_fn in glob_handle:
        print('VCF', str(count_main) + ':', vcf_fn)
        table = process_vcf(vcf_fn)  # конкретная таблица
        print('\ntable\n', table)
        tables.append(table)
        count_main += 1

    # Extract all variant names
    variant_names = [t['variant_name'].to_list() for t in tables]
    variant_names = [i for sub in variant_names for i in sub]
    variant_names = sorted(list(set(variant_names)))
    print('variant_names', type(variant_names), variant_names)

    # Compose a raw table with all variants and all new INFO fields
    raw_table = pd.DataFrame(variant_names, columns=['variant_name'])
    print('\nraw_table\n', raw_table)
    raw_table['variant_info'] = None
    for raw_index in range(len(raw_table.index)):
        print('raw_index', raw_index, raw_table.iloc[raw_index]['variant_name'])
        for t in tables:
            for t_index in range(len(t.index)):
                if t.iloc[t_index, t.columns.get_loc('variant_name')] == raw_table.iloc[raw_index, raw_table.columns.get_loc('variant_name')]:
                    variant_info = t.iloc[t_index, t.columns.get_loc('variant_info')]
                    raw_info = raw_table.iloc[raw_index, raw_table.columns.get_loc('variant_info')]
                    if raw_info is None:
                        raw_table.iloc[raw_index, raw_table.columns.get_loc('variant_info')] = variant_info
                    else:
                        # raw_table.iloc[raw_index, raw_table.columns.get_loc('variant_info')] = raw_info + variant_info
                        raw_table.iloc[raw_index, raw_table.columns.get_loc('variant_info')] = modify_target(raw_info, variant_info)
    print('\nraw_table\n', raw_table)

    # Preparing the raw table for the output, splitting & adding & moving columns
    raw_table[['#CHROM', 'POS', 'REF', 'ALT']] = raw_table['variant_name'].str.split(';', expand=True)
    raw_table['ID'] = '.'
    raw_table['QUAL'] = '.'
    raw_table['FILTER'] = '.'
    raw_table['FORMAT'] = '.'
    raw_table[input_folder.split('/')[-2]] = '.'
    raw_table.rename(columns={'variant_info': 'INFO'}, inplace=True)
    raw_table = raw_table[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', input_folder.split('/')[-2]]]
    pd.set_option('display.max_columns', None)
    print('\nraw_table\n', raw_table)

    # Writing the table with a comment section (VCF)
    if os.path.isfile(out_vcf_fn):
        os.remove(out_vcf_fn)
    write_vcf_as_table(raw_table, out_vcf_fn)
    compress_vcf(out_vcf_fn)
