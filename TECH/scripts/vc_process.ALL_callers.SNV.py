import glob
from cyvcf2 import VCF, Writer
import pandas as pd
from datetime import datetime
import sys
import os

# INPUT ARGUMENT: analysis dir, search for other sample_name.*.SNV.vcf and sample_name.*.SNV.vcf
# PARSE all these files and merge all variants, adding new INFO fields by variant callers
# CREATE tables / VCFs sample_name.raw.SNV.vcf and sample_name.raw.SNV.vcf with the merged variants
# -- check if raw VCFs exist first

# План - сканировать папку на предмет наличия этого файла и файлов коллеров,
# отдельно вытаскивать все варианты из последних и добавлять в общую таблицу с пометкой коллера.
# Вывод таблицы в vcf.


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
        variant_info = None
        variant_filter = 'PASS' if v.FILTER == None else str(v.FILTER)
        variant_filter = variant_filter.replace(';', ',')
        # if 'mutect2_to' in vcf_fn:
        #     variant_info = 'MUTECT2_TO_Q=' + v.QUAL + ';' + 'Mutect2_TO_FILTER=' + v.FILTER
        if 'vardict' in vcf_fn:
            variant_info = 'VD_Q=' + str(v.QUAL) + ';' + 'VD_FILTER=' + variant_filter + ';'
        elif 'sinvict' in vcf_fn:
            variant_info = 'SINVICT=1' + ';'
        elif 'scalpel' in vcf_fn:
            variant_info = 'SC_FILTER=' + variant_filter + ';' + \
                           'SC_CHI2=' + str(v.INFO.get('CHI2')) + ';' + 'SC_MINCOV=' + str(v.INFO.get('MINCOV')) + ';'

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


if __name__ == '__main__':
    # input_folder = '/home/aurora/Atlas/Illumina_pipeline_development/samples/bH43SnB9yVka9/analysis'
    input_folder = input_dialog()
    wildcard_fn = '*SNV.vcf'

    # Parse VCF files in a given input_folder
    tables = []
    count_main = 0
    for vcf_fn in glob.glob(input_folder + '/' + wildcard_fn):
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
                        raw_table.iloc[raw_index, raw_table.columns.get_loc('variant_info')] = raw_info + variant_info
    raw_table['variant_info'] = raw_table['variant_info'].str[:-1]
    print('\nraw_table\n', raw_table)

    # Preparing the raw table for the output, splitting & adding & moving columns
    raw_table[['#CHROM', 'POS', 'REF', 'ALT']] = raw_table['variant_name'].str.split(';', expand=True)
    raw_table['ID'] = '.'
    raw_table['QUAL'] = '.'
    raw_table['FILTER'] = '.'
    raw_table[input_folder.split('/')[-2]] = '.'
    raw_table.rename(columns={'variant_info': 'INFO'}, inplace=True)
    raw_table = raw_table[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', input_folder.split('/')[-2]]]
    pd.set_option('display.max_columns', None)
    print('\nraw_table\n', raw_table)

    # Writing the table with a comment section (VCF)
    curr_date = datetime.now().strftime("%A, %d. %B %Y %I:%M%p")
    comment = """##fileformat=VCFv4.1 
##fileDate=""" + curr_date
    comment += """\n##source=vc_process.py
##INFO=<ID=VD_Q,Number=1,Type=Float,Description="VarDict quality">
##INFO=<ID=VD_FILTER,Number=1,Type=Integer,Description="VarDict FILTER (PASS/excuse)">
##INFO=<ID=SINVICT,Number=1,Type=Integer,Description="Sinvict discovered the variant or not">
##INFO=<ID=SC_FILTER,Number=1,Type=String,Description="Scalpel FILTER (PASS/excuse)">\n"""
    contig = ''
    for c in sorted(set(raw_table['#CHROM'].to_list())):
        contig += '##contig=<ID=' + c + '>\n'
    comment += contig

    # raw_table.to_csv(input_folder + '/' + 'variant.raw._INDEL_.vcf', sep='\t', header=True, index=False)
    with open(input_folder + '/' + 'variant.raw.' + wildcard_fn.replace('*', '_').replace('.vcf', '_.vcf'), 'a') as f:
        f.write(comment)
        raw_table.to_csv(f, sep='\t', header=True, index=False)
