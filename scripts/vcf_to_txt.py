import glob
import numpy as np
from cyvcf2 import VCF, Writer
import pandas as pd
import sys
import os
import time


def input_dialog():
    if len(sys.argv) == 1:
        in_vcf_folder = input('Specify the input VCF wildcard (where to look for the files): ')
        out_txt_fn = input('Specify the output TSV table path: ')

    else:
        in_vcf_folder = sys.argv[1]
        out_txt_fn = sys.argv[2]

    return in_vcf_folder, out_txt_fn


def input_db_server():
    from func import read_resources
    res = "/pipeline/scripts/resources.csv"
    path = read_resources(res)
    exac_fn = path['exac']
    gnomad_genome_fn = path['gnomad_genome_abc']
    gnomad_exome_fn = path['gnomad_exome_abc']
    return exac_fn, gnomad_genome_fn, gnomad_exome_fn


def create_list_of_samples(variant_name, mega_table):
    mega_table = mega_table[mega_table['variant_name'] == variant_name]
    return ','.join(mega_table['sample_name'].to_list())


def get_gnomad_byvar_local(variant, gnomad_genome_fn):
    time_start = time.time()
    # variant = CHROMOSOME:POSITION:REFERENCE:VARIANT
    # print('\nExtracting MAF (gnomAD) for', variant)
    variant = variant.split(':')
    chrom = variant[0]
    pos = variant[1]
    ref = variant[2]
    alt = variant[3]

    gnomad = pd.read_table(gnomad_genome_fn)

    # gnomad = gnomad.loc[(gnomad['#Chr'] == chrom) and
    #                     (gnomad['Start'] == pos) and
    #                     (gnomad['Ref'] == ref) and
    #                     (gnomad['Alt'] == alt)]
    # gnomad = gnomad[(gnomad.loc[:, '#Chr'] == chrom[3:]) &
    #                 (gnomad.loc[:, 'Start'] == pos) &
    #                 (gnomad.loc[:, 'Ref'] == ref) &
    #                 (gnomad.loc[:, 'Alt'] == alt)]
    # gnomad = gnomad[gnomad.loc[:, '#Chr'] == chrom[3:]]
    gnomad = gnomad[(gnomad['#Chr'].astype(str) == str(chrom[3:])) &
                    (gnomad['Start'].astype(str) == str(pos)) &
                    (gnomad['Ref'].astype(str) == ref) &
                    (gnomad['Alt'].astype(str) == alt)]

    if gnomad.empty:
        af = 0
    else:
        if 'AF_genome' in gnomad.columns:
            af = gnomad['AF_genome'].to_list()[0]
        else:
            af = gnomad['AF_exome'].to_list()[0]

    # print('MAF (gnomAD) for', ':'.join(variant), 'is', af)
    # print('This took {:.3f} s'.format(time.time() - time_start))
    return af


def get_exac_byvar_local(variant, exac_fn):
    time_start = time.time()
    # variant = CHROMOSOME:POSITION:REFERENCE:VARIANT
    # print('\nExtracting MAF (ExAC) for', variant)
    variant = variant.split(':')
    chrom = variant[0]
    pos = variant[1]
    ref = variant[2]
    alt = variant[3]

    # exac_vcf = VCF('/home/aurora/Atlas/Illumina_pipeline_development/DB/ExAC/' + 'ExAC.r1.sites.vep.' + chrom[3:] + '.vcf.gz')
    # chrom = chrom[3:]
    # exac_vcf = VCF('/'.join(exac_fn.split('/')[:-1]) + '/' + 'ExAC.r1.sites.vep.' + chrom + '.vcf.gz')
    exac_vcf = VCF('/'.join(exac_fn.split('/')[:-1]) + '/' + 'ExAC.regions.AODABCV1.vcf.gz')

    # ac = 0
    # an = 0
    af = 0
    for v in exac_vcf:
        if v.CHROM == chrom and str(v.start + 1) == pos and v.REF == ref:  # v.start + 1  chrom[3:]
            # ac = v.INFO.get('AC')
            # an = v.INFO.get('AN')
            af = v.INFO.get('AF')
            if type(af) == tuple and alt in v.ALT:
                alt_index = v.ALT.index(alt)
                # ac = ac[alt_index]
                af = af[alt_index]
            break
    exac_vcf.close()
    # print('MAF (ExAC) for', ':'.join(variant), 'is', af)
    # print('This took {:.3f} s'.format(time.time() - time_start))
    return af


def change_variant_name(variant_name):
    variant_name = variant_name.split(':')
    return ':'.join(variant_name[:2]) + '>'.join(variant_name[2:])


def process_vcf(vcf_fn):
    in_vcf = VCF(vcf_fn)
    count_all = 0
    count_pass = 0
    variant_table = []
    for v in in_vcf:
        # if v.INFO.get('VAF') > 0.2:  # SGA_DP SGA_VAF
        # if v.INFO.get('DP') > 500 and v.INFO.get('VAF') > 0.2:  # SGA_DP SGA_VAF
        variant_name = v.CHROM + ':' + str(v.start + 1) + ':' + v.REF + ':' + v.ALT[0]
        # variant_info_all = str(v).split('\t')[7].split(';')
        # print(variant_name, variant_info)

        # GT:AD:AF:DP:F1R2:F2R1:SB        0/1:412,435:0.513:847:193,223:219,209:245,167,261,174
        # SB > 0.5 (берем из .vcf файла mutect, F1R2($2/$1)/F2R1($2/$1) -> reverse, если больше 1)
        variant_DP = sum(v.format('AD')[0])  # v.format('AD')[1]
        variant_AF = v.format('AF')[0][0]
        f1r2 = v.format('F1R2')[0]
        f2r1 = v.format('F2R1')[0]
        if any(x == 0 for x in f1r2):
            f1r2 += 1
        if any(x == 0 for x in f2r1):
            f2r1 += 1
        variant_SB = (f1r2[1]/f1r2[0]) / (f2r1[1]/f2r1[0])
        if variant_SB > 1:
            variant_SB = (f2r1[1] / f2r1[0]) / (f1r2[1] / f1r2[0])

        # print('DP', variant_DP, 'AF', variant_AF, 'variant_SB', variant_SB)
        
        # DP (берем из .vcf файла mutect, сумма по AD, поле FORMAT) > 500; VAF (берем из .vcf файла mutect, поле AF в format) > 20%; SB > 0.5 (берем из .vcf файла mutect, F1R2($2/$1)/F2R1($2/$1) -> reverse, если больше 1); MAF (в соответствии с ExAC) < 1%. 

        if variant_DP > 500 and variant_AF > 0.2 and variant_SB > 0.5:
            variant_table.append(
                {
                    'sample_name': in_vcf.samples[0],
                    'variant_name': variant_name,
                    'DP': variant_DP,
                }
            )
            count_pass += 1
        count_all += 1
    print('\tPassed', count_pass, 'variants of', count_all, 'parsed')
    variant_table = pd.DataFrame(variant_table)
    if variant_table.empty:
        variant_table = pd.DataFrame(columns=['variant_name'])
    in_vcf.close()
    return variant_table


if __name__ == '__main__':
    """
    Из имеющихся .vcf файлов по всем 233 образцам вытаскиваем варианты по следующим критериям:
    DP > 500; VAF > 20%; MAF (в соответствии с ExAC) < 1%.
    Уникальные варианты записываем в TXT
    """

    print('USAGE:', 'python vcf_to_txt.py path/to/*/vcf path/to/txt')
    # exac_fn = '/home/aurora/Atlas/Illumina_pipeline_development/samples/exac.vcf'
    exac_fn, gnomad_genome_fn, gnomad_exome_fn = input_db_server()
    in_vcf_folder, out_txt_fn = input_dialog()
    glob_handle = glob.glob(in_vcf_folder + '*' + '/variant.mutect2_to.vcf.gz')  # *raw*.vcf.gz
    # glob_handle = glob.glob(in_vcf_folder + '/*raw*.vcf.gz')
    print(in_vcf_folder, glob_handle)

    # pd.set_option('display.max_rows', None)
    mega_table = pd.DataFrame(columns=['variant_name'])
    for in_vcf_fn in glob_handle:
        # Parse VCF files in a given in_vcf_fn
        print('\nVCF:', in_vcf_fn)
        table = process_vcf(in_vcf_fn)
        print('\tList of variants:', table['variant_name'].to_list())
        # mega_table = pd.concat([mega_table, table], ignore_index=True, join='inner')
        mega_table = pd.concat([mega_table, table], axis=0)
        # print(mega_table)
        # print(len(mega_table['variant_name']), len(set(mega_table['variant_name'])))

    mega_table.sort_values(by=['variant_name'], inplace=True)
    mega_table.reset_index(inplace=True, drop=True)
    print('\n', mega_table.shape[0], 'non-unique variants meeting the criteria: DP > 500 and VAF > 0.2\n', mega_table)
    unique_table = mega_table.copy().groupby(mega_table.variant_name.tolist(), as_index=False).size()
    # unique_table = mega_table.groupby(mega_table.columns.tolist(), as_index=False).size()
    # mega_table.drop_duplicates(inplace=True, subset=['variant_name'])
    unique_table = unique_table.to_frame()
    unique_table.reset_index(inplace=True)
    unique_table.columns = ['variant_name', 'n_samples']
    print('\n', unique_table.shape[0], 'unique variants meeting the criteria: DP > 500 and VAF > 0.2\n', unique_table, '\n')

    # # Removing long variants
    unique_table[['tmp_rest', 'tmp_ref', 'tmp_alt']] = unique_table.variant_name.str.rsplit(':', n=2, expand=True)
    unique_table['tmp_ref'] = unique_table['tmp_ref'].apply(lambda x: len(x))
    unique_table['tmp_alt'] = unique_table['tmp_alt'].apply(lambda x: len(x))
    unique_table.drop(unique_table[unique_table.tmp_ref > 20].index, inplace=True)
    unique_table.drop(unique_table[unique_table.tmp_alt > 20].index, inplace=True)
    unique_table.drop(['tmp_rest', 'tmp_ref', 'tmp_alt'], axis=1, inplace=True)
    print('\n', unique_table.shape[0], 'unique variants with length of REF or ALT <= 20\n', unique_table, '\n')

    # Adding MAF from ExAC and gnomAD
    # unique_table = unique_table[:5]
    print('Extracting MAF (ExAC)')
    unique_table['MAF_ExAC'] = unique_table.apply(lambda x: get_exac_byvar_local(x['variant_name'], exac_fn), axis=1)
    print('Extracting MAF (gnomAD)')
    unique_table['MAF_gnomAD_E'] = unique_table.apply(lambda x: get_gnomad_byvar_local(x['variant_name'], gnomad_exome_fn), axis=1)
    unique_table['MAF_gnomAD_G'] = unique_table.apply(lambda x: get_gnomad_byvar_local(x['variant_name'], gnomad_genome_fn), axis=1)
    print('\nVariants with MAF from ExAC and gnomAD\n', unique_table)

    unique_table['MAF_gnomAD_max'] = unique_table.apply(lambda x: max(x['MAF_gnomAD_E'], x['MAF_gnomAD_G']), axis=1)
    unique_table.drop(unique_table[unique_table.MAF_gnomAD_max >= 0.01].index, inplace=True)
    # unique_table.drop(['MAF_gnomAD_max'], axis=1, inplace=True)
    print('\nVariants with MAF from MAF_gnomAD (G/E max) < 0.01\n', unique_table)

    variants_table = unique_table.copy()[['variant_name']]
    variants_table[['chr', 'pos', 'ref', 'alt']] = variants_table.variant_name.str.split(':', expand=True)
    print('\nVariant coordinates only\n', variants_table[['chr', 'pos']])
    variants_table[['chr', 'pos']].to_csv(out_txt_fn[:-4] + '.samtools.txt', sep='\t', header=False, index=False)
    print('Variant coordinates written to', out_txt_fn[:-4] + '.samtools.txt')

    # TODO Формируем список образцов, которые покрывают каждый вариант с DP более 200х:
    #  добавляем DP в таблицу
    #  уникальные варианты отправляем в отдельную таблицу (копию)
    #  каждый уникальный враиант сравниваем с большой таблицей и при совпадении вытаскиваем образцы с DP > 200
    # mega_table['variant_name'] = mega_table.apply(lambda x: change_variant_name(x['variant_name']), axis=1)
    mega_table.drop(mega_table[mega_table.DP < 200].index, inplace=True)
    unique_table['samples_DP200'] = unique_table.apply(lambda x: create_list_of_samples(x['variant_name'], mega_table), axis=1)
    unique_table['variant_name'] = unique_table.apply(lambda x: change_variant_name(x['variant_name']), axis=1)
    print('\nVariants with samples with DP >= 200\n', unique_table)

    unique_table.variant_name.to_csv(out_txt_fn, header=False, index=False)
    unique_table.to_csv(out_txt_fn[:-4] + '.full.tsv', sep='\t', header=True, index=False)
    print('List of variants written to', out_txt_fn)
    print('Table of variants written to', out_txt_fn[:-4] + '.full.tsv')




