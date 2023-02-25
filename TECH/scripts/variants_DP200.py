import glob
import numpy as np
from cyvcf2 import VCF, Writer
import pandas as pd
import sys
import os
import time
from subprocess import run


# TODO Формируем список образцов, которые покрывают этот вариант с DP более 200х:
#  - Для каждого образца прогоняем samtools depth по списку вариантов
#  - Для каждого варианта из списка-таблицы ловим глубину в каждом образце, если больше 200 то добавляем в колонку


def input_dialog():
    if len(sys.argv) == 1:
        samples_dir = input('Specify the input samples directory (where the individual sample dirs are): ')
        regions_fn = input('Specify the regions file path (TXT for samtools): ')
        variant_table_fn = input('Specify the variant table path (TSV): ')

    else:
        samples_dir = sys.argv[1]
        regions_fn = sys.argv[2]
        variant_table_fn = sys.argv[3]

    return samples_dir, regions_fn, variant_table_fn


def input_samtools():
    from func import read_resources
    res = "/pipeline/scripts/resources.csv"
    path = read_resources(res)
    return path['samtools']


def add_sample_name(x, depth_table, sample_name):
    depth = depth_table.iloc[x.index_fun][2]
    # print('depth', depth)

    # prev_sample_names = x['samples_from_samtools']
    # if prev_sample_names == '':

    if depth >= 200:
        return x['samples_from_samtools'] + sample_name + ','
    else:
        return x['samples_from_samtools'] + ''


def remove_last_char(final_str):
    if final_str[-1] == ',':
        final_str = final_str[:-1]
    return final_str


def modify_target(target_str, source_str):
    # print('target_str', target_str)
    # print('source_str', source_str)

    target = sorted(target_str.split(';'))
    source = sorted(source_str.split(';'))

    diff_source = set(source).difference(set(target))
    diff_target = set(target).difference(set(source))
    isec = set(source).intersection(set(target))

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
    # print('\nfinal_str', final_str)
    if final_str[-1] == ';':
        final_str = final_str[:-1]
    return final_str


if __name__ == '__main__':
    print('USAGE:', 'python variants_DP200.py path/to/samples_dir path/to/regions path/to/variant_table')

    samples_dir, regions_fn, variant_table_fn = input_dialog()
    variant_table = pd.read_table(variant_table_fn)
    variant_table['index_fun'] = variant_table.index.astype(int)
    variant_table['samples_from_samtools'] = ''

    glob_handle = glob.glob(samples_dir + '/*/' + 'analysis/provisional.bam')
    for sample in glob_handle:
        sample_name = sample.split('/')[-3]
        print(str(glob_handle.index(sample) + 1) + '/' + str(len(glob_handle)), ':', sample_name)

        regions_depth_fn = '/pipeline/scripts/depth_for_all_samples/' + 'regions.depth.' + sample_name + '.txt'
        run("%s depth -a -b %s %s > %s" % (input_samtools(), regions_fn, sample, regions_depth_fn), shell=True)

        depth_table = pd.read_table(regions_depth_fn, header=None)
        variant_table['samples_from_samtools'] = variant_table.apply(lambda x: add_sample_name(x, depth_table, sample_name), axis=1)
        print('\nvariant_table\n', variant_table)

    variant_table['samples_from_samtools'] = variant_table.apply(lambda x: remove_last_char(x['samples_from_samtools']), axis=1)
    variant_table.drop(['index_fun'], axis=1, inplace=True)
    variant_table.to_csv(variant_table_fn[:-4] + '.samples_from_samtools.tsv', sep='\t', index=False, header=True)
