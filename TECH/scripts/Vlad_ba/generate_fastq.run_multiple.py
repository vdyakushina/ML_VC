import glob
import pandas as pd
import numpy as np
import sys
import os
import time
from subprocess import run


def input_dialog():
    if len(sys.argv) == 1:
        xls_fn = input('Specify the input Excel file: ')
        while not os.path.isfile(xls_fn):
            xls_fn = input('Specify the input Excel file that exists: ')
        samples_dir = input('Specify the input samples directory (where the individual sample dirs are): ')
        while not os.path.isdir(samples_dir):
            samples_dir = input('Specify the input samples directory that exists: ')

    else:
        xls_fn = sys.argv[1]
        samples_dir = sys.argv[2]

        while not os.path.isfile(xls_fn):
            xls_fn = input('Specify the input Excel file that exists: ')
        while not os.path.isdir(samples_dir):
            samples_dir = input('Specify the input samples directory that exists: ')

    return xls_fn, samples_dir


def input_path():
    from func import read_resources
    res = "/pipeline/scripts/resources.csv"
    path = read_resources(res)
    return path


if __name__ == '__main__':
    print('USAGE:', 'python generate_fastq.run_multiple.changed_by_mistake.py path/to/xls')
    xls_fn, samples_dir = input_dialog()

    path = input_path()

    variant_name = xls_fn.split('/')[-1].split('.')[0]
    variant_name = ':'.join(variant_name.split('_')[:2]) + '>' + variant_name.split('_')[2]
    print('variant_name', variant_name)

    xls = pd.read_excel(xls_fn)
    print('\nxls\b', xls)
    samples = xls['Positive'].dropna().to_list()
    sample_paths = [samples_dir + '/' + str(sample) + '/' + 'analysis/provisional.bam' for sample in samples]

    print('samples')
    for sample in sample_paths:
        sample_name = sample.split('/')[-3]
        print(str(sample_paths.index(sample) + 1) + '/' + str(len(sample_paths)), ':', sample_name)
        print("\t%s/generate_fastq_positive.py %s %s %s" % (path['scripts'], variant_name, sample, '/home/gkhvorykh/train_positive/positive_fastq'))
        run("%s/generate_fastq_positive.py %s %s %s" % (path['scripts'], variant_name, sample, '/home/gkhvorykh/train_positive/positive_fastq'), shell=True)
