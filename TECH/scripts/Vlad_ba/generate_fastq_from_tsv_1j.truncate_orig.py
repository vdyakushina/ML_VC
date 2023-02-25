import glob
import pandas as pd
import numpy as np
import sys
import os
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools
from subprocess import run
import copy
import shutil
from pathlib import Path
from operator import itemgetter, attrgetter


def input_path():
    from func import read_resources
    res = "/pipeline/scripts/resources.csv"
    path = read_resources(res)
    return path


def generate_truncated_fqbam():
    time_start = time.time()

    orig_R1_fn = '/pipeline/data/samples_tsv/chr11_108098525G_A/raw_d0cNZaUxL5d7L/d0cNZaUxL5d7L_R1.fq'
    print('\nImporting reads from', orig_R1_fn)
    orig_R1 = list(SeqIO.parse(orig_R1_fn, "fastq"))

    orig_R2_fn = '/pipeline/data/samples_tsv/chr11_108098525G_A/raw_d0cNZaUxL5d7L/d0cNZaUxL5d7L_R2.fq'
    print('\nImporting reads from', orig_R2_fn)
    orig_R2 = list(SeqIO.parse(orig_R2_fn, "fastq"))

    target_R1_fn = '/pipeline/data/samples_tsv/chr11_108098525G_A/fqbam/chr11_108098525G_A.1.sequences_R1.fastq'
    print('\nImporting reads from', target_R1_fn)
    target_R1 = list(SeqIO.parse(target_R1_fn, "fastq"))

    target_R2_fn = '/pipeline/data/samples_tsv/chr11_108098525G_A/fqbam/chr11_108098525G_A.1.sequences_R2.fastq'
    print('\nImporting reads from', target_R2_fn)
    target_R2 = list(SeqIO.parse(target_R2_fn, "fastq"))

    out_fqbam_dir = '/pipeline/data/samples_tsv/chr11_108098525G_A/fqbam_d0cNZaUxL5d7L_s/'
    out_raw_dir = '/pipeline/data/samples_tsv/chr11_108098525G_A/raw_d0cNZaUxL5d7L_s/'

    sequences_R1 = []
    sequences_R2 = []

    for records1 in orig_R1:
        if records1.id in [x.id for x in target_R1]:
            sequences_R1.append(records1)
            print('R1:', records1.id, len(sequences_R1), '/', len(target_R1))
        if len(sequences_R1) == len(target_R1):
            break

    for records2 in orig_R2:
        if records2.id in [x.id for x in target_R2]:
            sequences_R2.append(records2)
            print('R2:', records2.id, len(sequences_R2), '/', len(target_R2))
        if len(sequences_R2) == len(target_R2):
            break

    print('Sorting R1')
    # sorted_id_R1 = sorted([x.id for x in sequences_R1])
    sequences_R1 = sorted(sequences_R1, key=attrgetter('id'))
    print('Sorting R2')
    # sorted_id_R2 = sorted([x.id for x in sequences_R2])
    sequences_R2 = sorted(sequences_R2, key=attrgetter('id'))

    print('\nWriting forward (R1) and reverse (R2) FASTQ sequences')
    sequences_R1_fn = out_fqbam_dir + 'd0cNZaUxL5d7L_s.sequences_R1.fastq'
    sequences_R2_fn = out_fqbam_dir + 'd0cNZaUxL5d7L_s.sequences_R2.fastq'
    with open(sequences_R1_fn, "w") as output_handle:
        SeqIO.write(sequences_R1, output_handle, "fastq")
        print('\t', len(sequences_R1), 'records written to', sequences_R1_fn)
    with open(sequences_R2_fn, "w") as output_handle:
        SeqIO.write(sequences_R2, output_handle, "fastq")
        print('\t', len(sequences_R2), 'records written to', sequences_R2_fn)


    # Running bwa and samtools to align, sort, index the sequences
    output_sam = out_fqbam_dir + "d0cNZaUxL5d7L_s.sequences.sam"
    output_bam = output_sam[:-4] + ".bam"
    output_bam_sorted = output_bam[:-4] + ".sorted.bam"
    output_bam_raw = out_raw_dir + "d0cNZaUxL5d7L_s.raw.bam"
    # output_sam = "'" + output_sam + "'"
    # output_bam = "'" + output_bam + "'"
    # output_bam_sorted = "'" + output_bam_sorted + "'"
    ok = []
    if os.path.isfile(sequences_R1_fn) and os.path.isfile(sequences_R2_fn):
        ok.append(True)
        # sequences_R1_fn = "'" + sequences_R1_fn + "'"
        # sequences_R2_fn = "'" + sequences_R2_fn + "'"
        print('\nAlignment:')
        print("%s mem -t 6 %s %s %s > %s" % (path['bwa'], path['ref'], sequences_R1_fn, sequences_R2_fn, output_sam))
        run("%s mem -t 6 %s %s %s > %s" % (path['bwa'], path['ref'], sequences_R1_fn, sequences_R2_fn, output_sam), shell=True)
        if os.path.isfile(output_sam):  # [1:-1]
            ok.append(True)
            print("\n%s view -@ 6 -b %s > %s" % (path['samtools'], output_sam, output_bam))
            run("%s view -@ 6 -b %s > %s" % (path['samtools'], output_sam, output_bam), shell=True)
            if os.path.isfile(output_bam):
                ok.append(True)
                print("\n%s sort -@ 6 %s > %s" % (path['samtools'], output_bam, output_bam_sorted))
                run("%s sort -@ 6 %s > %s" % (path['samtools'], output_bam, output_bam_sorted), shell=True)
                if os.path.isfile(output_bam_sorted):
                    ok.append(True)
                    # os.remove(output_sam)
                    # os.remove(output_bam)
                    print("\n%s index -@ 6 %s" % (path['samtools'], output_bam_sorted))
                    run("%s index -@ 6 %s" % (path['samtools'], output_bam_sorted), shell=True)
                    print('\nAlignment written to', output_bam_sorted)
                else:
                    ok.append(False)
                    print('ERROR:', output_bam_sorted, 'does not exist')
            else:
                ok.append(False)
                print('ERROR:', output_bam, 'does not exist')
        else:
            ok.append(False)
            print('ERROR:', output_sam, 'does not exist')
    else:
        ok.append(False)
        print('ERROR:', sequences_R1_fn, 'or', sequences_R2_fn, 'does not exist')

    time_end = time.time() - time_start
    print('Making FASTQs and alignments took {:.3f} s'.format(time_end) +
          ' or {:.3f} m'.format(time_end / 60.0) + ' or {:.3f} h'.format(time_end / 3600.0))
    if all(ok):
        print('SUCCESS')
        shutil.copy2(output_bam_sorted, output_bam_raw)
        shutil.copy2(output_bam_sorted + '.bai', output_bam_raw + '.bai')
        return 'SUCCESS'
    else:
        print('FAIL')
        return 'FAIL'


if __name__ == '__main__':
    path = input_path()
    generate_truncated_fqbam()
