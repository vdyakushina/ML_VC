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


def input_dialog():
    if len(sys.argv) == 1:
        pos_fn = input('Specify the input POSITIVE file: ')
        while not os.path.isfile(pos_fn):
            pos_fn = input('Specify the input POSITIVE file that exists: ')
        neg_fn = input('Specify the input NEGATIVE file: ')
        while not os.path.isfile(neg_fn):
            neg_fn = input('Specify the input NEGATIVE file that exists: ')
        samples_dir = input('Specify the input samples directory (where the individual sample dirs are): ')
        while not os.path.isdir(samples_dir):
            samples_dir = input('Specify the input samples directory that exists: ')

    else:
        pos_fn = sys.argv[1]
        neg_fn = sys.argv[2]
        samples_dir = sys.argv[3]

        while not os.path.isfile(pos_fn):
            pos_fn = input('Specify the input POSITIVE file that exists: ')
        while not os.path.isfile(neg_fn):
            neg_fn = input('Specify the input NEGATIVE file that exists: ')
        while not os.path.isdir(samples_dir):
            samples_dir = input('Specify the input samples directory that exists: ')

    return pos_fn, neg_fn, samples_dir


def input_path():
    from func import read_resources
    res = "/pipeline/scripts/resources.csv"
    path = read_resources(res)
    return path


def generate_fq_for_variant(table):
    time_start = time.time()
    variant_name = table['variant_name']
    variant_name_us = variant_name.replace(':', '_').replace('>', '_')
    out_fqbam_dir = samples_dir + '/' + variant_name_us + '/fqbam/'
    out_raw_dir = samples_dir + '/' + variant_name_us + '/raw/'
    Path(out_fqbam_dir).mkdir(parents=True, exist_ok=True)
    Path(out_raw_dir).mkdir(parents=True, exist_ok=True)

    output_sam = out_fqbam_dir + variant_name_us + ".nonum.sequences.sam"
    output_bam = output_sam[:-4] + ".bam"
    output_bam_sorted = output_bam[:-4] + ".sorted.bam"
    output_bam_raw = out_raw_dir + variant_name_us + ".nonum.raw.bam"
    # output_sam = "'" + output_sam + "'"
    # output_bam = "'" + output_bam + "'"
    # output_bam_sorted = "'" + output_bam_sorted + "'"
    
    # pos_samples = table['positive'].str.split(',')
    # neg_samples = table['negative'].str.split(',')
    # print('\nWorking on variant:', variant_name, 'with \n\t\tpositive samples', pos_samples, '\n\t\tnegative samples', neg_samples)

    adopted_reads_fn = '/pipeline/data/samples_tsv/chr11_108098525G_A/raw_d0cNZaUxL5d7L/d0cNZaUxL5d7L_R1.fq'
    print('\nAdopting reads from', adopted_reads_fn)
    adopted_reads = list(SeqIO.parse(adopted_reads_fn, "fastq"))

    this_positive_R1 = pos_fq_dir + variant_name + pos_R1_fq
    this_positive_R2 = pos_fq_dir + variant_name + pos_R2_fq
    this_negative_R1 = neg_fq_dir + variant_name + neg_R1_fq
    this_negative_R2 = neg_fq_dir + variant_name + neg_R2_fq
    print('\nWorking with variant', variant_name, 'files of interest:')
    print('\t', this_positive_R1, this_positive_R2, '\n\t', this_negative_R1, this_negative_R2)

    common_pos_records1 = list(SeqIO.parse(this_positive_R1, "fastq"))
    common_pos_records2 = list(SeqIO.parse(this_positive_R2, "fastq"))
    common_neg_records1 = list(SeqIO.parse(this_negative_R1, "fastq"))
    common_neg_records2 = list(SeqIO.parse(this_negative_R2, "fastq"))

    f = [0.05, 0.1, 0.25, 0.5]
    dp = [(500, 2001), (500, 2001), (250, 1501), (250, 1501)]

    f = f[0]
    dp = dp[0]

    f = 0.25
    dp = (300, 2000)
    print('Chosen f', f, '\nChosen DP range', dp)

    # Checking status
    ok = []

    for n in range(1, 2):  # 11
        time_start_loop = time.time()
        sequences_R1 = []
        sequences_R2 = []
        print('N:', n)
        # print('\t', [np.random.randint(x[0], x[1]) for x in dp])
        random_DP = range(1, np.random.randint(dp[0], dp[1]))  # range(1, np.random.randint(dp[0], dp[1])) range(1, 11)
        # random_DP = range(1, 100)
        print('random_DP', random_DP[-1])
        for i in random_DP:
            pos_records1 = list(common_pos_records1)
            pos_records2 = list(common_pos_records2)
            neg_records1 = list(common_neg_records1)
            neg_records2 = list(common_neg_records2)
            F = np.random.randint(1, 10000)
            if F < 10000*f:
                # print('positive')
                records1 = []
                records2 = []
                read_number = np.random.randint(1, len(pos_records1))  # min(len(pos_records1), len(pos_records2))
                print('\n\t', i, '/', random_DP[-1], ':', 'read_number', read_number,
                      'len(R1)', len(pos_records1), 'len(R2)', len(pos_records2))

                sought_id = pos_records1[read_number].id
                records1 = copy.deepcopy([x for x in pos_records1 if x.id == sought_id])  # copy.deepcopy(pos_records1)
                records2 = copy.deepcopy([x for x in pos_records2 if x.id == sought_id])  # copy.deepcopy(pos_records2)
                # print('\tnew len(R1)', len(records1), records1)
                # print('\tnew len(R2)', len(records2), records2)

                if records2:
                    print('\tF', F, 'SIMULATED POSITIVE R1', records1[0].id, 'R2', records2[0].id)
                    records1[0].id = adopted_reads[i].id
                    records2[0].id = adopted_reads[i].id
                    sequences_R1.append(records1[0])
                    sequences_R2.append(records2[0])

                    # random_id = str(np.random.randint(100000, 999999))
                    # records1[0].id = ':'.join(records1[0].id.split(':')[:-1]) + ':' + random_id
                    # records2[0].id = ':'.join(records2[0].id.split(':')[:-1]) + ':' + random_id
                    # records1[0].name = ':'.join(records1[0].name.split(':')[:-1]) + ':' + random_id
                    # records2[0].name = ':'.join(records2[0].name.split(':')[:-1]) + ':' + random_id
                    # records1[0].description = ':'.join(records1[0].description.split(':')[:-1]) + ':' + random_id
                    # records2[0].description = ':'.join(records2[0].description.split(':')[:-1]) + ':' + random_id
                    # records1[0].id += ':' + random_id  # '/1'
                    # records2[0].id += ':' + random_id  # '/2'
                    # TODO remove reads (read names) that occur more than twice
                    # print('\tlen([x.id for x in sequences_R1])', len([x.id for x in sequences_R1]),
                    #       'len(set([x.id for x in sequences_R1]))', len(set([x.id for x in sequences_R1])))
                    # if records1[0].id not in [x.id for x in sequences_R1]:
                    #     sequences_R1.append(records1[0])
                    # if records2[0].id not in [x.id for x in sequences_R2]:
                    #     sequences_R2.append(records2[0])
                    print('\tF', F, 'ADOPTED POSITIVE R1', records1[0].id, 'R2', records2[0].id)
                else:
                    print('\tOops the sought id', sought_id, 'was not found in R2')
                    pass
            else:
                # print('negative')
                records1 = []
                records2 = []
                read_number = np.random.randint(1, len(neg_records1))
                print('\n\t', i, '/', random_DP[-1], ':', 'read_number', read_number,
                      'len(R1)', len(neg_records1), 'len(R2)', len(neg_records2))

                sought_id = neg_records1[read_number].id
                records1 = copy.deepcopy([x for x in neg_records1 if x.id == sought_id])  # copy.deepcopy(neg_records1)
                records2 = copy.deepcopy([x for x in neg_records2 if x.id == sought_id])  # copy.deepcopy(neg_records2)
                # print('\tnew len(R1)', len(records1), records1)
                # print('\tnew len(R2)', len(records2), records2)

                if records2:
                    print('\tF', F, 'SIMULATED NEGATIVE R1', records1[0].id, 'R2', records2[0].id)
                    records1[0].id = adopted_reads[i].id
                    records2[0].id = adopted_reads[i].id
                    sequences_R1.append(records1[0])
                    sequences_R2.append(records2[0])

                    # random_id = str(np.random.randint(100000, 999999))
                    # records1[0].id = ':'.join(records1[0].id.split(':')[:-1]) + ':' + random_id
                    # records2[0].id = ':'.join(records2[0].id.split(':')[:-1]) + ':' + random_id
                    # records1[0].name = ':'.join(records1[0].name.split(':')[:-1]) + ':' + random_id
                    # records2[0].name = ':'.join(records2[0].name.split(':')[:-1]) + ':' + random_id
                    # records1[0].description = ':'.join(records1[0].description.split(':')[:-1]) + ':' + random_id
                    # records2[0].description = ':'.join(records2[0].description.split(':')[:-1]) + ':' + random_id
                    # records1[0].id += ':' + random_id  # '/1'
                    # records2[0].id += ':' + random_id  # '/2'
                    # TODO remove reads (read names) that occur more than twice
                    # print('\tlen([x.id for x in sequences_R1])', len([x.id for x in sequences_R1]),
                    #       'len(set([x.id for x in sequences_R1]))', len(set([x.id for x in sequences_R1])))
                    # if records1[0].id not in [x.id for x in sequences_R1]:
                    #     sequences_R1.append(records1[0])
                    # if records2[0].id not in [x.id for x in sequences_R2]:
                    #     sequences_R2.append(records2[0])
                    print('\tF', F, 'ADOPTED NEGATIVE R1', records1[0].id, 'R2', records2[0].id)
                else:
                    print('\tOops the sought id', sought_id, 'was not found in R2')
                    pass

        print('\nWriting forward (R1) and reverse (R2) FASTQ sequences')
        sequences_R1_fn = out_fqbam_dir + variant_name_us + '.' + str(n) + '.sequences_R1.fastq'
        sequences_R2_fn = out_fqbam_dir + variant_name_us + '.' + str(n) + '.sequences_R2.fastq'
        with open(sequences_R1_fn, "w") as output_handle:
            SeqIO.write(sequences_R1, output_handle, "fastq")
            print('\t', len(sequences_R1), 'records written to', sequences_R1_fn)
        with open(sequences_R2_fn, "w") as output_handle:
            SeqIO.write(sequences_R2, output_handle, "fastq")
            print('\t', len(sequences_R2), 'records written to', sequences_R2_fn)

        # print('\nWriting a single paired end FASTQ file')
        # output_fn = samples_dir + '/' + variant_name + '.' + str(n) + '.sequences_interleaved.fastq'
        # output_handle = open(output_fn, "w")
        # count = 0
        # f_iter = FastqGeneralIterator(open(sequences_R1_fn, 'r'))
        # r_iter = FastqGeneralIterator(open(sequences_R2_fn, 'r'))
        # for (f_id, f_seq, f_q), (r_id, r_seq, r_q) in zip(f_iter, r_iter):
        #     assert f_id == r_id
        #     count += 2
        #     # Write out both reads with "/1" and "/2" suffix on ID
        #     output_handle.write("@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n"
        #                         % (f_id, f_seq, f_q, r_id, r_seq, r_q))
        # output_handle.close()
        # print("\t%s records written to %s" % (count, output_fn))

        # Running bwa and samtools to align, sort, index the sequences
        output_sam = out_fqbam_dir + variant_name_us + "." + str(n) + ".sequences.sam"
        output_bam = output_sam[:-4] + ".bam"
        output_bam_sorted = output_bam[:-4] + ".sorted.bam"
        output_bam_raw = out_raw_dir + variant_name_us + "." + str(n) + ".raw.bam"
        # output_sam = "'" + output_sam + "'"
        # output_bam = "'" + output_bam + "'"
        # output_bam_sorted = "'" + output_bam_sorted + "'"
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
        time_end_loop = time.time() - time_start_loop
        print('Loop', n, 'for the variant', variant_name, 'took {:.3f} s'.format(time_end_loop) +
              ' or {:.3f} m'.format(time_end_loop / 60.0) + ' or {:.3f} h'.format(time_end_loop / 3600.0))

    time_end = time.time() - time_start
    print('Making FASTQs and alignments for the variant', variant_name, 'took {:.3f} s'.format(time_end) +
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
    print('USAGE:', 'python generate_fastq_from_tsv_1j.py path/to/pos.tsv path/to/neg.tsv path/to/output/samples')

    pos_fn, neg_fn, samples_dir = input_dialog()
    pos_fq_dir = '/'.join(pos_fn.split('/')[:-1]) + '/positive_fastq/'
    neg_fq_dir = '/'.join(neg_fn.split('/')[:-1]) + '/negative_fastq/'
    pos_R1_fq = '.positive.R1.fq'
    pos_R2_fq = '.positive.R2.fq'
    neg_R1_fq = '.negative.R1.fq'
    neg_R2_fq = '.negative.R2.fq'

    path = input_path()

    pos = pd.read_table(pos_fn, header=None)
    neg = pd.read_table(neg_fn, header=None)
    print('\npos\n', pos)
    print('\nneg\n', neg)

    if pos[0].to_list() == neg[0].to_list():
        pos_neg = pos.merge(neg, on=0)
        pos_neg.columns = ['variant_name', 'positive', 'negative']
        print('\npos_neg\n', pos_neg)
        pos_neg['status'] = pos_neg.apply(lambda x: generate_fq_for_variant(x), axis=1)
        # pos_neg['status_positive'] = pos_neg.apply(lambda x: generate_fq_for_variant(x['variant_name'], x['positive'], pos_fq_dir), axis=1)
        # pos_neg['status_negative'] = pos_neg.apply(lambda x: generate_fq_for_variant(x['variant_name'], x['negative'], neg_fq_dir), axis=1)
        print('\npos_neg\n', pos_neg)
