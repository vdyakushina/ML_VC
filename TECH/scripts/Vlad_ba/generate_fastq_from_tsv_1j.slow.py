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
    # pos_samples = table['positive'].str.split(',')
    # neg_samples = table['negative'].str.split(',')
    # print('\nWorking on variant:', variant_name, 'with \n\t\tpositive samples', pos_samples, '\n\t\tnegative samples', neg_samples)

    this_positive_R1 = pos_fq_dir + variant_name + pos_R1_fq
    this_positive_R2 = pos_fq_dir + variant_name + pos_R2_fq
    this_negative_R1 = neg_fq_dir + variant_name + neg_R1_fq
    this_negative_R2 = neg_fq_dir + variant_name + neg_R2_fq
    print('\nWorking with variant', variant_name, 'files of interest:')
    print('\t', this_positive_R1, this_positive_R2, '\n\t', this_negative_R1, this_negative_R2)

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
            F = np.random.randint(1, 10000)
            if F < 10000*f:
                # print('positive')
                print('records1 = list(SeqIO.parse(this_positive_R1, "fastq"))')
                records1 = list(SeqIO.parse(this_positive_R1, "fastq"))
                print('records2 = list(SeqIO.parse(this_positive_R2, "fastq"))')
                records2 = list(SeqIO.parse(this_positive_R2, "fastq"))
                print('read_number = np.random.randint(1, len(records1))')
                read_number = np.random.randint(1, len(records1))  # min(len(records1), len(records2))
                print('\n\t', i, '/', random_DP[-1], ':', 'read_number', read_number, 'len(R1)', len(records1), 'len(R2)', len(records2))

                print('sought_id = records1[read_number].id')
                sought_id = records1[read_number].id
                print('records2 = [x for x in records2 if x.id == sought_id]')
                records2 = [x for x in records2 if x.id == sought_id]
                # print('\tnew len(R2)', len(records2))

                if records2:
                    print('random_id = str(np.random.randint(1000000, 9999999))')
                    random_id = str(np.random.randint(1000000, 9999999))
                    print("records1[read_number].id += ':' + random_id ")
                    records1[read_number].id += ':' + random_id  # '/1'
                    print("records2[0].id += ':' + random_id")
                    records2[0].id += ':' + random_id  # '/2'
                    print("sequences_R1.append(records1[read_number])")
                    sequences_R1.append(records1[read_number])
                    print("sequences_R2.append(records2[0])")
                    sequences_R2.append(records2[0])
                    print('\tF', F, 'POSITIVE R1', records1[read_number].id, 'R2', records2[0].id)
                else:
                    print('\tOops the sought id', sought_id, 'was not found in R2')
                    pass
            else:
                # print('negative')
                records1 = list(SeqIO.parse(this_negative_R1, "fastq"))
                records2 = list(SeqIO.parse(this_negative_R2, "fastq"))
                read_number = np.random.randint(1, min(len(records1), len(records2)))
                print('\n\t', i, '/', random_DP[-1], ':', 'read_number', read_number, 'len(R1)', len(records1),
                      'len(R2)', len(records2))

                sought_id = records1[read_number].id
                records2 = [x for x in records2 if x.id == sought_id]
                # print('\tnew len(R2)', len(records2))

                if records2:
                    random_id = str(np.random.randint(1000000, 9999999))
                    records1[read_number].id += ':' + random_id  # '/1'
                    records2[0].id += ':' + random_id  # '/2'
                    sequences_R1.append(records1[read_number])
                    sequences_R2.append(records2[0])
                    print('\tF', F, 'NEGATIVE R1', records1[read_number].id, 'R2', records2[0].id)
                else:
                    print('\tOops the sought id', sought_id, 'was not found in R2')
                    pass

        print('\nWriting forward (R1) and reverse (R2) FASTQ sequences')
        sequences_R1_fn = 'xls_files/' + variant_name + '.' + str(n) + '.sequences_R1.fastq'
        sequences_R2_fn = 'xls_files/' + variant_name + '.' + str(n) + '.sequences_R2.fastq'
        with open(sequences_R1_fn, "w") as output_handle:
            SeqIO.write(sequences_R1, output_handle, "fastq")
            print('\t', len(sequences_R1), 'records written to', sequences_R1_fn)
        with open(sequences_R2_fn, "w") as output_handle:
            SeqIO.write(sequences_R2, output_handle, "fastq")
            print('\t', len(sequences_R1), 'records written to', sequences_R2_fn)

        # print('\nWriting a single paired end FASTQ file')
        # output_fn = 'xls_files/' + variant_name + '.' + str(n) + '.sequences_interleaved.fastq'
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
        output_sam = "xls_files/" + variant_name + "." + str(n) + ".sequences.sam"
        output_bam = output_sam[:-4] + ".bam"
        output_bam_sorted = output_bam[:-4] + ".sorted.bam"
        output_sam = "'" + output_sam + "'"
        output_bam = "'" + output_bam + "'"
        output_bam_sorted = "'" + output_bam_sorted + "'"
        if os.path.isfile(sequences_R1_fn) and os.path.isfile(sequences_R2_fn):
            ok.append(True)
            sequences_R1_fn = "'" + sequences_R1_fn + "'"
            sequences_R2_fn = "'" + sequences_R2_fn + "'"
            print('\nAlingment:')
            print("%s mem -t 6 %s %s %s > %s" % (path['bwa'], path['ref'], sequences_R1_fn, sequences_R2_fn, output_sam))
            run("%s mem -t 6 %s %s %s > %s" % (path['bwa'], path['ref'], sequences_R1_fn, sequences_R2_fn, output_sam), shell=True)
            if os.path.isfile(output_sam[1:-1]):
                ok.append(True)
                print("\n%s view -@ 6 -b %s > %s" % (path['samtools'], output_sam, output_bam))
                run("%s view -@ 6 -b %s > %s" % (path['samtools'], output_sam, output_bam), shell=True)
                if os.path.isfile(output_bam[1:-1]):
                    ok.append(True)
                    print("\n%s sort -@ 6 %s > %s" % (path['samtools'], output_bam, output_bam_sorted))
                    run("%s sort -@ 6 %s > %s" % (path['samtools'], output_bam, output_bam_sorted), shell=True)
                    if os.path.isfile(output_bam_sorted[1:-1]):
                        ok.append(True)
                        os.remove(output_sam[1:-1])
                        os.remove(output_bam[1:-1])
                        print("\n%s index -@ 6 %s" % (path['samtools'], output_bam_sorted))
                        run("%s index -@ 6 %s" % (path['samtools'], output_bam_sorted), shell=True)
                        print('\nAlingment written to', output_bam_sorted)
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
    print('Making FASTQs and alingments for the variant', variant_name, 'took {:.3f} s'.format(time_end) +
          ' or {:.3f} m'.format(time_end / 60.0) + ' or {:.3f} h'.format(time_end / 3600.0))
    if all(ok):
        return 'SUCCESS'
    else:
        return 'FAIL'


if __name__ == '__main__':
    print('USAGE:', 'python generate_fastq_from_tsv_1j.py path/to/pos.tsv path/to/neg.tsv path/to/samples/not/needed')

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


    # xls_fn = input_dialog()
    # xls = pd.read_excel(xls_fn)
    # print('\nxls\b', xls)

    # variants = ['chr11:108121509A>G', 'chr11:108168115CTA>C', 'chr13:32912337CTG>C']
    #
    # pos_fq_dir = '/home/gkhvorykh/train_positive/Valya/'
    # neg_fq_dir = '/home/gkhvorykh/train_negative/'
    #
    # pos_R1_fq = '.positive.R1.fq'
    # pos_R2_fq = '.positive.R2.fq'
    # neg_R1_fq = '.negative.R1.fq'
    # neg_R2_fq = '.negative.R2.fq'

    # for variant in variants:
    #
        # this_positive_R1 = pos_fq_dir + variant + pos_R1_fq
        # this_positive_R2 = pos_fq_dir + variant + pos_R2_fq
        # this_negative_R1 = neg_fq_dir + variant + neg_R1_fq
        # this_negative_R2 = neg_fq_dir + variant + neg_R2_fq
        # print('\nWorking with variant', variant, 'files of interest:')
        # print('\t', this_positive_R1, this_positive_R2, '\n\t', this_negative_R1, this_negative_R2)
        #
        # f = [0.05, 0.1, 0.25, 0.5]
        # dp = [(500, 2001), (500, 2001), (250, 1501), (250, 1501)]
        #
        # f = f[0]
        # dp = dp[0]
        #
        # f = 0.25
        # dp = (300, 2000)
        # print('Chosen f', f, '\nChosen DP range', dp)
        #
        # for n in range(1, 6):  # 11
        #     # time_start = time.time()
        #     sequences_R1 = []
        #     sequences_R2 = []
        #     print('N:', n)
        #     # print('\t', [np.random.randint(x[0], x[1]) for x in dp])
        #     random_DP = range(1, np.random.randint(dp[0], dp[1]))  # range(1, np.random.randint(dp[0], dp[1])) range(1, 11)
        #     print('random_DP', random_DP[-1])
        #     for i in random_DP:
        #         F = np.random.randint(1, 10000)
        #         if F < 10000*f:
        #             # print('positive')
        #             records1 = list(SeqIO.parse(this_positive_R1, "fastq"))
        #             records2 = list(SeqIO.parse(this_positive_R2, "fastq"))
        #             read_number = np.random.randint(1, len(records1))  # min(len(records1), len(records2))
        #             print('\n\t', i, '/', random_DP[-1], ':', 'read_number', read_number, 'len(R1)', len(records1), 'len(R2)', len(records2))
        #
        #             sought_id = records1[read_number].id
        #             records2 = [x for x in records2 if x.id == sought_id]
        #             # print('\tnew len(R2)', len(records2))
        #
        #             if records2:
        #                 # records1[read_number].id += '/1'
        #                 # records2[0].id += '/2'
        #                 sequences_R1.append(records1[read_number])
        #                 sequences_R2.append(records2[0])
        #                 print('\tF', F, 'POSITIVE R1', records1[read_number].id, 'R2', records2[0].id)
        #             else:
        #                 print('\tOops the sought id', sought_id, 'was not found in R2')
        #                 pass
        #         else:
        #             # print('negative')
        #             records1 = list(SeqIO.parse(this_negative_R1, "fastq"))
        #             records2 = list(SeqIO.parse(this_negative_R2, "fastq"))
        #             read_number = np.random.randint(1, min(len(records1), len(records2)))
        #             print('\n\t', i, '/', random_DP[-1], ':', 'read_number', read_number, 'len(R1)', len(records1),
        #                   'len(R2)', len(records2))
        #
        #             sought_id = records1[read_number].id
        #             records2 = [x for x in records2 if x.id == sought_id]
        #             # print('\tnew len(R2)', len(records2))
        #
        #             if records2:
        #                 # records1[read_number].id += '/1'
        #                 # records2[0].id += '/2'
        #                 sequences_R1.append(records1[read_number])
        #                 sequences_R2.append(records2[0])
        #                 print('\tF', F, 'NEGATIVE R1', records1[read_number].id, 'R2', records2[0].id)
        #             else:
        #                 print('\tOops the sought id', sought_id, 'was not found in R2')
        #                 pass
        #
        #     print('\nWriting forward (R1) and reverse (R2) FASTQ sequences')
        #     sequences_R1_fn = 'xls_files/' + variant + '.' + str(n) + '.sequences_R1.fastq'
        #     sequences_R2_fn = 'xls_files/' + variant + '.' + str(n) + '.sequences_R2.fastq'
        #     with open(sequences_R1_fn, "w") as output_handle:
        #         SeqIO.write(sequences_R1, output_handle, "fastq")
        #         print('\t', len(sequences_R1), 'records written to', sequences_R1_fn)
        #     with open(sequences_R2_fn, "w") as output_handle:
        #         SeqIO.write(sequences_R2, output_handle, "fastq")
        #         print('\t', len(sequences_R1), 'records written to', sequences_R2_fn)
        #
        #     print('\nWriting a single paired end FASTQ file')
        #     output_fn = 'xls_files/' + variant + '.' + str(n) + '.sequences_interleaved.fastq'
        #     output_handle = open(output_fn, "w")
        #     count = 0
        #     f_iter = FastqGeneralIterator(open(sequences_R1_fn, 'r'))
        #     r_iter = FastqGeneralIterator(open(sequences_R2_fn, 'r'))
        #     for (f_id, f_seq, f_q), (r_id, r_seq, r_q) in zip(f_iter, r_iter):
        #         assert f_id == r_id
        #         count += 2
        #         # Write out both reads with "/1" and "/2" suffix on ID
        #         output_handle.write("@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n"
        #                             % (f_id, f_seq, f_q, r_id, r_seq, r_q))
        #     output_handle.close()
        #     print("\t%s records written to %s" % (count, output_fn))
        #
        #     # Running bwa and samtools to align, sort, index the sequences
        #     output_sam = "xls_files/" + variant + "." + str(n) + ".sequences.sam"
        #     output_bam = output_sam[:-4] + ".bam"
        #     output_bam_sorted = output_bam[:-4] + ".sorted.bam"
        #     output_sam = "'" + output_sam + "'"
        #     output_bam = "'" + output_bam + "'"
        #     output_bam_sorted = "'" + output_bam_sorted + "'"
        #     if os.path.isfile(sequences_R1_fn) and os.path.isfile(sequences_R2_fn):
        #         sequences_R1_fn = "'" + sequences_R1_fn + "'"
        #         sequences_R2_fn = "'" + sequences_R2_fn + "'"
        #         print('\nAlingment:')
        #         print("%s mem -t 6 %s %s %s > %s" % (path['bwa'], path['ref'], sequences_R1_fn, sequences_R2_fn, output_sam))
        #         run("%s mem -t 6 %s %s %s > %s" % (path['bwa'], path['ref'], sequences_R1_fn, sequences_R2_fn, output_sam), shell=True)
        #         if os.path.isfile(output_sam[1:-1]):
        #             print("\n%s view -@ 6 -b %s > %s" % (path['samtools'], output_sam, output_bam))
        #             run("%s view -@ 6 -b %s > %s" % (path['samtools'], output_sam, output_bam), shell=True)
        #             if os.path.isfile(output_bam[1:-1]):
        #                 print("\n%s sort -@ 6 %s > %s" % (path['samtools'], output_bam, output_bam_sorted))
        #                 run("%s sort -@ 6 %s > %s" % (path['samtools'], output_bam, output_bam_sorted), shell=True)
        #                 if os.path.isfile(output_bam_sorted[1:-1]):
        #                     os.remove(output_sam[1:-1])
        #                     os.remove(output_bam[1:-1])
        #                     print("\n%s index -@ 6 %s" % (path['samtools'], output_bam_sorted))
        #                     run("%s index -@ 6 %s" % (path['samtools'], output_bam_sorted), shell=True)
        #                     print('\nAlingment written to', output_bam_sorted)
        #                 else:
        #                     print('ERROR:', output_bam_sorted, 'does not exist')
        #             else:
        #                 print('ERROR:', output_bam, 'does not exist')
        #         else:
        #             print('ERROR:', output_sam, 'does not exist')
