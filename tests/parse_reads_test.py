#make sure we can import from package directory
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) 

#use this to load parse_reads_cy properly
import pyximport; pyximport.install()

import unittest
import string, math
import amplimap.common
import amplimap.reader
from amplimap.parse_reads import *

class ParseReads(unittest.TestCase):
    def make_qual(self, qs, base = 33):
        return ''.join([ chr(base + max(0, min(40, q))) for q in qs ])
    def make_qual_from_p(self, ps, base = 33):
        return self.make_qual([ math.floor(-10 * math.log10(p)) for p in ps ], base)

    def test_mismatch_count_with_max(self):
        #perfect match
        self.assertEqual(mismatch_count_with_max(a=b"X"*10, b=b"X"*10, max_mismatches=0), 0)
        self.assertEqual(mismatch_count_with_max(a=b"X"*10, b=b"X"*10, max_mismatches=10), 0)
        #one mismatch
        self.assertEqual(mismatch_count_with_max(a=b"X"*10, b=b"X"*9 + b"B"*1, max_mismatches=0), 1)
        self.assertEqual(mismatch_count_with_max(a=b"X"*10, b=b"X"*9 + b"B"*1, max_mismatches=1), 1)
        self.assertEqual(mismatch_count_with_max(a=b"X"*10, b=b"X"*9 + b"B"*1, max_mismatches=10), 1)
        #two mismatches
        self.assertEqual(mismatch_count_with_max(a=b"A"+b"X"*9, b=b"X"*9 + b"B"*1, max_mismatches=0), 1) #we stop after 1
        self.assertEqual(mismatch_count_with_max(a=b"A"+b"X"*9, b=b"X"*9 + b"B"*1, max_mismatches=1), 2)
        self.assertEqual(mismatch_count_with_max(a=b"A"+b"X"*9, b=b"X"*9 + b"B"*1, max_mismatches=2), 2)
        self.assertEqual(mismatch_count_with_max(a=b"A"+b"X"*9, b=b"X"*9 + b"B"*1, max_mismatches=10), 2)

        #dots
        self.assertEqual(mismatch_count_with_max(a=b"X" + b"."*9, b=b"X"*9 + b"B"*1, max_mismatches=10), 0)
        self.assertEqual(mismatch_count_with_max(a=b"A" + b"."*9, b=b"X"*9 + b"B"*1, max_mismatches=10), 1)
        self.assertEqual(mismatch_count_with_max(a=b"."*10, b=b"X"*9 + b"B"*1, max_mismatches=10), 0)
        self.assertEqual(mismatch_count_with_max(a=b"A"*10, b=b"X"*9 + b".", max_mismatches=10), 9)

        #different length
        self.assertEqual(mismatch_count_with_max(a=b"A"*100, b=b"A"*10, max_mismatches=0), 0)
        self.assertEqual(mismatch_count_with_max(a=b"A"*100, b=b"X"*9 + b"B"*1, max_mismatches=0), 1)
        self.assertEqual(mismatch_count_with_max(a=b"A"*10, b=b"A"*100, max_mismatches=0), 90)
        self.assertEqual(mismatch_count_with_max(a=b"A"*10, b=b"A"*100, max_mismatches=90), 90)
        self.assertEqual(mismatch_count_with_max(a=b"A"*10, b=b"A"*100, max_mismatches=100), 90)

        #special cases
        self.assertEqual(mismatch_count_with_max(a=b"", b=b"", max_mismatches=10), 0)
        self.assertEqual(mismatch_count_with_max(a=b"", b=b"ABC", max_mismatches=10), 3)

        #errors
        with self.assertRaises(TypeError):
            mismatch_count_with_max(a=b"ABC", b=None, max_mismatches=10)
        with self.assertRaises(TypeError):
            mismatch_count_with_max(a=None, b=b"ABC", max_mismatches=10)

    def test_quality_trim_read(self):
        seq = (string.ascii_uppercase * 4)[0:100]

        qual_all_good = self.make_qual_from_p([0.001] * 100)
        qual_all_bad = self.make_qual_from_p([0.5] * 100)

        qual_all_good_base64 = self.make_qual_from_p([0.001] * 100, base = 64)

        #all or nothing
        self.assertEqual(quality_trim_read(seq, qual_all_good), (100, seq, qual_all_good))
        self.assertEqual(quality_trim_read(seq, qual_all_bad), (0, '', ''))

        #different phred bases
        self.assertEqual(quality_trim_read(seq, qual_all_good_base64, phred_base = 64), (100, seq, qual_all_good_base64))
        with self.assertRaises(AssertionError):
            quality_trim_read(seq, qual_all_good_base64)
        with self.assertRaises(AssertionError):
            quality_trim_read(seq, qual_all_good, phred_base = 64)

        #example from http://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Quality_trimming.html
        #see quality_trim_test_spreadsheet.xls for these numbers
        qual_example = self.make_qual([4, 8, 15, 30, 32, 23, 10, 31, 31, 20, 28, 29, 30, 10, 31, 31, 20, 28, 32, 31, 32, 33, 29, 30, 31, 32, 32, 31, 32, 33, 33, 32, 33, 34, 30, 31, 32, 33, 33, 32, 33, 34, 23, 7, 7, 12])
        seq_example = seq[:len(qual_example)]
        for limit_example, exp_start, exp_end in [
            (0.0005, 33, 41),
            (0.001, 4, 41),
            (0.05, 2, 42),
            (0.2, 1, 45),
            (0.4, 0, 45),
        ]:
            self.assertEqual(
                quality_trim_read(seq_example, qual_example, p_threshold = limit_example),
                (exp_end - exp_start + 1, seq_example[exp_start:exp_end+1], qual_example[exp_start:exp_end+1]))

    def test_make_extended_read_name(self):
        for read_id in ['', 'asdfdasf::', '::1::1::1']:
            for probe in ['', 'probe-1', '_probe']:
                for umi in ['', 'ACG', 'ACGTTTTTTT']:
                    if len(read_id) == 0 or len(probe) == 0:
                        with self.assertRaises(Exception):
                            amplimap.common.parse_extended_read_name(amplimap.common.make_extended_read_name(read_id, probe, umi))
                    else:
                        self.assertEqual(
                            (read_id, probe, umi),
                            amplimap.common.parse_extended_read_name(amplimap.common.make_extended_read_name(read_id, probe, umi))
                        )

    def test_make_trimmed_read(self):
        for do_find_second_arm in [False]: #should add this -- find_arm_sequence and find_arm_sequence_mismatches
            for trim_min_length in [None, 20, 200]:
                for trim_primers in [True, False]:
                    for umi in ["99999", "", "1111111111111"]:
                        target_len = 100
                        arm_len = 10
                        other_arm_len = 20
                        umi_len = len(umi)

                        read_seq_trimmed = (string.ascii_uppercase * 4)[0:target_len]
                        read_qual_trimmed = self.make_qual_from_p([0.001] * target_len)

                        read_seq_maybe_primers = read_seq_trimmed if trim_primers else "X"*arm_len + read_seq_trimmed + "Y"*other_arm_len
                        read_qual_maybe_primers = read_qual_trimmed if trim_primers else self.make_qual_from_p([0.1] * (arm_len)) + read_qual_trimmed + self.make_qual_from_p([0.2] * (other_arm_len))

                        read_seq = umi + "X"*arm_len + read_seq_trimmed + "Y"*other_arm_len + umi
                        read_qual = self.make_qual_from_p([0.1] * (umi_len+arm_len)) + read_qual_trimmed + self.make_qual_from_p([0.2] * (other_arm_len+umi_len))
                        read_id = 'READ:ID'
                        read = ["%s adsfdasfsadfsadfsdaf" % read_id, read_seq, read_qual]

                        probe = "a-zA-Z0-9._-"

                        primer_trim_end = umi_len+arm_len+target_len+(0 if trim_primers else other_arm_len)
                        quality_trimmed_difference = 0

                        infos = 'pr:Z:%s\tum:Z:%s\tol:i:%d\tts:i:%d\tte:i:%d' % (probe, umi if umi_len > 0 else 'X', len(read_seq), umi_len+(arm_len if trim_primers else 0), primer_trim_end)
                        fastq_lines = (
                            '@%s\t%s\n' % (amplimap.common.make_extended_read_name(read_id, probe, umi), infos) +
                            '%s\n' % read_seq_maybe_primers +
                            '+\n' +
                            '%s\n' % read_qual_maybe_primers
                        )

                        if trim_min_length == 200:
                            expected = (None, primer_trim_end, False, 0, len(read_seq_maybe_primers))
                        else:
                            expected = (fastq_lines, primer_trim_end, True, quality_trimmed_difference, len(read_seq_maybe_primers))

                        #print(do_find_second_arm, trim_min_length, trim_primers, umi, '-> expected:')
                        #print(expected)

                        output = make_trimmed_read(
                            read_id,
                            read,
                            probe,
                            target_len,
                            umi,
                            umi_len,
                            arm_len,
                            other_arm_len,
                            trim_primers = trim_primers,
                            trim_min_length = trim_min_length
                        )

                        #print(do_find_second_arm, trim_min_length, trim_primers, umi, '-> actual:')
                        #print(output)

                        self.assertEqual(output, expected)

    def test_find_probe_from_reads(self):
        #arm sequences to simulate and search for
        probes_dict = {
            'first_primer_5to3': {
                'p1': 'A'*10, 'p2': 'A'*10, 'p3': 'B'*15, 'p4': 'C'*20, 'p5': 'Y'*10, 'p6': 'Z'*20
            },
            #in read orientation (will be reversed once to make fragment seq and back again to make read seq)
            'second_primer_5to3': {
                'p1': 'm'*15, 'p2': 'n'*5, 'p3': 'o'*15, 'p4': 'p'*1, 'p5': 'q'*15, 'p6': 'r'*10
            }
        }

        #convert input primers into expected format
        first_arms_probes, first_arms_seq = zip(*probes_dict['first_primer_5to3'].items()) #split into lists

        #convert to utf-8
        first_arms_seq = tuple([seq.encode('utf-8') for seq in first_arms_seq])
        first_arms_dict = dict([(key, val.encode('utf-8')) for key, val in probes_dict['first_primer_5to3'].items()])
        second_arms_dict = dict([(key, val.encode('utf-8')) for key, val in probes_dict['second_primer_5to3'].items()])

        #length of the target region (between arms)
        for target_len in [100, 10]:
            #first arm
            for probe_first_arm in first_arms_dict.keys():
                #second arm
                for probe_second_arm in second_arms_dict.keys():
                    #length of first umi
                    for umi_one in [0, 5]:
                        #length of second umi
                        for umi_two in [0, 5]:
                            #number of mismatches to allow
                            for mismatches in [0, 2, 5]:
                                #where to add mismatches
                                for add_mismatch_to in ['first', 'both']:
                                    #how many mismatches to add
                                    for add_mismatches in [0, 1, 2, 3]:
                                        #read length to simulate
                                        for read_len in [10, 11, 75, 200]:
                                            #arms
                                            first_arm = first_arms_dict[probe_first_arm]
                                            second_arm = second_arms_dict[probe_second_arm]

                                            #target sequence
                                            fragment_seq_trimmed = (string.ascii_uppercase * 4)[0:target_len].encode('utf8')

                                            #add mismatches to arms
                                            for mismatch_ix in range(add_mismatches):
                                                if mismatch_ix < len(first_arm):
                                                    first_arm = first_arm[:mismatch_ix] + b'#' + first_arm[(mismatch_ix+1):]
                                                if add_mismatch_to == 'both' and mismatch_ix < len(second_arm):
                                                    second_arm = second_arm[:mismatch_ix] + b'#' + second_arm[(mismatch_ix+1):]

                                            #build whole region seq
                                            fragment_seq_with_primers = first_arm + fragment_seq_trimmed + second_arm[::-1] #note reverse of second arm
                                            fragment_seq_with_primers_umis = b'1' * umi_one + fragment_seq_with_primers + b'2' * umi_two

                                            #simulate reads from this sequence
                                            first_read_str = fragment_seq_with_primers_umis[:read_len:1]
                                            second_read_str = fragment_seq_with_primers_umis[:-(read_len+1):-1] #:-1 reverses, but need to go to -(len+1) then

                                            expected = []
                                            #if first arm is 0 then second arm len can actually be just
                                            if read_len >= max(len(first_arm)+umi_one, len(second_arm)+umi_two):
                                                first_second_hybrid = (probe_first_arm == 'p1' and probe_second_arm == 'p2')
                                                second_first_hybrid = (probe_first_arm == 'p2' and probe_second_arm == 'p1')
                                                if (probe_first_arm == probe_second_arm) or first_second_hybrid or second_first_hybrid:
                                                    if add_mismatches <= mismatches:
                                                        expected.append([probe_first_arm if not first_second_hybrid and not second_first_hybrid else probe_second_arm,
                                                            add_mismatches, min(add_mismatches, len(second_arm)) if add_mismatch_to == 'both' else 0])

                                            #in this case we can still find the second arm pair
                                            #second probes matches first with five mismatches, since first arm is same and second arm just has 5 bp
                                            if probe_first_arm in ['p1', 'p2'] and probe_second_arm != 'p2':
                                                if read_len >= max(len(first_arm)+umi_one, len(second_arms_dict['p2'])+umi_two):
                                                    if mismatches >= 5 and add_mismatches <= mismatches:
                                                        expected.append(['p2', min(add_mismatches, len(first_arm)), 5])

                                            if probe_first_arm in ['p4'] and probe_second_arm != 'p4':
                                                if read_len >= max(len(first_arm)+umi_one, len(second_arms_dict['p4'])+umi_two):
                                                    if mismatches >= 1 and add_mismatches <= mismatches:
                                                        expected.append(['p4', min(add_mismatches, len(first_arm)), 1])


                                            output = find_probe_from_reads(
                                                first_arms_probes,
                                                first_arms_seq,
                                                second_arms_dict,
                                                first_read_str,
                                                second_read_str,
                                                umi_one = umi_one,
                                                umi_two = umi_two,
                                                mismatches = mismatches
                                            )

                                            try:
                                                self.assertEqual(
                                                    sorted(output, key=operator.itemgetter(0)),
                                                    sorted(expected, key=operator.itemgetter(0))
                                                )
                                            except AssertionError:
                                                #rerun with debug
                                                find_probe_from_reads(
                                                    first_arms_probes,
                                                    first_arms_seq,
                                                    second_arms_dict,
                                                    first_read_str,
                                                    second_read_str,
                                                    umi_one = umi_one,
                                                    umi_two = umi_two,
                                                    mismatches = mismatches,
                                                    debug = True
                                                )

                                                print('RL=', read_len, 'arms=', probe_first_arm, probe_second_arm, 'max mm =', mismatches, 'add mm =', add_mismatches, 'to', add_mismatch_to, '->')
                                                print(first_arm)
                                                print(second_arm)
                                                print('len normal=', max(len(first_arm), len(second_arm)) + max(umi_one, umi_two))
                                                print('len hybrid=', max(len(first_arm), len(second_arms_dict['p2'])) + max(umi_one, umi_two))
                                                print('total=', fragment_seq_with_primers_umis)
                                                print('1st=', first_read_str)
                                                print('2nd=', second_read_str)
                                                print('exp:')
                                                print(expected)
                                                print('obs:')
                                                print(output)
                                                raise


if __name__ == '__main__':
    unittest.main()