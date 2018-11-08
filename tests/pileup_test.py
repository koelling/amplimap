#make sure we can import from package directory
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) 

import unittest
import string, math
import pysam
import amplimap.common
import amplimap.reader
from amplimap.pileup import *

class SimulatedPileupRead():
    pass

class Pileups(unittest.TestCase):
    def test_process_pileup_read(self):
        for mapq in [0, 20, 100]:
            for vtype in ['snp', 'ins', 'del', 'skip']:
                for naive_umi in ['', 'NNNN', 'ACGTACGT']:
                    for add_softclip in [True, False]:
                        for test_validate  in [0, 1, 2, 3, 4, 5]:
                            should_fail = False

                            pr = None #see below
                            is_del = vtype == 'del'
                            is_refskip = vtype == 'skip'
                            indel = -1 if vtype == 'del' else 1 if vtype == 'ins' else 0

                            probes_dict = None
                            reference_name = 'chrZ'
                            reference_pos = 999
                            min_mapq = 20
                            ignore_groups = False
                            use_naive_groups = False
                            validate_probe_targets = False

                            query_position = 0 if vtype == 'snp' else None

                            read_sequence = 'ACGTA' * 10
                            read_qualities = [40] * len(read_sequence)
                            cigarstring = '%dM' % len(read_sequence) if not add_softclip else '5S%dM' % (len(read_sequence) - 5)
                            umi_group = 'group1'
                            probe_id = 'probe1'
                            original_qname = 'read1'
                            qname = amplimap.common.make_extended_read_name(original_qname, probe_id, naive_umi)

                            #test probe validation with different possibilities
                            if test_validate > 0:
                                validate_probe_targets = True
                                if test_validate == 1:
                                    probes_dict = { 'chr': { }, 'target_start0': {}, 'target_end': {} }
                                    should_fail = True
                                elif test_validate == 1:
                                    probes_dict = { 'chr': { probe_id: reference_name }, 'target_start0': { probe_id: reference_pos }, 'target_end': { probe_id: reference_pos+10 } }
                                    should_fail = False
                                elif test_validate == 2:
                                    probes_dict = { 'chr': { probe_id: reference_name }, 'target_start0': { probe_id: reference_pos - 10 }, 'target_end': { probe_id: reference_pos } }
                                    should_fail = False
                                elif test_validate == 3:
                                    probes_dict = { 'chr': { probe_id: reference_name }, 'target_start0': { probe_id: reference_pos - 10 }, 'target_end': { probe_id: reference_pos - 1 } }
                                    should_fail = True
                                elif test_validate == 4:
                                    probes_dict = { 'chr': { probe_id: 'FAIL' }, 'target_start0': { probe_id: reference_pos }, 'target_end': { probe_id: reference_pos+10 } }
                                    should_fail = True
                                elif test_validate == 5:
                                    probes_dict = { 'chr': { probe_id: reference_name, 'OTHER': reference_name },
                                        'target_start0': { probe_id: reference_pos - 10, 'OTHER': reference_pos-10 },
                                        'target_end': { probe_id: reference_pos - 1, 'OTHER': reference_pos+10 } }
                                    should_fail = True

                            pr = SimulatedPileupRead()
                            pr.is_del = is_del
                            pr.is_refskip = is_refskip
                            pr.indel = indel
                            pr.query_position = query_position #position of the read base at the pileup site, 0-based.

                            pr.alignment = pysam.AlignedSegment()

                            #should be returned
                            pr.alignment.qname = qname

                            #will compare to min_mapq
                            pr.alignment.mapping_quality = mapq

                            #will be checking this for 'S'
                            pr.alignment.cigarstring = cigarstring

                            #read sequence bases, including soft clipped bases (None if not present).
                            pr.alignment.query_sequence = read_sequence
                            #Quality scores are returned as a python array of unsigned chars. Note that this is not the ASCII-encoded value typically seen in FASTQ or SAM formatted files. Thus, no offset of 33 needs to be subtracted.
                            pr.alignment.query_qualities = read_qualities

                            tags = {'pr': probe_id, 'UG': umi_group, 'um': naive_umi}
                            pr.alignment.set_tags(tags.items())

                            if mapq < min_mapq or 'N' in naive_umi or add_softclip or should_fail:
                                with self.assertRaises(PileupRowFilterException): 
                                    process_pileup_read(
                                        pr,
                                        probes_dict,
                                        reference_name,
                                        reference_pos,
                                        min_mapq,
                                        ignore_groups,
                                        use_naive_groups,
                                        validate_probe_targets,
                                        debug = False)
                            else:
                                my_call = None
                                my_phred = None
                                if vtype == 'del':
                                    my_call = 'DEL'
                                    my_phred = 40
                                elif vtype == 'skip':
                                    my_call = 'SKIP'
                                    my_phred = 40
                                elif vtype == 'ins':
                                    my_call = 'INS'
                                    my_phred = 40
                                elif vtype == 'snp':
                                    my_call = read_sequence[query_position]
                                    my_phred = read_qualities[query_position]                            

                                expect = (my_call, my_phred, original_qname, umi_group, probe_id, naive_umi)

                                output = process_pileup_read(
                                    pr,
                                    probes_dict,
                                    reference_name,
                                    reference_pos,
                                    min_mapq,
                                    ignore_groups,
                                    use_naive_groups,
                                    validate_probe_targets,
                                    debug = False)

                                self.assertEqual(expect, output)

    def test_process_pileup_read(self):
        pass

    def test_record_read_in_group(self):
        read_calls = {}
        read_calls_expect = {}

        my_call = 'A'
        my_phred = 10
        my_umi = 'ACGTA'
        qname = 'read1'
        group = 'group1'

        #should be added
        read_calls_expect[qname] = [my_call, my_phred, my_umi]
        output = record_read_in_group(read_calls, my_call, my_phred, my_umi, qname)
        self.assertEqual(output, False)
        self.assertEqual(read_calls, read_calls_expect)
        self.assertEqual(len(read_calls), 1)

        #should be added as well
        qname = 'read2'
        read_calls_expect[qname] = [my_call, my_phred, my_umi]
        output = record_read_in_group(read_calls, my_call, my_phred, my_umi, qname)
        self.assertEqual(output, False)
        self.assertEqual(read_calls, read_calls_expect)
        self.assertEqual(len(read_calls), 2)

        #should be ignored
        my_phred = 5
        output = record_read_in_group(read_calls, my_call, my_phred, my_umi, qname)
        self.assertEqual(output, True)
        self.assertEqual(read_calls, read_calls_expect)

        #phred should be changed
        my_phred = 15
        read_calls_expect[qname] = [my_call, my_phred, my_umi]
        output = record_read_in_group(read_calls, my_call, my_phred, my_umi, qname)
        self.assertEqual(output, True)
        self.assertEqual(read_calls, read_calls_expect)

        #call should be changed
        my_call = 'DEL'
        my_phred = 25
        read_calls_expect[qname] = [my_call, my_phred, my_umi]
        output = record_read_in_group(read_calls, my_call, my_phred, my_umi, qname)
        self.assertEqual(output, True)
        self.assertEqual(read_calls, read_calls_expect)

        #call should be N
        my_call = 'INS'
        my_phred = 25
        read_calls_expect[qname] = ['N', my_phred, my_umi]
        output = record_read_in_group(read_calls, my_call, my_phred, my_umi, qname)
        self.assertEqual(output, True)
        self.assertEqual(read_calls, read_calls_expect)

        #call should be set back to the call now, since we got a higher phred
        my_call = 'G'
        my_phred = 30
        read_calls_expect[qname] = [my_call, my_phred, my_umi]
        output = record_read_in_group(read_calls, my_call, my_phred, my_umi, qname)
        self.assertEqual(output, True)
        self.assertEqual(read_calls, read_calls_expect)

        #should raise
        qname = 'read2'
        my_umi = 'NNNN'
        with self.assertRaises(AssertionError):
            record_read_in_group(read_calls, my_call, my_phred, my_umi, qname)

        #check we didn't mess anything up
        self.assertEqual(len(read_calls), 2)

    def test_get_group_consensus(self):
        self.assertEqual( ('DEL', 1, 30), get_group_consensus( [('DEL', 30, None)], min_consensus_count = 1, ignore_groups = False ))
        self.assertEqual( ('SKIP', 1, 30), get_group_consensus( [('SKIP', 30, None)], min_consensus_count = 1, ignore_groups = False ))
        self.assertEqual( ('A', 1, 30), get_group_consensus( [('A', 30, None)], min_consensus_count = 1, ignore_groups = False ))

        self.assertEqual( ('A', 10, 30), get_group_consensus( [('A', 30, None)]*10, min_consensus_count = 2, ignore_groups = False ) )
        self.assertEqual( ('A', 10, 30), get_group_consensus( [('A', 30, None)]*10, min_consensus_count = 10, ignore_groups = False ) )

        self.assertEqual( ('A', 12, 31), get_group_consensus( [('A', 30, None)]*10 + [('G', 40, None), ('A', 31, None), ('A', 20, None)], min_consensus_count = 2, ignore_groups = False ))
        self.assertEqual( ('G', 11, 40), get_group_consensus( [('G', 20, None), ('A', 30, None)]*10 + [('G', 40, None)], min_consensus_count = 2, ignore_groups = False ))
        self.assertEqual( ('G', 11, 40), get_group_consensus( [('G', 20, None), ('A', 30, None)]*10 + [('G', 40, None)], min_consensus_count = 2, ignore_groups = False ))

        #no base call
        with self.assertRaises(PileupGroupFilterException):
            get_group_consensus( [('N', 30, None)], min_consensus_count = 1, ignore_groups = False )
        with self.assertRaises(PileupGroupFilterException):
            get_group_consensus( [('N', 30, None)], min_consensus_count = 1, ignore_groups = True )

        #no consensus (10 and 10)
        with self.assertRaises(PileupGroupFilterException): #10+10
            get_group_consensus( [('G', 20, None), ('A', 30, None)]*10, min_consensus_count = 2, ignore_groups = False )
        with self.assertRaises(PileupGroupFilterException): #15+10+10
            get_group_consensus( [('G', 20, None)]*15 + [('A', 30, None)]*10 + [('C', 30, None)]*10, min_consensus_count = 2, ignore_groups = False )
        #but fine if we ignore groups
        #but note: we just take the first value here, which happens to be G. not ideal, but should be a rare case
        #however, this should never happen in practise because if ignore_groups is true then we should only get one call per group (group=read)
        self.assertEqual( ('G', 10, 20), get_group_consensus( [('G', 20, None), ('A', 30, None)]*10, min_consensus_count = 2, ignore_groups = True ) )

        #min consensus count
        with self.assertRaises(PileupGroupFilterException):
            get_group_consensus( [('A', 30, None)], min_consensus_count = 2, ignore_groups = False )
        #but fine if we ignore groups
        self.assertEqual( ('A', 1, 30), get_group_consensus( [('A', 30, None)], min_consensus_count = 100, ignore_groups = True ))

        #min consensus fraction
        #check default - min 51%
        self.assertEqual( ('A', 2, 30), get_group_consensus( [('A', 25, None), ('G', 20, None), ('A', 30, None)] )) #1+2
        self.assertEqual( ('A', 3, 35), get_group_consensus( [('A', 25, None), ('G', 20, None), ('A', 30, None), ('A', 35, None)] )) #1+3
        with self.assertRaises(PileupGroupFilterException):
            get_group_consensus( [('G', 20, None), ('A', 30, None)] ) #1+1
        with self.assertRaises(PileupGroupFilterException):
            get_group_consensus( [('G', 20, None), ('A', 30, None)]*2 ) #2+2
        self.assertEqual( ('G', 11, 20), get_group_consensus( [('G', 20, None), ('A', 30, None)]*10 + [('G', 20, None)] )) #11+10
        with self.assertRaises(PileupGroupFilterException):
            get_group_consensus( [('G', 20, None), ('A', 30, None)]*10 ) #10+10
        #check with 75%
        self.assertEqual( ('A', 30, 31), get_group_consensus( [('G', 20, None), ('A', 30, None)]*10 + [('A', 31, None)]*20, min_consensus_fraction = 0.75 )) #30+10
        with self.assertRaises(PileupGroupFilterException):
            get_group_consensus( [('G', 20, None), ('A', 30, None)]*10 + [('A', 30, None)]*20, min_consensus_fraction = 0.76 ) #30+10
        #check with 0% (just get the first one for ties)
        self.assertEqual( ('G', 1, 20), get_group_consensus( [('G', 20, None), ('A', 30, None)], min_consensus_fraction = 0.0 )) #1+1
        #...but still the correct one for two
        self.assertEqual( ('A', 2, 30), get_group_consensus( [('A', 25, None), ('G', 20, None), ('A', 30, None)], min_consensus_fraction = 0.0 )) #1+2
        #check with 100%
        self.assertEqual( ('G', 1, 20), get_group_consensus( [('G', 20, None)], min_consensus_fraction = 1.0 )) #1+1

    def umi_groups_check(self, output, groups):
        seen = []
        for group in groups:
            my = None
            for item in group:
                if my is None:
                    my = output[item]
                    assert my not in seen
                    seen.append(my)

                assert output[item] == my, '%s %s %s' % (str(my), str(item), str(output[item]))
        return True

    def test_find_umi_groups(self):
        #easy - just two different UMIs
        self.umi_groups_check( amplimap.common.find_umi_groups({ b'AAAAA': 10, b'GGGGG': 10 }),
            [[b'AAAAA'], [b'GGGGG']])
        #one mismatch, different sizes
        self.umi_groups_check( amplimap.common.find_umi_groups({ b'AAAAA': 10, b'AAAAG': 1, b'GGGGG': 10 }),
            [[b'AAAAA', b'AAAAG'], [b'GGGGG']])
        #one N, different sizes
        self.umi_groups_check( amplimap.common.find_umi_groups({ b'AAAAA': 10, b'NAAAA': 1, b'GGGGG': 10 }),
            [[b'AAAAA', b'NAAAA'], [b'GGGGG']])
        #one mismatch but half size
        self.umi_groups_check( amplimap.common.find_umi_groups({ b'AAAAA': 10, b'AAAAG': 5, b'GGGGG': 10 }),
            [[b'AAAAA', b'AAAAG'], [b'GGGGG']])
        #one mismatch but equal size
        self.umi_groups_check( amplimap.common.find_umi_groups({ b'AAAAA': 10, b'AAAAG': 10, b'GGGGG': 10 }),
            [[b'AAAAA'], [b'AAAAG'], [b'GGGGG']])
        #one N equal size
        self.umi_groups_check( amplimap.common.find_umi_groups({ b'AAAAA': 10, b'NAAAA': 10, b'GGGGG': 10 }),
            [[b'AAAAA'], [b'NAAAA'], [b'GGGGG']]) #note: N is consiered like any other nucleotide, probably correct
        #two mismatches, but one edit to next
        self.umi_groups_check( amplimap.common.find_umi_groups({ b'AAAAA': 10, b'AAAAG': 1, b'TAAAG': 1, b'GGGGG': 10 }),
            [[b'AAAAA', b'AAAAG', b'TAAAG'], [b'GGGGG']])
        #two mismatches, but without the step in-between
        self.umi_groups_check( amplimap.common.find_umi_groups({ b'AAAAA': 10, b'TAAAG': 1, b'GGGGG': 10 }),
            [[b'AAAAA'], [b'TAAAG'], [b'GGGGG']])
        
if __name__ == '__main__':
    unittest.main()