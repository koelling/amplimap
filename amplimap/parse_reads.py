# -*- coding: utf-8 -*-
"""
This module contains methods for reading raw FASTQ files, finding primers in them and trimming the reads.
"""

#python 3 compat
#http://python-future.org/compatible_idioms.html
from __future__ import print_function

import sys
import os

import gzip

import numpy as np
import re
import Bio.SeqIO
#import Bio.Seq
#import Bio.pairwise2
#import itertools #for efficient zipping

import glob
import shutil #for copy

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
log.addHandler(sh)

import time

#use pandas to read the CSV file and write output files
import pandas as pd
#to compute string distances
#import distance
#for defaultdict + sorting
import collections
import operator
#for debug
import random

#use python's argumentparser to support command line parameters like --probe=62
import argparse

#functions for reading input files
from .common import make_extended_read_name, find_umi_groups, parse_extended_read_name
from .reader import read_new_probe_design

#import cython module. if we have run setup.py this should have been built already
#but we may get import errors (eg. for readthedocs)
try:
    from .parse_reads_cy import find_probe_from_reads, mismatch_count_with_max
except ImportError:
    sys.stderr.write('Failed to import parse_reads_cy, trying again with pyximport\n')
    
    #try again with pyximport
    import pyximport
    pyximport.install()
    from .parse_reads_cy import find_probe_from_reads, mismatch_count_with_max


def quality_trim_read(read_seq, read_qual_raw, phred_base: int = 33, p_threshold: float = 0.01) -> tuple:
    """Use Richard Mott's algorithm implemented in bwa and others to quality-trim reads.

    Returns:
        tuple: length of read after quality trimming, trimmed sequence, trimmed quality scores

    - see description: http://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Quality_trimming.html
    - see sample code (but not sure if correct!): https://www.biostars.org/p/1923/
    - see version in bio.seqio: http://biopython.org/DIST/docs/api/Bio.SeqIO.AbiIO-module.html
    - see implementation (may be inefficient, keeps list of all cumulative scores...): http://biopython.org/DIST/docs/api/Bio.SeqIO.AbiIO-pysrc.html#_abi_trim

    NOTE: our implementation will always trim from first high-q base (sum > 0) to the best part of the read
    (max of the running sum), even if there are some low-quality bases in-between."""
    assert len(read_seq) == len(read_qual_raw)

    running_sum = 0
    running_sum_max = 0
    quality_trim_start = None
    quality_trim_end = None

    for ix in range(len(read_qual_raw)):
        #calculate phred score from character
        phred = ord(read_qual_raw[ix]) - phred_base
        assert phred >= 0 and phred <= 40, ('Invalid phred quality encountered: %s -- %s -> %d. Incorrect phred base?' % (read_qual_raw, read_qual_raw[ix], phred))

        #calculate error probability from phred
        p_error = 10 ** (-phred / 10)

        #calculate offset from limit (will be negative if error prob > threshold)
        offset = p_threshold - p_error

        #add offset to running sum
        running_sum += offset

        if running_sum > 0:
            #set quality_trim_start to current index the first time we get positive sum
            #NOTE: this condition means that we will always cut from the first good base until the last max base.
            #it's possible that we don't want this, since it does mean we include low-quality bases in-between
            #not sure what's better, the description of the algorithm is a bit ambivalent there.
            #unlikely to matter much anyway...
            if quality_trim_start is None:
                quality_trim_start = ix

            #record max position of the running sum
            if running_sum >= running_sum_max:
                running_sum_max = running_sum
                quality_trim_end = ix + 1
        #if sum < 0, set to 0
        elif running_sum < 0:
            running_sum = 0

    #did we find a start and an end?
    if quality_trim_start is not None and quality_trim_end is not None:
        assert quality_trim_end > quality_trim_start

        #note quality_trim_end should be index of last character + 1, so end-start = length
        return quality_trim_end - quality_trim_start, read_seq[quality_trim_start:quality_trim_end], read_qual_raw[quality_trim_start:quality_trim_end]
    else:
        #return empty result, this should be filtered out downstream
        return 0, '', ''

def make_trimmed_read(
        read_id,
        read, probe,
        target_len,
        umi, umi_len,
        arm_len, other_arm_len,
        trim_primers = False, trim_min_length = None,
        check_sequence = None, check_sequence_mismatches = None,
        find_arm_sequence = None, find_arm_sequence_mismatches = None,
        quality_trim_threshold = None, quality_trim_phred_base = None):
    """
    Take a fastq read tuple and generate new fastq string after removing umi and potential opposite primer arm.

    Returns:
        tuple:
            - str: FastQ lines
            - int: Index of first base trimmed off from the end (to trim opposite primer)
            - bool: was trimmed len great than trim_min_len
            - int: change after quality trimming
            - int: trimmed length
    """

    if check_sequence is not None:
        raise Exception('not implemented fully')
        other_arm_seq = read[1][umi_len+arm_len+target_len:umi_len+arm_len+target_len+other_arm_len]
        print(other_arm_seq)
        print(check_sequence)
        dist = mismatch_count_with_max(other_arm_seq, check_sequence, check_sequence_mismatches)
        print(dist)
        sys.exit(1)
        if not dist <= check_sequence_mismatches:
            return ''

    assert target_len > 0, 'invalid (negative) target length for probe %s' % probe
    if trim_primers:
        primer_trim_start = (umi_len + arm_len)
        primer_trim_end = (umi_len + arm_len + target_len)
    else:
        primer_trim_start = (umi_len)
        primer_trim_end = (umi_len + arm_len + target_len + other_arm_len)

    #actually look for the arm sequence instead of blindly cutting
    #find_arm_sequence = 'NNNNNNNNNNNNNNNN'
    #find_arm_sequence = 'GTGTCCGGTGTGAAAAAAAAAA' #at end of read
    if find_arm_sequence is not None:
        assert find_arm_sequence_mismatches is not None
        if True:
            #this is quite a dumb approach, but should work
            search_start = umi_len + arm_len #always search after the first arm
            rest_len = len(read[1]) - search_start

            primer_trim_ends = {}
            for offset in range(rest_len - 5):
                #incomplete matches at the end - don't allow mismatches here anymore
                match_len = min(len(find_arm_sequence), rest_len - offset)
                max_mismatches = find_arm_sequence_mismatches if match_len >= len(find_arm_sequence) else 0

                #print(read[1][search_start + offset:], '@', offset, 'l =', match_len)
                dist = mismatch_count_with_max(read[1][search_start + offset:search_start + offset + match_len], find_arm_sequence[:match_len], max_mismatches = max_mismatches)
                #print(find_arm_sequence[:match_len], '=', dist)

                if dist <= max_mismatches:
                    #keep offset for each possible dist, overwriting later with earlier matches
                    #but only if the match len is at least as long as the current one
                    if not dist in primer_trim_ends or primer_trim_ends[dist][0] <= match_len:
                        primer_trim_ends[dist] = (match_len, offset)

            #print(primer_trim_ends)

            if len(primer_trim_ends) == 0:
                #sys.exit('fail')
                #return ''
                primer_trim_end = None
            else:
                trim_match = primer_trim_ends[min(primer_trim_ends.keys())]
                #print(trim_match)
                primer_trim_end = search_start + trim_match[1]
                #sys.exit('ok')
        else:
            alignment = Bio.pairwise2.align.localxx(read[1][primer_trim_start:], find_arm_sequence, one_alignment_only=True)[0]
            print(find_arm_sequence)
            print(alignment[0])
            print(alignment[1])
            print(alignment[2:])
            _, _, al_score, al_begin, _ = alignment

            match_len = len(read[1]) - al_begin
            print(len(read[1]),  match_len)
            print(al_score - match_len)
            if al_score - match_len <= args.mismatches:
                sys.exit('ok')
                primer_trim_end = primer_trim_start + al_begin
            else:
                sys.exit('fail')
                return ''

    #make sure these are ints (they should be, but the probes.csv might contain floats)
    primer_trim_start = int(primer_trim_start)
    primer_trim_end = int(primer_trim_end)

    #trim seq and quality
    read_seq_trimmed = read[1][primer_trim_start:primer_trim_end]
    read_qual_trimmed = read[2][primer_trim_start:primer_trim_end]
    primer_trimmed_len = len(read_seq_trimmed)

    if trim_primers and find_arm_sequence is None:
        assert len(read_seq_trimmed) <= target_len
        assert len(read_qual_trimmed) <= target_len

    infos = 'pr:Z:%s\tum:Z:%s\tol:i:%d\tts:i:%d\tte:i:%d' % (probe, umi if len(umi) > 0 else 'X', len(read[1]), primer_trim_start, -1 if primer_trim_end is None else primer_trim_end)

    #perform additional quality trimming (after arm trimming, can cut off both 3p and 5p!)
    if quality_trim_threshold is not None:
        quality_trimmed_len, read_seq_trimmed, read_qual_trimmed = quality_trim_read(read_seq_trimmed, read_qual_trimmed,
            p_threshold = quality_trim_threshold,
            phred_base = quality_trim_phred_base)
        quality_trimmed_difference = primer_trimmed_len - quality_trimmed_len
        infos += '\ttq:i:%d' % (quality_trimmed_difference)
        trimmed_len = quality_trimmed_len
    else:
        quality_trimmed_difference = 0
        trimmed_len = primer_trimmed_len

    #check length
    if trim_min_length is not None and trimmed_len < trim_min_length:
        return (None, primer_trim_end, False, 0, trimmed_len)
    else:
        r = ''
        r += '@%s\t%s\n' % (make_extended_read_name(read_id, probe, umi), infos)
        r += '%s\n' % read_seq_trimmed
        r += '+\n'
        r += '%s\n' % read_qual_trimmed

        return (r, primer_trim_end, True, quality_trimmed_difference, trimmed_len)


def parse_read_pairs(
    sample,
    files1,
    files2,
    fastq_out1,
    fastq_out2,
    probes_dict,
    stats_samples,
    stats_reads,
    unknown_arms_directory,
    umi_one, umi_two, mismatches,
    trim_primers, trim_min_length, trim_primers_strict, trim_primers_smart,
    quality_trim_threshold, quality_trim_phred_base,
    allow_multiple_probes,
    consensus_fastqs = None,   
    min_consensus_count = 1,
    min_consensus_fraction = 0.51,
    debug = False, debug_probe = False
):
    """
    Read pairs of reads from FASTQ files and write new FASTQ files with trimmed reads.

    Recognizes probes based on primer sequences, trims off UMI and probe arms,
    optionally performs qualtiy trimming and writes consensus FASTQs based on UMI groups.

    Args:
        sample (str): Sample name.
        files1 (list of str): Input path to paired FASTQ files (Read1)
        files2 (list of str): Input path to paired FASTQ files (Read2)
        fastq_out1 (list of str): Output path to paired, trimmed FASTQ files (Read1)
        fastq_out2 (list of str): Output path to paired, trimmed FASTQ files (Read2)
        probes_dict (dict): Dictionary of probe data (see amplimap.reader.read_probes).
        stats_samples (dict of lists): Sample stats. Can be empty, in which case it will be created.
        stats_reads (list of dicts): Read stats will be appended to this list.
        unknown_arms_directory (str): Directory to write unknown arms files to.
        umi_one, umi_two (int): Length of UMIs on read1 and read2.
        mismatches (int): Number of mismatches to allow in primer sequences.
    """

    if len(stats_samples) == 0:
        counter_columns = [ 'sample', 'files', 'pairs_total', 'pairs_unknown_arms' ]
        if allow_multiple_probes:
            counter_columns.append('pairs_multiple_matching_arms')
        counter_columns.append('pairs_good_arms')
        counter_columns += ['pairs_r1_too_short', 'pairs_r2_too_short']
        if trim_primers_smart:
            counter_columns += ['pairs_r1_second_arm_not_found', 'pairs_r2_second_arm_not_found']
        if quality_trim_threshold is not None:
            counter_columns += ['sum_bp_cut']

        for cc in counter_columns:
            stats_samples[cc] = []

    #convert input primers into expected format
    first_arms_probes, first_arms_seq = zip(*probes_dict['first_primer_5to3'].items()) #split into lists
    #convert to utf-8
    first_arms_seq = tuple([seq.encode('utf-8') for seq in first_arms_seq])
    second_arms_dict = dict([(key, val.encode('utf-8')) for key, val in probes_dict['second_primer_5to3'].items()])
    log.info('Converted input primer sequences for optimised code')

    nDebugged = 0
    counters = dict([(name, 0) for name in stats_samples.keys() if name != 'sample'])
    max_readlen = 0
    umi_per_probe_counter = collections.defaultdict(collections.Counter)
    unknown_arms = (collections.defaultdict(int), collections.defaultdict(int))
    
    tStart = time.time()

    fastqs = None
    if fastq_out1 is not None and fastq_out2 is not None:
        fastqs = {
            'R1': gzip.open(fastq_out1, 'wt'),
            'R2': gzip.open(fastq_out2, 'wt')
        }

    assert len(files1) == len(files2)
    for ixfile in range(len(files1)):
        file1 = files1[ixfile]
        file2 = files2[ixfile]
        log.info('Processing %s and %s (#%d)', file1, file2, ixfile)

        counters['files'] += 1
        opener = gzip.open if file1.endswith('.gz') else open
        with opener(file1, 'rt') as hdl1, opener(file2, 'rt') as hdl2:
            #https://news.open-bio.org/2009/09/25/biopython-fast-fastq/
            for first_read, second_read in zip( Bio.SeqIO.QualityIO.FastqGeneralIterator(hdl1),
                    Bio.SeqIO.QualityIO.FastqGeneralIterator(hdl2) ):
                first_read_str = first_read[1]
                second_read_str = second_read[1]

                counters['pairs_total'] += 1
                if counters['pairs_total'] % 50000 == 0:
                    log.info("processed %d read pairs of which %d have matching arms (%.2f%%) - %.f sec elapsed, %.4f sec/pair, %.1f pairs/sec",
                        counters['pairs_total'],
                        counters['pairs_good_arms'],
                        100.0 * counters['pairs_good_arms'] / counters['pairs_total'],
                        time.time() - tStart,
                        (time.time() - tStart) / counters['pairs_total'],
                        counters['pairs_total'] / (time.time() - tStart))

                debug_this = False
                if debug:
                    if nDebugged >= 5 and not debug_probe:
                        #print(read_details)
                        raise Exception('Stopping (debug)')

                    if random.randint(1, 1000) == 1:
                        debug_this = True
                        nDebugged += 1

                        print('')
                        print('Read #', counters['pairs_total'], ':')
                        print('R1:', first_read_str)
                        print('R2:', second_read_str)

                matches = find_probe_from_reads(
                    first_arms_probes,
                    first_arms_seq,
                    second_arms_dict,
                    first_read_str.encode('utf-8'),
                    second_read_str.encode('utf-8'),
                    umi_one = umi_one,
                    umi_two = umi_two,
                    mismatches = mismatches #,
                    #debug_probe = debug_probe,
                    #debug_this = debug_this
                )

                if len(matches) == 0:
                    if counters['pairs_unknown_arms'] < 10000: #only record the first 10k unknown arms
                        unknown_arms[0][first_read_str[umi_one:umi_one + 20]] += 1
                        unknown_arms[1][second_read_str[umi_two:umi_two + 20]] += 1

                    counters['pairs_unknown_arms'] += 1
                    continue
                else:
                    #if we're still here we actually found a match
                    counters['pairs_good_arms'] += 1

                    #but maybe it's multi?
                    if len(matches) > 1:
                        if allow_multiple_probes:
                            counters['pairs_multiple_matching_arms'] += 1
                        else:
                            log.error('Multiple probes for read pair #%d (%s) -- %s', counters['pairs_total'], first_read[0], str(matches))
                            raise Exception('ABORTING: Multiple probes encountered!')

                    #extract sequences
                    first_umi = first_read_str[:umi_one]
                    second_umi = second_read_str[:umi_two]

                    #extract read id
                    first_id = first_read[0]
                    if ' ' in first_id:
                        first_id = first_id[:first_id.index(' ')]
                    second_id = second_read[0]
                    if ' ' in second_id:
                        second_id = second_id[:second_id.index(' ')]
                    assert len(first_id) > 0, 'Encountered read with missing read name, stopping! Read #%d' % counters['pairs_total']
                    #trim off /1 and /2
                    if first_id.endswith('/1') and second_id.endswith('/2'):
                        first_id = first_id[:-2]
                        second_id = second_id[:-2]
                    assert first_id == second_id, 'Encountered mismatching read IDs between R1 and R2, stopping! Read #%d R1=%s / R2=%s' % (counters['pairs_total'], first_id, second_id)
                    pair_id = first_id

                    for match in matches:
                        probe = match[0]

                        if fastqs is not None:
                            first_len = len(probes_dict['first_primer_5to3'][probe])
                            second_len = len(probes_dict['second_primer_5to3'][probe])

                            r1_trimmed, r1_primer_trim_end, r1_long_enough, r1_quality_trimmed_difference, r1_trimmed_len = make_trimmed_read(
                                pair_id,
                                first_read, probe, probes_dict['target_length'][probe], first_umi+second_umi,
                                umi_one,
                                first_len, second_len,
                                trim_primers = trim_primers, trim_min_length = trim_min_length,
                                check_sequence = None if not trim_primers_strict else probes_dict['second_primer_5to3__rc'][probe],
                                find_arm_sequence = None if not trim_primers_smart else probes_dict['second_primer_5to3__rc'][probe],
                                find_arm_sequence_mismatches = mismatches,
                                quality_trim_threshold = quality_trim_threshold,
                                quality_trim_phred_base = quality_trim_phred_base
                            )

                            r2_trimmed, r2_primer_trim_end, r2_long_enough, r2_quality_trimmed_difference, r2_trimmed_len = make_trimmed_read(
                                pair_id,
                                second_read, probe, probes_dict['target_length'][probe], first_umi+second_umi, 
                                umi_two,
                                second_len, first_len,
                                trim_primers = trim_primers, trim_min_length = trim_min_length,
                                check_sequence = None if not trim_primers_strict else probes_dict['first_primer_5to3__rc'][probe],
                                find_arm_sequence = None if not trim_primers_smart else probes_dict['first_primer_5to3__rc'][probe],
                                find_arm_sequence_mismatches = mismatches,
                                quality_trim_threshold = quality_trim_threshold,
                                quality_trim_phred_base = quality_trim_phred_base
                            )

                            #keep track of max readlen
                            max_readlen = max(max_readlen, r1_trimmed_len)
                            max_readlen = max(max_readlen, r2_trimmed_len)

                            if r1_long_enough and r2_long_enough:
                                fastqs['R1'].write(r1_trimmed)
                                fastqs['R2'].write(r2_trimmed)

                            if quality_trim_threshold is not None:
                                counters['sum_bp_cut'] += r1_quality_trimmed_difference
                                counters['sum_bp_cut'] += r2_quality_trimmed_difference

                            if not r1_long_enough:
                                counters['pairs_r1_too_short'] += 1
                            if not r2_long_enough:
                                counters['pairs_r2_too_short'] += 1

                            if r1_primer_trim_end is None:
                                counters['pairs_r1_second_arm_not_found'] += 1
                            if r2_primer_trim_end is None:
                                counters['pairs_r2_second_arm_not_found'] += 1

                        umi_per_probe_counter[probe][(first_umi+second_umi).encode('utf-8')] += 1

                    if debug_probe and counters['pairs_good_arms'] >= 3:
                        raise Exception('Debugging single probe, only showing first 3')

    #fastq output - make sure file handles are closed properly
    if fastqs is not None:
        for fastq in fastqs.values():
            fastq.close()
        log.info('FASTQ files closed')

    #write table of unknown probes
    for iarms, arms in enumerate(unknown_arms):
        n_unknown_arms_recorded = 0
        for count in arms.values():
            n_unknown_arms_recorded += count

        arms_sorted = sorted(arms.items(), key=operator.itemgetter(1), reverse = True)

        outname = os.path.join(unknown_arms_directory, 'unknown_arms_%s_R%d.csv' % (sample, iarms+1))
        with open(outname, 'w') as unknown_file:
            unknown_file.write('seq_20bp,raw_count,fraction\n')
            for iarm, (arm, count) in enumerate(arms_sorted):
                unknown_file.write('%s,%d,%f\n' % (arm, count, 1.0*count/n_unknown_arms_recorded))

                #only write top 100
                if iarm > 100:
                    break
    log.info('Unknown arm files written')

    #read statistics - with or without umi stats
    if umi_one + umi_two > 0:
        #prepare arrays for consensus calling
        max_umi_group_index = 0
        umi_groups_global = dict()
        umi_probes_global = dict()

        #calculate umi count stats
        for probe, umi_raw_counts in umi_per_probe_counter.items():
            umi_to_group = find_umi_groups(umi_raw_counts)

            #now do a new count
            umi_group_counts = collections.Counter()
            umi_groups_global[probe] = dict()
            for umi, group in umi_to_group.items():
                umi_groups_global[probe][umi] = group + max_umi_group_index
                umi_probes_global[group + max_umi_group_index] = probe.encode('utf-8') #this needs to be bytes

                #we want the number of reads supporting the umi group here, which is
                #the sum of the numbers of reads supporting each grouped umi
                umi_group_counts[group] += umi_raw_counts[umi]
            #update the max group index
            max_umi_group_index = max(umi_groups_global[probe].values()) + 1

            umi_counter_series = pd.Series(list(umi_group_counts.values()))
            umi_counter_series_quantiles = umi_counter_series.quantile([.05, .25, .5, .75, .95])

            stats_reads.append({
              'sample': sample,
              'probe': probe,
              'read_pairs': umi_counter_series.sum(),
              'umis_total': len(umi_counter_series),
              'umis_coverage_min': umi_counter_series.min(), 
              'umis_coverage_mean': umi_counter_series.mean(), 
              'umis_coverage_max': umi_counter_series.max(), 
              'umis_coverage_05pct': umi_counter_series_quantiles[0.05], 
              'umis_coverage_25pct': umi_counter_series_quantiles[0.25], 
              'umis_coverage_median': umi_counter_series_quantiles[0.5], 
              'umis_coverage_75pct': umi_counter_series_quantiles[0.75],
              'umis_coverage_95pct': umi_counter_series_quantiles[0.95],
              'umis_coverage_ge_2': (umi_counter_series > 2).sum(),
              'umis_coverage_ge_3': (umi_counter_series > 3).sum(),
              'umis_coverage_ge_5': (umi_counter_series > 5).sum(),
              'umis_coverage_ge_10': (umi_counter_series > 10).sum()
            })

        #do consensus calling
        if consensus_fastqs is not None:
            #TODO: could actually do this separately for read1/read2 to cut memory req in half
            #TODO: could also handle N differently here!
            consensus_shape = ( 2, max_umi_group_index+1, max_readlen, 5 )

            log.info('Allocating consensus arrays with shape %s', str(consensus_shape))
            base_counts = np.zeros( shape = consensus_shape, dtype = np.unsignedinteger )
            base_max_qual_chars = np.zeros( shape = consensus_shape, dtype = np.unsignedinteger )
            log.info('Consensus arrays allocated')

            #dict to translate indices back to character
            call_to_char = {
                1: ord(b'A'),
                2: ord(b'C'),
                3: ord(b'G'),
                4: ord(b'T'),
                0: ord(b'N'),
            }

            #reopen the trimmed fastq files that we just wrote
            log.info('Re-reading trimmed fastq files...')
            n_reads = 0
            umi_group_read_names = collections.defaultdict(set)

            opener = gzip.open if fastq_out1.endswith('.gz') else open
            with opener(fastq_out1, 'rt') as hdl1, opener(fastq_out2, 'rt') as hdl2:
                #https://news.open-bio.org/2009/09/25/biopython-fast-fastq/
                for first_read, second_read in zip( Bio.SeqIO.QualityIO.FastqGeneralIterator(hdl1),
                        Bio.SeqIO.QualityIO.FastqGeneralIterator(hdl2) ):        
                    n_reads += 1
                    qname = first_read[0].split('\t')[0]
                    read_name, read_probe, read_umi = parse_extended_read_name(qname)
                    read_umi_bytes = read_umi.encode('utf-8')
                    
                    umi_group_index = umi_groups_global[read_probe][read_umi_bytes]
                    if len(umi_group_read_names[umi_group_index]) < 10:
                        umi_group_read_names[umi_group_index].add(read_name.encode('utf-8'))

                    for read_index in [0, 1]:
                        my_calls = np.array(list( (first_read[1] if read_index == 0 else second_read[1]).encode('utf-8') ))
                        my_quals = np.array(list( (first_read[2] if read_index == 0 else second_read[2]).encode('utf-8') ))

                        #TODO: this also gets messy with read len < max read len (need to pad)
                        #for base_index, base_byte in call_to_char.items():
                            #increase count at base_index for all positions (=call_index) where our call is the given base
                            #base_counts[read_index, umi_group_index, my_calls == base_byte, base_index] += 1

                            #TODO: should do this with qual chars as well, but complicated...
                            #my_qual_chars = my_quals[my_calls == base_byte]
                            #base_max_qual_chars[read_index, umi_group_index, np.and(my_calls == base_byte, my_quals > base_max_qual_chars[read_index, umi_group_index, my_calls == base_byte, base_index] , base_index]

                        #still need to loop over this for now
                        for call_index in range(len(my_calls)):
                            my_call = my_calls[call_index]

                            base_index = 0 #== N
                            if my_call == ord(b'A'):
                                base_index = 1
                            elif my_call == ord(b'C'):
                                base_index = 2
                            elif my_call == ord(b'G'):
                                base_index = 3
                            elif my_call == ord(b'T'):
                                base_index = 4

                            base_counts[read_index, umi_group_index, call_index, base_index] += 1
                            if my_quals[call_index] > base_max_qual_chars[read_index, umi_group_index, call_index, base_index]:
                                base_max_qual_chars[read_index, umi_group_index, call_index, base_index] = my_quals[call_index]

            log.info('Processed %d reads', n_reads)
            # for i in range(10):
            #     print('Group ', i)
            #     print(base_counts[0, i])
            #     print(base_max_qual_chars[0, i])
            # sys.exit(1)

            #figure out how many reads we have per base, and what the max per mate is, and finally what the min across mates per umi group is
            n_reads_per_base = base_counts.sum(axis = 3)
            n_reads_per_mate = n_reads_per_base.max(axis=2)
            n_reads_per_group = n_reads_per_mate.min(axis=0)          
            #make sure we have the same value in 
            #call consensus
            max_base = base_counts.argmax(axis = 3)
            max_base_count = base_counts.max(axis = 3)
            #set to N if number of reads supporting this base is < 50% of all reads
            no_consensus = max_base_count < n_reads_per_base * min_consensus_fraction
            log.info('Found %d bases with no consensus, setting to N', no_consensus.sum())
            max_base[no_consensus] = 0
            log.info('Called consensus bases')

            #loop over r1/r2
            log.info('Writing consensus fastq files...')
            assert min_consensus_count >= 1
            for read_index in [0, 1]:
                log.info('Writing to: %s', consensus_fastqs[read_index])
                n_consensus_reads = 0
                with gzip.open(consensus_fastqs[read_index], 'wb') as consensus_fastq:
                    #write one read entry per umi group
                    for umi_group_index in range(max_umi_group_index):
                        umi_group_read_count = n_reads_per_group[umi_group_index]
                        #this should be at least 1 -- skip groups which have a read count of 0, implying that at least one read had no mates
                        if umi_group_read_count >= min_consensus_count:                     
                            bases_with_reads = n_reads_per_base[read_index, umi_group_index] > 0 #should only happen at end!
                            n_bases_disagreement = (n_reads_per_base[read_index, umi_group_index] > max_base_count[read_index, umi_group_index]).sum()

                            #generate name with group index, number of supporting reads and probe name
                            my_name = b'id_%d_reads_%d:::pr_%s' % (umi_group_index, umi_group_read_count, umi_probes_global[umi_group_index])
                            my_tags = b'cr:Z:%s\tcn:i:%d\tcd:i:%d' % (b','.join(umi_group_read_names[umi_group_index]), umi_group_read_count, n_bases_disagreement)
                            #generate seq based on max_base
                            my_seq = bytes([ call_to_char[base_index] \
                                for call_index, base_index in enumerate(max_base[read_index, umi_group_index, bases_with_reads]) ])
                            #generate quality based on max qual for each max_base (NB: 33 assumes phred33 encoding!)
                            my_quals = bytes([ base_max_qual_chars[read_index, umi_group_index, call_index, base_index] if base_index > 0 else 33 \
                                for call_index, base_index in enumerate(max_base[read_index, umi_group_index, bases_with_reads]) ])

                            read = b'@%s\t%s\n%s\n+\n%s\n' % (
                                my_name, my_tags,
                                my_seq,
                                my_quals
                            )

                            consensus_fastq.write(read)
                            n_consensus_reads += 1
                log.info('Wrote %d consensus reads', n_consensus_reads)

            log.info('Consensus FASTQ files written')
    else:
        #we don't have any umis, so for each probe we just have a dict of length one
        for probe, umi_raw_counts in umi_per_probe_counter.items():
            assert len(umi_raw_counts) == 1

            stats_reads.append({
              'sample': sample,
              'probe': probe,
              'read_pairs': sum(umi_raw_counts.values())
            })
    log.info('Read statistics calculated')

    stats_samples['sample'].append(sample)        
    for counter_name, counter_value in counters.items():
        stats_samples[counter_name].append(counter_value)

    log.info('%s done -- sample stats: Total=%d, No Probe=%d (%.1f%%), Multi Probe=%d => Good=%d' % (
        sample, counters['pairs_total'], counters['pairs_unknown_arms'], 100.0*counters['pairs_unknown_arms']/counters['pairs_total'] if counters['pairs_total'] > 0 else -1,
        counters['pairs_multiple_matching_arms'] if 'pairs_multiple_matching_arms' in counters else -1,
        counters['pairs_good_arms']))

    if trim_primers_smart or quality_trim_threshold is not None:
        log.info('%s smart/quality trimming -- %d (%.1f%%) R1 no second arm, %d (%.1f%%) R2 no second arm, %d (%.1f%%) R1 too short, %d (%.1f%%) R2 too short' % (
            sample,
            counters['pairs_r1_second_arm_not_found'] if 'pairs_r1_second_arm_not_found' in counters else -1,
                100.0*counters['pairs_r1_second_arm_not_found']/counters['pairs_good_arms'] if counters['pairs_good_arms'] > 0 and 'pairs_r1_second_arm_not_found' in counters else -1,
            counters['pairs_r2_second_arm_not_found'] if 'pairs_r2_second_arm_not_found' in counters else -1,
                100.0*counters['pairs_r2_second_arm_not_found']/counters['pairs_good_arms'] if counters['pairs_good_arms'] > 0 and 'pairs_r2_second_arm_not_found' in counters else -1,
            counters['pairs_r1_too_short'], 100.0*counters['pairs_r1_too_short']/counters['pairs_good_arms'] if counters['pairs_good_arms'] > 0 else -1,
            counters['pairs_r2_too_short'], 100.0*counters['pairs_r2_too_short']/counters['pairs_good_arms'] if counters['pairs_good_arms'] > 0 else -1))

def output_stats_samples(outname, stats_samples):
    if len(stats_samples) > 0:
        df = pd.DataFrame(stats_samples)
        if 'sum_bp_cut' in stats_samples:
            stats_samples['pairs_mean_bp_cut'] = 1.0 * stats_samples['sum_bp_cut'] / stats_samples['pairs_good_arms']
        df.to_csv(outname, index=False)
    else:
        log.warn('No sample info to output, creating empty file!')
        open(outname, 'a').close()
    print(outname)

def output_stats_reads(outname, stats_reads, umi_one, umi_two):
    if len(stats_reads) > 0:
        columns = [
            'sample',
            'probe',
            'read_pairs'
        ]

        if umi_one + umi_two > 0:
            columns += [
                'umis_total',
                'umis_coverage_min', 
                'umis_coverage_mean', 
                'umis_coverage_max', 
                'umis_coverage_05pct', 
                'umis_coverage_25pct', 
                'umis_coverage_median', 
                'umis_coverage_75pct',
                'umis_coverage_95pct',
                'umis_coverage_ge_2',
                'umis_coverage_ge_3',
                'umis_coverage_ge_5',
                'umis_coverage_ge_10'
            ]

        df = pd.DataFrame(stats_reads, columns = columns)
        df.to_csv(outname, index=False)
    else:
        log.warn('No stats reads to output, creating empty file!')
        open(outname, 'a').close()
    print(outname)

def main():
    log.info('Called with arguments: "%s"', '" "'.join(sys.argv))
    
    #parse the arguments, which will be available as properties of args (e.g. args.probe)
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #specify parameters
    parser.add_argument("-i", "--input", help="input directory containing the fastq files")
    parser.add_argument("-o", "--output", help="working directory", default="out_reads/")
    parser.add_argument("--file-r1", help="name of the read 1 fastq file from a single sample to process")
    parser.add_argument("--files", help="txt file containing the (beginning of) filenames to analyse, one per line (default: files.txt in working directory)")
    parser.add_argument("--files-all", help="add this option to simply read all fastq files in the directory",
                        action='store_true')
    parser.add_argument("--file-mask", help="regular expression to use for finding sample ID etc", default="^(.*_S[0-9]+)_(L[0-9]+)_(R[0-9]+)_([0-9]+).fastq(.gz)?$")
    parser.add_argument("--design", help="CSV file with probes and their arm sequences (default: probes.csv in working directory)")

    parser.add_argument("--index", help="index of the file to process (starting from 1)",
                        type=int)
    parser.add_argument("--output-fastq", help="generate fastq files with the trimmed reads for each sample and probe",
                        action='store_true')
    parser.add_argument("--trim-primers", help="for fastq output, trim off primers (probe arms) in addition to umis",
                        action='store_true')
    #parser.add_argument("--trim-primers-strict", help="TODO: for fastq output, only trim primers if second primer is in expected location",
    #                    action='store_true')
    parser.add_argument("--trim-primers-smart", help="UNTESTED: search for second primer sequence in first read and vice versa to find point from where to trim",
                        action='store_true')
    parser.add_argument("--quality-trim-threshold", help="quality trimming p-value threshold - None = no trimming, values around 0.01 recommended",
                        default=None, type=float)
    parser.add_argument("--quality-trim-phred-base", help="quality trimming phred base",
                        default=33, type=int)
    parser.add_argument("--trim-min-length", help="minimum length of final read after flexible trimming",
                        default=20, type=int)

    parser.add_argument("--umi-one", help="length of the UMI at the beginning of read one, 0=none", type=int, default=5)
    parser.add_argument("--umi-two", help="length of the UMI at the beginning of read two, 0=none", type=int, default=5)
    parser.add_argument("--mismatches", help="maximum number of mismatches between expected and observed arm sequence",
                        type=int, default=2)
    parser.add_argument("--allow-missing-samples", help="continue even if we find fewer samples than expected",
                        action='store_true')
    parser.add_argument("--allow-multiple-probes", help="continue even if we find multiple probes matching a single read (will lead to duplication of reads)",
                        action='store_true')
    parser.add_argument("--allow-duplicate-probes", help="continue even if we find duplicate rows for the same probe",
                        action='store_true')
    parser.add_argument("--allow-wrong-columns", help="continue even if the column names in the probe design table are incorrect",
                        action='store_true')

    parser.add_argument("--check-probes", help="Just check the probes file and quit",
                        action='store_true')
    parser.add_argument("--debug", help="For testing: enable debug output for first few reads",
                        action='store_true')
    parser.add_argument("--debug-probe", help="For testing: probe to debug")
    parser.add_argument("--debug-only-first-100k", help="For testing: only use the first 100,000 reads",
                        action='store_true')
    parser.add_argument("--no-rc-arm-r1", help="For testing: do NOT reverse complement first probe",
                        action='store_true')
    parser.add_argument("--rc-arm-r2", help="For testing: DO reverse complement the second probe",
                        action='store_true')
    args = parser.parse_args()
    args.trim_primers_strict = None #TODO

    if not os.path.isdir(args.output):
        raise Exception('Working directory does not exist: %s' % args.output)

    if not args.design:
        #try probes.csv first, followed by design.csv
        args.design = os.path.join(args.output, 'probes.csv')

    if not args.files_all:
        if not args.files:
            args.files = os.path.join(args.output, 'files.txt')

    if args.output_fastq:
        fastq_directory = os.path.join(args.output, 'fastq')
        try:
            os.makedirs(fastq_directory)
        except OSError as e:
            #log.info('Got OSError while trying to create FASTQ directory %s: %s' % (fastq_directory, str(e)))
            pass

    #read design
    design = read_new_probe_design(args.design, recalculate_targets = False, allow_wrong_columns = args.allow_wrong_columns, allow_duplicate_probes = args.allow_duplicate_probes)
    #if we made it till here and are only checking probes, we can quit
    if args.check_probes:
        return 0

    # #reverse-complement the first arm sequence
    # if not args.no_rc_arm_r1:
    #     log.info('Reverse-complementing read 1 probe sequence in the design table to match read orientation...')
    #     #print(design.head(3))

    #     #swap columns
    #     old = design['first_arm']
    #     design.loc[:, 'first_arm'] = design['first_arm__rc']
    #     design.loc[:, 'first_arm__rc'] = old

    # #rc the second arm sequence
    # if args.rc_arm_r2:
    #     log.info('Reverse-complementing read 2 probe sequencein the design table to ????? (does not match 5p -> 3p orientation)...')
    #     #print(design.head(3))

    #     #swap columns
    #     old = design['second_arm']
    #     design.loc[:, 'second_arm'] = design['second_arm__rc']
    #     design.loc[:, 'second_arm__rc'] = old

    if args.debug_probe:
        log.warn('Only debugging probe: %s', args.debug_probe)
        design = design[ design['id'] == args.debug_probe ]
        assert len(design) == 1                    

    #much faster to access than using DataFrame.loc every time
    probes_dict = design.to_dict()

    #figure out which filenames we need to read
    sample_file_regex = args.file_mask
    files_to_process = []
    input_dir = args.input
    if args.file_r1:
        input_dir = os.path.dirname(args.file_r1)
        input_file = os.path.basename(args.file_r1)
        if re.search(sample_file_regex, input_file):
            n_samples_expected = 1
            files_to_process = [input_file, input_file.replace('_R1', '_R2')]
            log.info('Only processing single sample. Files: %s' % (str(files_to_process)))
        else:
            log.error('Sample file name does not match expected pattern: %s', args.file_r1)
    else:
        file_prefixes = None
        n_samples_expected = None
        if args.files:
            log.info('Looking for file name prefixes in %s' % args.files)
            with open(args.files, 'Ur') as f:
                file_prefixes = [l for l in [l.strip() for l in f.readlines()] if len(l) > 0]

            #remove non-unique
            file_prefixes = set(file_prefixes)

            if len(file_prefixes) == 1 and file_prefixes.pop() == '*':
                log.info('Searching ALL files because files.txt contained only an asterisk!')
                file_prefixes = None
            else:
                n_samples_expected = len(file_prefixes)
                log.info('Searching for %d file prefixes: %s', len(file_prefixes), ', '.join(file_prefixes))

        log.info('Searching for fastq(.gz) files in %s', args.input)
        nFiles = 0
        for file in os.listdir(args.input):  
            nFiles += 1

            match = re.search(sample_file_regex, file)

            found_prefix = True
            if file_prefixes is not None:
                found_prefix = False
                for fp in file_prefixes:
                    if file.startswith(fp):
                        found_prefix = True
                        break

            if match and found_prefix:
                files_to_process.append(file)
            elif not found_prefix:
                pass #we just didn't find the prefix we wanted
            else:
                log.warn('Unexpected file name: %s' % file)

        log.info('Directory contains %d files, %d match pattern' % (nFiles, len(files_to_process)))

    #now process filenames and make a list of samples
    sample_fastqs = {}
    for file in files_to_process:
        match = re.search(sample_file_regex, file)

        sample, lane, read, _, _ = match.groups()

        if not sample in sample_fastqs:
            sample_fastqs[sample] = {}

        while read in sample_fastqs[sample]:
            log.warn('Duplicate sample/read index combination! file=%s sample=%s existing=%s' % (file, sample, str(sample_fastqs[sample])))
            sample += '_DUP'
            if not sample in sample_fastqs:
                sample_fastqs[sample] = {}

        # #make sure we match the sample number
        # other_read = 'R2' if read == 'R1' else 'R1' if read == 'R2' else sys.exit('invalid read?!')
        # if other_read in sample_fastqs[sample]:
        #     assert sample_fastqs[sample][other_read]['Sample'] == sample

        sample_fastqs[sample][read] = {
            'File': file,
            #'Sample': sample,
            'Lane': lane
        }

    log.info('Found %d samples: %s', len(sample_fastqs), ', '.join([
        '%d=%s' % (ix+1, sample) for ix, sample in enumerate(sample_fastqs.keys())]))

    if n_samples_expected is not None:
        if len(sample_fastqs) == n_samples_expected:
            log.info('Found the expected number of samples.')
        else:
            log.error('Did not find the expected number of samples! Incorrect/missing file prefix?')

        if not args.allow_missing_samples:
            assert len(sample_fastqs) == n_samples_expected, 'set --allow-missing-samples to allow this!'

    output_code = 'merged'
    if args.file_r1:
        (sample, info) = list(sample_fastqs.items())[0]
        output_code = 'sample_'+sample #+'_'+info['R1']['Sample']
        #we already have our one file
    elif args.index:
        output_code = 'index__'+str(args.index)

        #this is a bit hacky, but allows us to keep the loop code below
        my_item = sample_fastqs.items()[args.index - 1]
        sample_fastqs = dict([my_item])
        log.info('Only processing sample #%d', args.index)
    else:
        log.warn('No index provided, reading ALL samples!')

    if len(sample_fastqs) == 0:
        raise Exception('No files found!')

    log.info('Using output code: %s', output_code)

    log.info('Processing %d samples: %s', len(sample_fastqs), ', '.join(sample_fastqs.keys()))

    stats_samples = collections.OrderedDict()
    stats_reads = []

    nSamples = 0
    for (sample, read_files) in sample_fastqs.items():
        nSamples += 1
        if re.search("Undetermined", sample): 
            log.info('Ignoring sample: %s' % sample)
            continue

        if not 'R1' in read_files or not 'R2' in read_files:
            raise Exception('Did not get exactly one read file for each R1 and R2 for sample: %s' % (sample, ))

        log.info('%s (%d/%d) - Read 1 = %s; Read 2 = %s' % (sample, nSamples, len(sample_fastqs),
            read_files['R1']['File'], read_files['R2']['File']))

        file1 = os.path.join(input_dir, read_files['R1']['File'])
        file2 = os.path.join(input_dir, read_files['R2']['File'])

        fastq_out1 = None
        fastq_out2 = None
        if args.output_fastq:
            fastq_out1 = os.path.join(fastq_directory, '%s__MIP_TRIMMED__%s_R%d_001.fastq.gz' % (sample, read_files['R1']['Lane'], 1))
            fastq_out2 = os.path.join(fastq_directory, '%s__MIP_TRIMMED__%s_R%d_001.fastq.gz' % (sample, read_files['R2']['Lane'], 2))

        min_consensus_fraction = 0.51 #just set this here for now, this code path isn't used anymore anyway
        parse_read_pairs(
            sample,
            [file1],
            [file2],
            fastq_out1,
            fastq_out2,
            probes_dict,
            stats_samples,
            stats_reads,
            args.output,
            args.umi_one, args.umi_two, args.mismatches,
            args.trim_primers, args.trim_min_length, args.trim_primers_strict, args.trim_primers_smart,
            args.quality_trim_threshold, args.quality_trim_phred_base,
            args.allow_multiple_probes,
            None,
            min_consensus_count,
            min_consensus_fraction,
            args.debug, args.debug_probe
        )

    outname = os.path.join(args.output, 'stats_samples__%s.csv' % output_code)
    output_stats_samples(outname, stats_samples)

    outname = os.path.join(args.output, 'stats_reads__%s.csv' % output_code)
    output_stats_reads(outname, stats_reads, args.umi_one, args.umi_two)

    print('DONE')
    return 0

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    sys.exit(main())