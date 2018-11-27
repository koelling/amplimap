# -*- coding: utf-8 -*-
"""
This module contains methods for generating pileups (per-base allele counts) using pysam/samtools.
"""

#python 3 compat
#http://python-future.org/compatible_idioms.html
from __future__ import print_function

import sys
import os
import re

import argparse

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
log.addHandler(sh)

import time

#use pandas to read the CSV file and write output files
import pandas as pd
#for defaultdict + sorting
import collections, itertools
import operator

#for pileup
import pysam
#for getting ref base 
import pyfaidx
#for debugging umi distance
import distance

#functions for reading input files
from .reader import read_snps_txt, read_new_probe_design, read_targets
from .common import parse_extended_read_name, find_umi_groups

#need these in various places
call_types_bases = ['A', 'C', 'G', 'T']
call_types = call_types_bases + ['INS', 'DEL', 'SKIP']

def get_al_mate_starts(al: pysam.AlignedSegment):
    """Get set of mate starts for AlignedSegment, always giving read1 first and read2 second"""
    return (al.reference_start, al.next_reference_start) if al.is_read1 else (al.next_reference_start, al.reference_start)

def get_pileup_row(chrom, pos_0, raw_coverage = 0, target_id = None, target_type = None, ref = None, validate_probe_targets = False):
    """Get ordered dict of columns for pileup table."""
    row = collections.OrderedDict()
    row['chr'] = chrom
    row['pos'] = pos_0 + 1
    if target_id is not None:
        row['target_id'] = target_id
    if target_type is not None:
        row['target_type'] = target_type
    if ref:
        row['ref'] = ref[chrom][pos_0:pos_0+1]
        row['alts'] = ''
    row['raw_coverage'] = raw_coverage

    row['umi_groups'] = 0
    row['unique_raw_umis'] = 0
    row['unique_mate_starts'] = 0
    row['group_below_min'] = 0
    row['group_no_majority'] = 0

    row['group_no_call'] = 0
    row['number_called'] = 0
    row['number_called_hq'] = 0
    
    row['filtered_umi_n'] = 0
    row['filtered_low_quality'] = 0
    row['filtered_flagged'] = 0
    row['filtered_softclipped'] = 0
    if validate_probe_targets:
        row['outside_probe_target'] = 0
    row['probes'] = ''

    row['overlapping_mates'] = 0

    #row['count_indel_pos'] = 0
    #row['count_indel_neg'] = 0
    for call in call_types:
        row['count_%s' % call] = 0
    for call in call_types:
        row['count_hq_%s' % call] = 0

    return row   

class PileupRowFilterException(Exception):
    """Raised when row fails a pileup filter, with filter_column being the column to count this in."""
    def __init__(self, filter_column, skip_read_pair = True):
        self.filter_column = filter_column
        self.skip_read_pair = skip_read_pair

class PileupGroupFilterException(Exception):
    """Raised when UMI group fails a pileup filter, with filter_column being the column to count this in."""
    def __init__(self, filter_column):
        self.filter_column = filter_column

def process_pileup_read(
        pr,
        probes_dict,
        reference_name,
        reference_pos_0,
        read_metadata,
        min_mapq,
        ignore_groups,
        group_with_mate_positions,
        validate_probe_targets,
        filter_softclipped,
        no_probe_data,
        debug):
    """
    Process single pileup read.
    """
    #http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment
    #http://pysam.readthedocs.io/en/latest/api.html#pysam.PileupColumn
    #http://pysam.readthedocs.io/en/latest/api.html#pysam.PileupRead

    #filter out low mapping quality
    if pr.alignment.mapping_quality < min_mapq:
        raise PileupRowFilterException('filtered_low_quality', skip_read_pair = False) #we don't want to skip the entire pair (for now, to keep compat with old branch)

    #ignore non-proper pair, qc fail, secondary alignments, ...
    #removed: (not pr.alignment.is_proper_pair) or
    #removed: (pr.alignment.is_secondary) or (for secondary alignments from bowtie2, which we would want to keep. hopefully should not matter otherwise since bwa never gives us secondary (?))
    if (pr.alignment.is_qcfail) or \
    (pr.alignment.mate_is_unmapped) or \
    (pr.alignment.is_supplementary):
        raise PileupRowFilterException('filtered_flagged')

    #filter out reads with softclipped bases
    if filter_softclipped:
        if 'S' in pr.alignment.cigarstring:
            raise PileupRowFilterException('filtered_softclipped', skip_read_pair = False) #we don't want to skip the entire pair!

    #get read metadata        
    read_name = pr.alignment.query_name
    mate_starts = get_al_mate_starts(pr.alignment) if group_with_mate_positions else None
    if ignore_groups:
        #we haven't pre-parsed the read yet, so let's do it now
        read_umi = None

        #parse info from read name (should be bowtie2 compatible)
        if no_probe_data:
            read_probe = None
        else:
            _, read_probe, _ = parse_extended_read_name(read_name)
    else:
        #get read metadata from before (will raise if we didn't see this before, but we should always have)
        alignment_tuple = (read_name, mate_starts)
        read_name, read_probe, read_umi, _ = read_metadata[alignment_tuple]

    #filter out reads where the original UMI contains an N
    if not ignore_groups and read_umi is not None:
        if b'N' in read_umi:
            raise PileupRowFilterException('filtered_umi_n')

    if validate_probe_targets and read_probe is not None:
        if reference_pos_0 < probes_dict['target_start_0'][read_probe] \
        or reference_pos_0 >= probes_dict['target_end'][read_probe] \
        or reference_name != probes_dict['chr'][read_probe]:
            raise PileupRowFilterException('outside_probe_target', skip_read_pair = False)

    #figure out call and phred qual
    try:
        my_call = None
        my_phred = None
        if pr.indel > 0: #this peaks ahead!
            my_call = 'INS' #note: this replaces the basecall at the base preceding the INS!
            my_phred = 40
        elif pr.is_del: #pr.indel < 0 peaks ahead, so we really just want to use pr.is_del here!
            my_call = 'DEL'
            my_phred = 40 #may want to use query_position_or_next here to get qual of next base?
        elif pr.is_refskip:
            my_call = 'SKIP'
            my_phred = 40
        else:
            pos = pr.query_position
            query_sequence = pr.alignment.query_sequence

            if pos is None:
                print(pr)
                print(pr.indel)
                print(pr.is_del)
                raise Exception('query_position is None?!')
            elif pos >= len(query_sequence):
                raise Exception('pos = %d but len = %d' % (pos, len(query_sequence)))
            elif pos < 0:
                raise Exception('query position < 0: %d' % pos)
            
            my_call = query_sequence[pos]
            my_phred = pr.alignment.query_qualities[pos]
    except:
        print(reference_name, reference_pos_0)
        my_call = 'ERROR'
        raise

    return my_call, my_phred, read_name, read_probe, read_umi, mate_starts

def record_read_in_group(read_calls, my_call, my_phred, my_umi, read_name):
    """
    Add call data from read to read_calls dictionary of name -> call.
    """
    my_call_data = [my_call, my_phred, my_umi] #, my_start

    #do we already have a call from this read pair?
    if read_name in read_calls:
        assert my_umi == read_calls[read_name][2] #should always match!

        # #always take lowest start (or overwrite if None, although in that case we're probably None too)
        # if read_calls[read_name][3] is None or my_call_data[3] > read_calls[read_name][3]:
        #     my_call_data[3] = read_calls[read_name][3]

        #resolve overlap by taking higher-quality call
        if read_calls[read_name][1] < my_phred:
            read_calls[read_name] = my_call_data
        #handle case where we have the same quality for two different calls (get suspicious)
        elif read_calls[read_name][1] == my_phred and read_calls[read_name][0] != my_call:
            read_calls[read_name][0] = 'N'

        return True
    else:
        read_calls[read_name] = my_call_data
        return False

def get_group_consensus(group_calls, min_consensus_count = 1, min_consensus_fraction = 0.51, ignore_groups = False, debug = False):
    """
    Calculate consensus call, count and phred for UMI group.
    """
    group_call_counts = collections.defaultdict(int)
    group_calls_maxphred = collections.defaultdict(int)
    
    #count how often we saw each call
    for (my_call, my_phred, _) in group_calls:
        group_call_counts[my_call] += 1
        group_calls_maxphred[my_call] = max(group_calls_maxphred[my_call], my_phred)

    #find call with highest count
    #note: as long as min_consensus_fraction is >0.5 we should never end up with a duplicate call
    #in the end, since if we have a duplicate its fraction has to be <= 50%
    group_consensus_call, group_consensus_count = None, 0    
    for my_call, my_count in group_call_counts.items():
        if my_count > group_consensus_count:
            group_consensus_call = my_call
            group_consensus_count = my_count
    group_consensus_phred = group_calls_maxphred[group_consensus_call]

    if debug:
        print('Group', '( len =', len(group_calls), ')', ' calls:',
            ', '.join(['%dx %s (q%d)' % (group_call_counts[k], k, group_calls_maxphred[k]) for k in group_call_counts.keys()]))
        print(group_call_counts.most_common(1))

    #filter out non-consensus calls by raising a PileupGroupFilterException
    if not ignore_groups and (group_consensus_count < min_consensus_count):
        if debug:
            print('-> Below min count -- ', group_consensus_call, ' (', group_consensus_count, 'x ), max Q =', group_consensus_phred)
        raise PileupGroupFilterException('group_below_min')
    elif not ignore_groups and (group_consensus_count < len(group_calls) * min_consensus_fraction):
        if debug:
            print('-> No majority call -- ', group_consensus_call, ' (', group_consensus_count, 'x ), max Q =', group_consensus_phred)
        raise PileupGroupFilterException('group_no_majority')
    elif group_consensus_call == 'N':
        if debug:
            print('-> No base call (N) -- ', group_consensus_call, ' (', group_consensus_count, 'x ), max Q =', group_consensus_phred)
        raise PileupGroupFilterException('group_no_call')
    else:
        if debug:
            print('-> Consensus: ', group_consensus_call, ' (', group_consensus_count, 'x ), max Q =', group_consensus_phred)

    return group_consensus_call, group_consensus_count, group_consensus_phred   

def process_pileup_row(
        row,
        seen_probes,
        call_groups,
        snps_dict,
        ignore_groups,
        min_consensus_count,
        min_consensus_fraction,
        min_baseq,
        ref,        
        debug = False):
    """
    Generate a pileup table row from processed read data and calculate some additional stats and annotation.
    """
    #make a list of probes
    row['probes'] = ';'.join(sorted(seen_probes))

    phreds = list()
    call_phreds = collections.defaultdict(list)

    #process dictionary of best calls per read pair
    nonref_call_groups = []
    nonref_call_supporting_reads = []
    row['umi_groups'] = len(call_groups)
    for (group_name, my_group) in call_groups.items():
        try:
            #find consensus for this group
            group_consensus_call, group_consensus_count, group_consensus_phred = get_group_consensus(
                my_group.values(),
                min_consensus_count = min_consensus_count,
                min_consensus_fraction = min_consensus_fraction,
                ignore_groups = ignore_groups,
                debug = debug
            )

            #record phreds
            phreds.append(group_consensus_phred)
            call_phreds[group_consensus_call].append(group_consensus_phred)

            #make lists of the groups supporting the nonref call, and their number of consensus reads
            if ref and group_consensus_call != row['ref']:
                if len(nonref_call_groups) < 5:
                    nonref_call_groups.append(str(group_name))
                if len(nonref_call_supporting_reads) < 25:
                    nonref_call_supporting_reads.append(group_consensus_count)

            #record the call
            row['number_called'] += 1
            row['count_%s' % group_consensus_call] += 1

            if group_consensus_phred >= min_baseq:
                row['number_called_hq'] += 1
                row['count_hq_%s' % group_consensus_call] += 1
        except PileupGroupFilterException as ex:
            #if the row is filtered, we raise a PileupRowFilterException with message being the corresponding filter column
            row[ex.filter_column] += 1
    
    #mean phred score for position 
    if len(phreds) > 0:
        row['phred_mean'] = 1.0 * sum(phreds) / len(phreds)    
    else:
        row['phred_mean'] = None

    #mean phred score per count
    for call in call_types_bases:
        if call in call_phreds:
            row['phred_%s' % call] = 1.0 * sum(call_phreds[call]) / len(call_phreds[call])
        else:
            row['phred_%s' % call] = None

    #extract call counts again
    call_counts = [row['count_hq_%s' % call] for call in call_types]

    #ref/alts count
    if ref:
        if row['ref'] in call_types:
            row['ref_hq_count'] = row['count_hq_%s' % row['ref']]
            row['nonref_hq_count'] = sum([row['count_hq_%s' % call] for call in call_types if call != row['ref']])
            row['alts'] = ';'.join([call_types[i] for i in range(len(call_types)) if call_types[i] != row['ref'] and call_counts[i] > 0])
        else:
            row['ref_hq_count'] = None
            row['nonref_hq_count'] = None
            row['alts'] = None
    
        #make a list of umi groups with the nonref call
        row['nonref_call_groups_5'] = ';'.join(nonref_call_groups)
        if row['nonref_hq_count'] is not None and row['nonref_hq_count'] > 5:
            row['nonref_call_groups_5'] += ';...'
        #make a sorted list of umi group converage
        row['nonref_call_supporting_reads_25'] = ';'.join([str(x) for x in sorted(nonref_call_supporting_reads)])
        if row['nonref_hq_count'] is not None and row['nonref_hq_count'] > 25:
            row['nonref_call_supporting_reads_25'] += ';...'

    #snp ref/alt count
    if snps_dict is not None and ('snp_has_genotypes' in snps_dict):
        if snps_dict['snp_has_genotypes'][row['target_id']]:
            row['snp_ref_hq_count'] = row['count_hq_%s' % snps_dict['snp_ref'][row['target_id']]]
            row['snp_alt_hq_count'] = row['count_hq_%s' % snps_dict['snp_alt'][row['target_id']]]
        else:
            row['snp_ref_hq_count'] = None
            row['snp_alt_hq_count'] = None
    
    #maj/second/nonmaj count
    row['maj_hq_count'] = max(call_counts)
    if row['maj_hq_count'] > 0:
        #if we have a maj count we need to have called something
        assert row['number_called_hq'] > 0

        row['maj_hq_percent'] = 100.0 * row['maj_hq_count'] / row['number_called_hq']
        row['maj_hq_calls'] = ';'.join([call_types[i] for i in range(len(call_types)) if call_counts[i] == row['maj_hq_count']])
        row['maj_hq_phreds'] = ';'.join([str(row['phred_%s' % call_types_bases[i]]) for i in range(len(call_types_bases)) if call_counts[i] == row['maj_hq_count']])
        
        row['non_maj_hq_count'] = sum([cc for cc in call_counts if cc != row['maj_hq_count']])
        #row['non_maj_hq_percent'] = 100.0 * row['non_maj_hq_count'] / row['number_called_hq']
        
        non_maj = [i for i in range(len(call_types)) if call_counts[i] > 0 and call_counts[i] != row['maj_hq_count']]
        non_maj_counts = [call_counts[i] for i in non_maj]
        non_maj_phreds = list(itertools.chain.from_iterable(( call_phreds[call_types[i]] for i in non_maj )))            
        row['non_maj_hq_calls'] = ';'.join([call_types[i] for i in non_maj])
        #row['non_maj_mean_phreds'] = ';'.join([str(row['phred_%s' % call_types[i]]) if 'phred_%s' % call_types[i] in row else 'NA' for i in non_maj])
        #row['non_maj_hq_phreds_mean'] = 1.0 * sum(non_maj_phreds) / len(non_maj_phreds) if len(non_maj_phreds) > 0 else None

        if len(non_maj_counts) > 0:
            row['second_hq_count'] = max(non_maj_counts)
            row['second_hq_percent'] = 100.0 * row['second_hq_count'] / row['number_called_hq']

            if row['second_hq_count'] > 0:
                second = [i for i in range(len(call_types)) if call_counts[i] == row['second_hq_count']]
                row['second_hq_calls'] = ';'.join([call_types[i] for i in second])
            else:
                row['second_hq_calls'] = None
        else:
            row['second_hq_count'] = None
            row['second_hq_percent'] = None
            row['second_hq_calls'] = None
    else:
        row['maj_hq_percent'] = None
        row['maj_hq_calls'] = None
        row['maj_hq_phreds'] = None                    
        row['non_maj_hq_count'] = None
        row['non_maj_hq_calls'] = None                        
        row['second_hq_count'] = None
        row['second_hq_percent'] = None
        row['second_hq_calls'] = None

    return row

def process_pileup_base(
    ref,
    probes_dict,
    targets_dict,
    snps_dict,
    region_index,
    pileup_base,
    min_consensus_count,
    min_consensus_fraction,
    min_mapq,
    min_baseq,
    ignore_groups,
    group_with_mate_positions,
    validate_probe_targets,
    filter_softclipped,
    no_probe_data,
    read_metadata,
    filtered_pair_qnames,
    umi_to_group,
    debug
):
    """
    Process a single basepair of pileup, generating a pileup table row.
    """
    #pysam is zero-based! http://pysam.readthedocs.io/en/latest/glossary.html#term-region
    row = get_pileup_row(pileup_base.reference_name, pileup_base.reference_pos,
        pileup_base.nsegments,
        targets_dict['id'][region_index] if targets_dict is not None and region_index < len(targets_dict['id']) else None,
        targets_dict['type'][region_index] if targets_dict is not None and region_index < len(targets_dict['id']) else None,
        ref,
        validate_probe_targets = validate_probe_targets)

    #loop over pileups, make a dict of the best-phred call per group per read id
    call_groups = collections.defaultdict(dict)
    seen_probes = set()
    seen_raw_umis = set()
    seen_mate_starts = set()
    for pr in pileup_base.pileups:
        #if we have already filtered this read just skip it immediately
        if pr.alignment.query_name in filtered_pair_qnames:
            row[filtered_pair_qnames[pr.alignment.query_name]] += 1
        else:
            try:
                read_call, read_phred, read_name, read_probe, read_umi, mate_starts = process_pileup_read(
                    pr,
                    probes_dict,
                    pileup_base.reference_name,
                    pileup_base.reference_pos,
                    read_metadata,
                    min_mapq = min_mapq,
                    ignore_groups = ignore_groups,
                    group_with_mate_positions = group_with_mate_positions,
                    validate_probe_targets = validate_probe_targets,
                    filter_softclipped = filter_softclipped,
                    no_probe_data = no_probe_data,
                    debug = debug)

                #record which probes we have seen
                if read_probe is not None and not read_probe in seen_probes:
                    seen_probes.add(read_probe)

                #record which mate starts we have seen
                seen_mate_starts.add(mate_starts)

                #figure out umi (needed for grouping)
                if not ignore_groups:
                    #determine group
                    my_group = umi_to_group[mate_starts][read_umi]
                    seen_raw_umis.add(read_umi)
                else:
                    #just group by read name (every read is its own group)
                    my_group = read_name
         
                #record read straight in umi group - should be fast
                is_overlap = record_read_in_group(call_groups[my_group], read_call, read_phred, read_umi, read_name)
                if is_overlap:
                    row['overlapping_mates'] += 1        
            except PileupRowFilterException as ex:
                #if the row is filtered, we raise a PileupRowFilterException with message being the corresponding filter column
                row[ex.filter_column] += 1

                #this tells us whether we should always skip this read, or just skip this position of this read
                if ex.skip_read_pair:
                    filtered_pair_qnames[pr.alignment.query_name] = ex.filter_column

    row['unique_raw_umis'] = len(seen_raw_umis)
    row['unique_mate_starts'] = len(seen_mate_starts)

    row = process_pileup_row(
        row,
        seen_probes,
        call_groups,
        snps_dict = snps_dict,
        ignore_groups = ignore_groups,
        min_consensus_count = min_consensus_count,
        min_consensus_fraction = min_consensus_fraction,
        min_baseq = min_baseq,
        ref = ref,
        debug = debug
        )

    return row


def aggregate(folder, snps_file, ref, generate_calls):
    log.info('Aggregating tables in %s', folder)

    #load SNPs to add SNP info
    snps = None
    if snps_file is not None and len(snps_file) > 0:
        snps = read_snps_txt(snps_file)

    agg = None
    agg_long = None
    coverage_long = None
    coverage_agg = None

    files = [file for file in os.listdir(folder) if file.endswith('.pileup.csv')]
    log.info('Found %d files', len(files))
    assert len(files) > 0

    files.sort()
    ixfile = 0
    for file in files:
        ixfile += 1
        sample = file
        sample = sample.replace('.pileup.csv', '')
        sample = sample.replace('_L001_R1_001', '')
        sample = sample.replace('__MIP_TRIMMED_', '')

        try:
            df = pd.read_csv(os.path.join(folder, file), low_memory=False, #otherwise get type warning, probably mis-detection of T/F or chr
                dtype = {'chr': str, 'pos': int, 'ref': str, 'alts': str})
            log.info('%d/%d: Found %d pileup rows in %s -> %s', ixfile, len(files), len(df), file, sample)
            if len(df) > 0:
                dup = df.duplicated(subset = ['chr', 'pos'], keep = 'first')
                if dup.any():
                    log.error('Found %d duplicate chr/pos pairs in file, cancelling!', dup.sum())
                    raise Exception()

                #add sample name and add to long table
                df['sample'] = sample

                #set filter status to NA
                df['filter'] = '--'
                    
                #columns to join on
                cols = None
                #columns we actually want to output in the end
                out_cols = ['probes'] #always get the probes
                #columns that contain the actual counts
                count_cols = None #columns to output, but not join on

                #add snp info
                if snps is not None:
                    df = snps.merge(df, on = ('chr', 'pos'), how = 'left')
                    df.loc[df['snp_ref'].isnull(), 'snp_ref'] = ''
                    df.loc[df['snp_alt'].isnull(), 'snp_alt'] = ''

                    #use the snp columns as our baseline
                    cols = snps.columns.tolist()
                    count_cols = ['snp_ref_hq_count', 'snp_alt_hq_count']

                    #add filter stats
                    df['filter'] = ''

                    if ref:
                        #compare ref/observed alts to snp alleles (empty snp_ref/alt = pass)
                        df.loc[(df['snp_ref'] != '') & (df['ref'] != df['snp_ref']), 'filter'] += 'REF_mismatch;'
                        df.loc[(df['snp_alt'] != '') & (~df['alts'].isnull()) \
                            & (df['alts'] != df['snp_alt']), 'filter'] += 'unexpected_alt;'
                    else:
                        #compare observed alleles to snp alleles
                        df['filter'] = ['.' if pd.Series([row.maj_hq_calls] + row.second_hq_calls.split(';')).isin(
                            ['', row.snp_ref, row.snp_alt]).all() else 'unknown_allele' for row in df.itertuples()]

                        #warn if multiple major/second
                        df.loc[df['maj_hq_calls'].str.contains(';'), 'filter'] += 'multiple_major;'
                        df.loc[df['second_hq_calls'].str.contains(';'), 'filter'] += 'multiple_second;'

                    df.loc[df['filter'] == '', 'filter'] = '.' #set empty to dot
                else:
                    #columns to include in every row
                    cols = ['chr', 'pos']
                    #add target id, if we have it
                    if 'target_id' in df.columns:
                        cols.append('target_id')
                    if 'target_type' in df.columns:
                        cols.append('target_type')

                    if ref:
                        count_cols = ['ref_hq_count', 'nonref_hq_count']
                        cols.append('ref') #keep the ref column
                        cols.append('alts') #keep the alt column
                    else:
                        df.loc[df['maj_hq_calls'].isnull(), 'maj_hq_calls'] = ''
                        df.loc[df['second_hq_calls'].isnull(), 'second_hq_calls'] = ''

                        count_cols = ['maj_hq_count', 'second_hq_count']

                #also get the total number of called umis
                count_cols.append('number_called_hq')

                #make sure these are properly zero-filled
                for col in count_cols:
                    df[col] = df[col].fillna(0).astype(int)

                #add counts to output
                out_cols += count_cols

                #calculate fractions
                df[count_cols[0]+'_fraction'] = 1.0 * df[count_cols[0]] / df['number_called_hq']
                df[count_cols[1]+'_fraction'] = 1.0 * df[count_cols[1]] / df['number_called_hq']

                #either make this the first one or add to existing
                if agg_long is None:
                    agg_long = df
                else:
                    agg_long = pd.concat([agg_long, df], ignore_index=True, copy=True)

                #subset to merge_cols + output cols + fractions for wide table
                df = df.loc[:, cols + out_cols + [count_cols[0]+'_fraction', count_cols[1]+'_fraction', 'filter']]

                #drop alts for wide table
                if 'alts' in df.columns:
                    df.drop('alts', axis=1, inplace=True)

                #rename columns
                df.rename(columns = {'filter': sample+'_filter',
                    'maj_hq_count': sample + '_maj',
                    'second_hq_count': sample + '_sec',
                    'maj_hq_count_fraction': sample + '_fraction_maj',
                    'second_hq_count_fraction': sample + '_fraction_sec',
                    'ref_hq_count': sample + '_ref',
                    'nonref_hq_count': sample + '_alt',
                    'ref_hq_count_fraction': sample + '_fraction_ref',
                    'nonref_hq_count_fraction': sample + '_fraction_alt',
                    'snp_ref_hq_count': sample + '_ref',
                    'snp_alt_hq_count': sample + '_alt',
                    'snp_ref_hq_count_fraction': sample + '_fraction_ref',
                    'snp_alt_hq_count_fraction': sample + '_fraction_alt',
                    'number_called_hq': sample + '_total'},
                    inplace = True)

                if agg is None:
                    agg = df
                    agg['total_count'] = agg[sample + '_total'].fillna(0)
                else:
                    agg = agg.merge(df, how='outer', on = [c for c in cols if c != 'alts'], suffixes = ('', '_2'))

                    #add up total count, drop extra column
                    agg['total_count'] = agg['total_count'].fillna(0) + agg[sample + '_total'].fillna(0)
                    #if no probe yet set probes to our probe, check probes, drop extra column
                    agg.loc[agg['probes'].isnull(), 'probes'] = agg.loc[agg['probes'].isnull(), 'probes_2']
                    agg.loc[~agg['probes_2'].isnull() & (agg['probes'] != agg['probes_2']), sample+'_filter'] = 'probes_mismatch'
                    agg.drop('probes_2', axis=1, inplace=True)

                #stats
                log.info('%d/%d: Pileup rows processed, agg shape now %s, agg long shape now %s', ixfile, len(files), str(agg.shape), str(agg_long.shape))
        except pd.io.common.EmptyDataError:
            log.exception("No data found: %s", file)

        coverage_file = file.replace('.pileup.csv', '.targets.csv')
        if os.path.exists(os.path.join(folder, coverage_file)):
            try:
                coverage_df = pd.read_csv(os.path.join(folder, coverage_file), low_memory=False)
                log.info('%d/%d: Found %d coverage rows in %s -> %s', ixfile, len(files), len(coverage_df), coverage_file, sample)

                if len(coverage_df) > 0:
                    coverage_df = coverage_df[ coverage_df['type'] == 'target' ]
                    log.info('%d/%d: %d coverage rows after removing non-targets', ixfile, len(files), len(coverage_df))

                    coverage_df['sample'] = sample
                    coverage_long = pd.concat([coverage_long, coverage_df], ignore_index=True, copy=True)

                    #remove some default columns that we don't need 
                    coverage_df = coverage_df[ [c for c in coverage_df.columns if not c in ['chr', 'start_0', 'end', 'length', 'type']] ]
                    coverage_df.columns = [ c if c == 'id' else sample + '_' + c for c in coverage_df.columns ]

                    if coverage_agg is None:                
                        coverage_agg = coverage_df
                    else:
                        coverage_agg = coverage_agg.merge(coverage_df, how='outer', on = 'id')
            except pd.io.common.EmptyDataError:
                log.exception("No data found: %s", coverage_file)

    if agg is None:
        raise Exception('No data found (in %d files), aborting!' % len(files))

    log.info('Agg Wide: Found %d rows total, %d columns', len(agg), agg.shape[1])
    if snps is None:
        agg = agg[agg['total_count'] > 0]
        log.info('%d rows left after removing zero-counts', len(agg))

    log.info('Agg Long: Found %d rows total, %d columns', len(agg_long), agg_long.shape[1])

    if coverage_agg is not None:
        log.info('Cov Wide: Found %d coverage rows total, %d columns', len(coverage_agg), coverage_agg.shape[1])

    if coverage_long is not None:
        log.info('Cov Long: Found %d coverage rows total, %d columns', len(coverage_long), coverage_long.shape[1])

    agg['chr_num'] = [ 31 if chr[3:] == 'X' else 32 if chr[3:] == 'Y' else 33 if chr[3:] == 'M' else int(chr[3:]) if chr[3:].isdigit() else 99 for chr in agg['chr'] ]
    agg.sort_values(['chr_num', 'chr', 'pos'], inplace=True)
    agg.drop('chr_num', axis=1, inplace=True)

    outname = os.path.join(folder, '%spileups_wide.csv' % ('' if snps is None else 'target_snps_'))
    agg.to_csv(outname, index=False)
    print(outname)

    #figure out the baseline columns
    agg_columns = ['sample']
    if snps is not None:
        agg_columns += snps.columns.tolist()
    else:
        #columns to include in every row
        agg_columns += ['chr', 'pos']

        #add target id, if we have it
        if 'target_id' in df.columns:
            agg_columns.append('target_id')
        if 'target_type' in df.columns:
            agg_columns.append('target_type')

        if ref:
            agg_columns.append('ref') #keep the ref column
            agg_columns.append('alts') #keep the alt column
    agg_columns.append('probes')

    #make two separate lists of columns for the two output files
    agg_long_columns = agg_columns.copy()
    agg_long_detailed_columns = agg_columns.copy()

    #figure out columns for long output
    if snps is not None:
        agg_long_count_cols = ['snp_ref_hq_count', 'snp_alt_hq_count']
    elif ref:
        agg_long_count_cols = ['ref_hq_count', 'nonref_hq_count']
    else:
        agg_long_count_cols = ['maj_hq_count', 'second_hq_count']
    agg_long_columns += agg_long_count_cols
    agg_long_columns.append('number_called_hq')
    agg_long_columns += ['%s_fraction' % count_col for count_col in agg_long_count_cols]
    agg_long_columns.append('filter')

    #figure out columns for detailed output
    agg_long_detailed_columns += ['count_hq_%s' % call for call in call_types]
    agg_long_detailed_columns.append('number_called_hq')
    agg_long_detailed_columns += ['%s_fraction' % count_col for count_col in agg_long_count_cols]
    agg_long_detailed_columns.append('filter')

    #output long table
    outname = os.path.join(folder, '%spileups_long.csv' % ('' if snps is None else 'target_snps_'))
    agg_long[agg_long_columns].to_csv(outname, index=False)
    print(outname)

    #output nucleotide count table
    outname = os.path.join(folder, '%spileups_long_detailed.csv' % ('' if snps is None else 'target_snps_'))
    agg_long[agg_long_detailed_columns].to_csv(outname, index=False)
    print(outname)

    #output coverage, but not if we are running SNPs
    if coverage_agg is not None and snps is None:
        outname = os.path.join(folder, 'target_coverage.csv')
        coverage_agg.to_csv(outname, index=False)
        print(outname)

    #output long coverage
    if coverage_long is not None and snps is None:
        outname = os.path.join(folder, 'target_coverage_long.csv')
        coverage_long.to_csv(outname, index=False)
        print(outname)

    if generate_calls:
        #ideas:
        # http://archive.broadinstitute.org/cancer/cga/mutect
        # http://www.nature.com/nbt/journal/v31/n3/abs/nbt.2514.html
        # http://dkoboldt.github.io/varscan/germline-calling.html

        if snps is not None:
            sys.exit('Cannot call variants when analysing known SNPs')
        if not ref:
            sys.exit('Cannot call variants without reference genome')
            
        def call_variants(pileup_rows):
            #use: 'ref_hq_count', 'nonref_hq_count, alts
            pass

        grouped_pileups = agg_long.groupby(['chr', 'pos', 'ref'])
        variant_calls = grouped_pileups.aggregate(call_variants)
        
        outname = os.path.join(folder, 'variant_calls.csv')
        variant_calls.to_csv(outname, index=False)
        print(outname)

def process_file(input, output, probes_file, snps_file, targets_file, validate_probe_targets,
    fasta_file,
    min_mapq, min_baseq, ignore_groups,
    min_consensus_count,
    min_consensus_fraction = 0.51,
    group_with_mate_positions = False,
    filter_softclipped = True,
    ignore_duplicates = False,
    no_probe_data = False,
    umi_tag_name = None,
    reference_type = 'genome',
    subsample_reads = None,
    debug_umi_groups = False, debug_pos = None, debug = False
):
    if subsample_reads:
        log.info('Subsampling reads to at most %d per position!', subsample_reads)

    #load probes for target validation
    probes_dict = None
    if validate_probe_targets:
        assert not no_probe_data, 'you cannot validate probe targets with --no-probe-data'
        if probes_file is not None and len(probes_file) > 0:
            probes = read_new_probe_design(probes_file)
            probes_dict = probes.to_dict()

    #load SNPs
    snps = None
    snps_dict = None
    if snps_file is not None and len(snps_file) > 0:
        snps = read_snps_txt(snps_file, reference_type = reference_type)
        snps_dict = snps.to_dict()

    #load target regions
    targets = None
    if targets_file is not None and len(targets_file) > 0:
        assert snps is None, 'cannot pileup on target regions and SNPs at the same time!'
        targets = read_targets(targets_file, reference_type = reference_type, file_type = 'bed')
    elif snps is not None:
        #only pileup on SNPs. otherwise SNP in target regions have the wrong target id
        targets = pd.DataFrame({'chr': snps.chr, 'start_0': snps.pos - 1, 'end': snps.pos, 'id': snps.id, 'type': 'snp'})
        targets['start_0'] = targets['start_0'].astype(int)
        targets['end'] = targets['end'].astype(int)
        targets['length'] = targets['end'] - targets['start_0']
        targets.index = range(len(targets))

    #load reference genome
    ref = None
    if fasta_file is not None and len(fasta_file) > 0:
        ref = pyfaidx.Fasta(fasta_file, read_ahead=1000, as_raw=True, sequence_always_upper=True, rebuild=False)
        log.info('Loaded reference genome fasta from %s', fasta_file)

    #figure out actual ends
    #TODO: this is actually no longer in use, since reference_type is forced to 'genome' now.
    if reference_type == 'transcriptome' and targets is not None:
        if ref is None:
            raise Exception('Error: need reference fasta file to do pileup on transcriptome!')

        try:
            targets['end'] = [ len(ref[row.chr]) for row in targets.itertuples() ]
        except KeyError:
            raise
        targets['length'] = targets['end'] - targets['start_0']
        log.info('Fixed target end and lengths')
        print(targets)

    #figure out pileup regions
    pileup_regions = [None]
    if targets is not None:
        #pysam start is 0-based, so we can keep this as is
        pileup_regions = [ (row.chr, row.start_0, row.end, row.id) for row in targets.itertuples() ]

    #no make targets dict (which may actually be SNPs)
    targets_dict = None
    if targets is not None:
        targets_dict = targets.to_dict()

    if len(pileup_regions) > 0 and pileup_regions[0] is not None:
        log.info('Got %d pileup regions: %s ...', len(pileup_regions), str(pileup_regions[0:3]))
    else:
        log.info('Got no specific pileup regions - piling up everywhere!')

    #prepare coverage columns
    if not targets is None:
        targets['coverage_min'] = None
        targets['coverage_avg'] = None
        targets['coverage_zero_percentage'] = None
        targets['called_hq_min'] = None
        targets['called_hq_avg'] = None
        targets['called_hq_zero_percentage'] = None

    log.info('Processing file: %s', input)
    if not os.path.isfile(input):
        sys.exit('input does not exist: %s' % input)

    with pysam.AlignmentFile(input, "rb") as samfile:
        t_shown = t_start = time.time()

        rows = []
        for region_index, region in enumerate(pileup_regions):
            region_length = None
            max_depth = 1e7 if not subsample_reads else subsample_reads

            #our region
            assert region is not None
            region_id = region[3]
            region_length = region[2] - region[1]
            region_coverage = [0] * region_length
            region_coverage_hq = [0] * region_length      

            #first figure out umi groups
            umi_to_group = dict()
            read_metadata = dict()
            if no_probe_data and not umi_tag_name:
                log.info('No UMI info, not grouping reads!')
            elif ignore_groups:
                log.info('Ignore groups is set, not grouping reads!')
            else:
                log.info('Grouping reads in region #%d: %s (length = %d)', region_index, region, region_length)
                iter_grouping = samfile.fetch(region[0], region[1], region[2])
                
                if debug_umi_groups:
                    mate_start_counts = collections.defaultdict(collections.Counter)
                    read_alignment_counter = collections.Counter()

                umi_counts = collections.defaultdict(collections.Counter)
                n_alignments = 0
                for al in iter_grouping:
                    n_alignments += 1
                    raw_read_name = al.query_name

                    #are we dealing with reads that don't have probe data attached?
                    if no_probe_data:
                        original_read_name = raw_read_name
                        read_probe = None
                        read_umi = None

                        #get UMI from tag then!
                        if not umi_tag_name is None:
                            try:
                                read_umi = al.get_tag(umi_tag_name)
                            except KeyError:
                                raise Exception('Encountered alignment with missing UMI tag: "%s", was looking for "%s"' % (raw_read_name, umi_tag_name))
                    else:
                        #parse info from read name (should be bowtie2 compatible)
                        original_read_name, read_probe, read_umi = parse_extended_read_name(raw_read_name)

                    #count this
                    if debug_umi_groups:
                        read_alignment_counter[original_read_name] += 1

                    #make sure we have bytes for read_umi
                    assert read_umi is not None and len(read_umi) > 0
                    read_umi = read_umi.encode('utf-8')

                    #use starts of both mates to uniquely identify alignment
                    mate_starts = get_al_mate_starts(al) if group_with_mate_positions else None

                    #only do this once per read pair and start
                    alignment_tuple = (raw_read_name, mate_starts)
                    read_metadata_tuple = (original_read_name, read_probe, read_umi, [n_alignments])
                    if alignment_tuple in read_metadata:
                        #check metadata (except for n_alignments)
                        for ix in range(3):
                            assert read_metadata[alignment_tuple][ix] == read_metadata_tuple[ix]
                        #add read index
                        read_metadata[alignment_tuple][3].append(n_alignments)
                    else:
                        #record metadata
                        read_metadata[alignment_tuple] = read_metadata_tuple

                        #record umi (needed for grouping)
                        umi_counts[mate_starts][read_umi] += 1
                        if debug_umi_groups:
                            mate_start_counts[read_umi][mate_starts] += 1

                if debug_umi_groups:
                    aligned_reads_per_read_pair = collections.Counter()
                    for v in read_alignment_counter.values():
                        aligned_reads_per_read_pair[v] += 1
                    log.info('Seen alignments per read pair (expect 1/2): %s', str(aligned_reads_per_read_pair))

                    reads_per_alignment = collections.Counter()
                    for v in read_metadata.values():
                        assert len(v[3]) in [1, 2] #we expect to see either 1 or 2 reads for each read_name/mate_positions combo
                        reads_per_alignment[len(v[3])] += 1
                    log.info('Seen reads per alignment (expect 1/2): %s', str(reads_per_alignment))

                    starts_per_umi = collections.Counter()
                    for read_umi, mate_counts_for_umi in mate_start_counts.items():
                        starts_per_umi[len(mate_counts_for_umi)] += 1
                        if len(mate_counts_for_umi) > 1:
                            log.info('Seen %s with multiple mate starts: %s', str(read_umi), str(mate_counts_for_umi.keys()))
                    log.info('Seen mate positions per UMI (expect 1): %s', str(starts_per_umi))

                #free up memory
                log.info('%s: Processed UMIs from %d read pairs with %d alignments.', region_id, len(read_metadata), n_alignments)
                 
                umi_id_offset = 0
                target_raw_umi_count = 0
                if debug_umi_groups:
                    umi_count_histogram = collections.Counter()
                    umi_count_histogram_umis = collections.defaultdict(list)
                #log.info('%s: Got %d different subgroups (starts?) for UMIs.', region_id, len(umi_counts))
                for mate_starts, umi_counts_for_start in umi_counts.items():
                    target_raw_umi_count += len(umi_counts_for_start)

                    if debug_umi_groups:
                        for read_umi, pairs_for_umi in umi_counts_for_start.items():
                            umi_count_histogram[pairs_for_umi] += 1
                            if pairs_for_umi > 5:
                                umi_count_histogram_umis[pairs_for_umi].append(read_umi)

                    #use umi_tools to group UMIs together
                    umi_to_group[mate_starts] = find_umi_groups(umi_counts_for_start, umi_id_offset)
                    umi_id_offset = max(umi_to_group[mate_starts].values()) + 1

                    #log.info('%s: Recorded groups for %d raw UMIs with start = %s, new offset: %d', region_id, len(umi_counts_for_start), str(mate_starts), umi_id_offset)
                log.info('%s: Grouped %d raw UMIs with %d different mate start combinations', region_id, target_raw_umi_count, len(umi_counts))

                if debug_umi_groups:
                    log.info('Raw UMI count hist:')
                    for count in sorted(umi_count_histogram.keys()):
                        log.info('  %d\t%d\t%s', count, umi_count_histogram[count], str(umi_count_histogram_umis[count][:5]) if count > 5 else '')

                #free up memory
                del umi_counts

            #then do actual pileup
            #note: do not filter at this stage since we will do this ourselves
            iter_pileup = samfile.pileup(region[0], region[1], region[2], max_depth = max_depth, stepper = 'nofilter') 
            log.info('Piling up in region #%d: %s (length = %d)', region_index, region, region_length)

            seen_positions = set()
            filtered_pair_qnames = dict()
            for pileup_base in iter_pileup:
                #make sure we don't venture outside the regions, even if we have reads there
                if region is not None:
                    if pileup_base.reference_pos < region[1] or pileup_base.reference_pos >= region[2]:
                        #log.info('Skipping calls at 0-based position %d (outside of target region)', pileup_base.reference_pos)
                        continue

                #if we have a debug pos, ignore everything else!
                if debug_pos is not None:
                    chr_pos = debug_pos.split(':')
                    if pileup_base.reference_name != chr_pos[0] or str(pileup_base.reference_pos + 1) != chr_pos[1]:
                        continue

                #we need to set a max_depth above, but we can make it high
                if not subsample_reads:
                    assert pileup_base.nsegments < max_depth, 'maximum coverage %d reached, increase limit!' % max_depth

                #record row
                row = process_pileup_base(
                    ref,
                    probes_dict,
                    targets_dict,
                    snps_dict,
                    region_index,
                    pileup_base,
                    min_consensus_count,
                    min_consensus_fraction,
                    min_mapq,        
                    min_baseq,              
                    ignore_groups,
                    group_with_mate_positions,
                    validate_probe_targets,
                    filter_softclipped,
                    no_probe_data,
                    read_metadata,
                    filtered_pair_qnames,
                    umi_to_group,
                    debug
                )
                rows.append(row)
                seen_positions.add(pileup_base.reference_pos)

                #record coverage
                if region is not None:
                    coverage_index = pileup_base.reference_pos - region[1]
                    if coverage_index >= 0 and coverage_index < region_length:
                        region_coverage[coverage_index] = row['raw_coverage']
                        region_coverage_hq[coverage_index] = row['number_called_hq']

                t_now = time.time()
                if len(rows) == 1 or len(rows) == 10 or len(rows) % 10000 == 0 or (t_now - t_shown > 60):
                    t_shown = t_now
                    log.info("processed %d positions - %.f sec elapsed, %.4f sec/position, %.1f positions/sec",
                        len(rows),
                        t_now - t_start,
                        (t_now - t_start) / len(rows),
                        len(rows) / (t_now - t_start))
                    #log.info(str(row))

            if region is not None:
                #make sure we have pileup rows for every single nucleotide in the region
                #log.info('%s: Got %d pileup rows so far - filling in missing positions between %d and %d (0-based)', region_id, len(rows), region[1], region[2])
                for pos_0 in range(region[1], region[2]):
                    if not pos_0 in seen_positions:
                        #log.info('Row missing for: %d (0-based)', pos_0)
                        #TODO: could refactor this into just doing get/process in one go, now that we call this anywy
                        row = get_pileup_row(region[0], pos_0,
                            0,
                            targets_dict['id'][region_index] if targets_dict is not None and region_index < len(targets_dict['id']) else None,
                            targets_dict['type'][region_index] if targets_dict is not None and region_index < len(targets_dict['id']) else None,
                            ref,
                            validate_probe_targets = validate_probe_targets)

                        row = process_pileup_row(
                            row,
                            seen_probes = set(),
                            call_groups = dict(),
                            snps_dict = snps_dict,
                            ignore_groups = ignore_groups,
                            min_consensus_count = min_consensus_count,
                            min_consensus_fraction = min_consensus_fraction,
                            min_baseq = min_baseq,
                            ref = ref,
                            debug = debug
                            )
                        rows.append(row)
                log.info('%s: %d pileup rows', region_id, len(rows))

                #calculate overall coverage
                coverage_min = min(region_coverage)
                coverage_avg = 1.0 * sum(region_coverage) / region_length
                coverage_zeros = region_coverage.count(0)
                log.info('%s coverage (raw): min=%d, avg=%f, zeros=%d', region_id, coverage_min, coverage_avg, coverage_zeros)

                coverage_hq_min = min(region_coverage_hq)
                coverage_hq_avg = 1.0 * sum(region_coverage_hq) / region_length
                coverage_hq_zeros = region_coverage_hq.count(0)
                log.info('%s coverage (called + baseq > %d): min=%d, avg=%f, zeros=%d', region_id, min_baseq, coverage_hq_min, coverage_hq_avg, coverage_hq_zeros)

                if targets is not None and region_index < len(targets):
                    targets.loc[region_index, 'coverage_min'] = coverage_min
                    targets.loc[region_index, 'coverage_avg'] = coverage_avg
                    targets.loc[region_index, 'coverage_zero_percentage'] = 100.0 * coverage_zeros / region_length
                    targets.loc[region_index, 'called_hq_min'] = coverage_hq_min
                    targets.loc[region_index, 'called_hq_avg'] = coverage_hq_avg
                    targets.loc[region_index, 'called_hq_zero_percentage'] = 100.0 * coverage_hq_zeros / region_length
                else:
                    log.info('Not recording coverage for region #%d because it is not a target region (SNP?)', region_index)
                
        log.info('Pileups resulted in %d rows', len(rows))
        outname = output + '.pileup.csv'
        if len(rows) > 0:
            df = pd.DataFrame(rows)

            #remove duplicate rows
            dup = df.duplicated(subset = ['chr', 'pos'], keep = 'first')
            if dup.any():
                if ignore_duplicates:
                    log.warn('Dropping %d/%d duplicated rows', dup.sum(), len(df))
                    df = df[~dup]
                    log.info('Now have %d rows after dropping dupes', len(df))
                else:
                    print(df[dup, 'target_id'].value_counts())
                    sys.exit('Found duplicate rows, cancelling! Make sure none of the regions in targets.bed overlap.')
            else:
                log.info('No duplicates found')

            #sort rows (may have become unsorted due to targest)
            df['chr_num'] = [ 31 if chr[3:] == 'X' else 32 if chr[3:] == 'Y' else 33 if chr[3:] == 'M' else int(chr[3:]) if chr[3:].isdigit() else 99 for chr in df['chr'] ]
            df.sort_values(['chr_num', 'chr', 'pos'], inplace=True)
            df.drop('chr_num', axis=1, inplace=True)

            #for call in call_types:
            #    df['percent_%s' % call] = 100.0 * df['count_%s' % call] / df['number_called']
            #for call in call_types:
            #    df['percent_hq_%s' % call] = 100.0 * df['count_hq_%s' % call] / df['number_called_hq']

            df.to_csv(outname, index=False)
        else:
            log.warn('No data to output, creating empty file!')
            open(outname, 'a').close()

        print(outname)

        if snps is None:
            if targets is not None:
                outname = output + '.targets.csv'
                targets.to_csv(outname, index=False)
                print(outname)

def main():
    print('Called with arguments: "%s"' % '" "'.join(sys.argv))
    
    #parse the arguments, which will be available as properties of args (e.g. args.probe)
    parser = argparse.ArgumentParser()
    #specify parameters
    parser.add_argument("--aggregate", help="folder with processed CSV files to aggregate")

    parser.add_argument("--input", help="input bam file name")
    parser.add_argument("--output", help="output file prefix")
    parser.add_argument("--probes", help="CSV file with probes and their arm sequences")
    parser.add_argument("--targets", help="target regions bed file")
    parser.add_argument("--snps", help="target SNPs file")

    parser.add_argument("--ref", help="find reference base at each position (requires --fasta)", action='store_true')
    parser.add_argument("--fasta", help="reference genome fasta (needs faidx)", default='/home/crangen/zding/projects/common/ref/old/hg19_GRCh37/hg19.fa')

    #parser.add_argument("--generate-calls", help="generate variant calls", action='store_true') TODO

    parser.add_argument("--ignore-groups", help="ignore the UMI group entirely", action='store_true')
    parser.add_argument("--group-with-mate-positions", help="only group read pairs if the aligned start positions of both mates are the same", action='store_true')
    parser.add_argument("--min-consensus-count", help="minimum number of consistent read pairs supporting a UMI", default=2, type=int)
    parser.add_argument("--min-consensus-percentage", help="minimum fraction of reads supporting the consensus call in a UMI group", default=51, type=int)

    parser.add_argument("--min-mapq", help="minimum mapping quality (remove read pairs with either read having mapq < X)", default=20, type=int)
    parser.add_argument("--min-baseq", help="minimum base quality (ignore calls from groups below baseq < X at this position)", default=30, type=int)
    parser.add_argument("--ignore-duplicates", help="ignore duplicate positions caused by overlapping targets", action='store_true')
    parser.add_argument("--validate-probe-targets", help="filter base calls outside of the actual probe target regions (eg. if primers have not been trimmed already)", action='store_true')

    parser.add_argument("--no-filter-softclipped", help="do not filter reads that contain softclipped bases", action='store_true')
    parser.add_argument("--no-probe-data", help="do not try to read probe data", action='store_true')
    parser.add_argument("--umi-tag-name", help="name of the BAM tag containing UMIs")

    parser.add_argument("--debug", help="enable debug output", action='store_true')
    parser.add_argument("--debug-pos", help="only do pileup at this position (chr:pos)")
    parser.add_argument("--debug-umi-groups", help="debug umi groups", action='store_true')
    args = parser.parse_args()
    args.generate_calls = None #TODO

    if args.aggregate:
        aggregate(args.aggregate, args.snps, args.ref, args.generate_calls)            
    else:
        assert not (args.ref and args.fasta is None)
        assert args.input is not None
        assert args.output is not None

        process_file(args.input, args.output,
            args.probes, args.snps, args.targets, args.validate_probe_targets,
            args.fasta if args.ref else None,
            args.min_mapq, args.min_baseq,
            ignore_groups = args.ignore_groups,
            group_with_mate_positions = args.group_with_mate_positions,
            min_consensus_count = args.min_consensus_count,
            min_consensus_fraction = args.min_consensus_percentage / 100.0,
            filter_softclipped = not args.no_filter_softclipped,
            ignore_duplicates = args.ignore_duplicates,
            no_probe_data = args.no_probe_data,
            umi_tag_name = args.umi_tag_name,
            debug_umi_groups = args.debug_umi_groups,
            debug_pos = args.debug_pos,
            debug = args.debug)

    return 0