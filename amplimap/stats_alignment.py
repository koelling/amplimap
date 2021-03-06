# -*- coding: utf-8 -*-
"""
This module contains methods for generating statistics about the alignments generated by amplimap.
This includes the number of reads and read families per probe, per sample as well as data on off-target
alignments.
"""

# python 3 compat
# http://python-future.org/compatible_idioms.html
from __future__ import print_function

import sys
import os
import re
import pprint

from .common import parse_extended_read_name, find_umi_groups

"""
NOTES:
The BAM contains alignments of reads pairs. We can have these situations:
- Multiple possible alignments with different properties, each can be one of:
    - Alignment of both mates (both have is_unmapped & mate_is_unmapped set to False)
    - Alignment of one mate, other unmapped (one has is_unmapped set, one has mate_is_unmapped set)

- Or: Both unmapped (both have is_unmapped and mate_is_unmapped set) - if this is the case, it should be the only entry for this read pair

APPROACH:
- The BAM is coordinate sorted, so we need to read through it and cache, keeping track of mate positions to find the right mates for each pair.
- For each new mate, we can either:
    - If self is mapped:
        - If mate is mapped:
            - Already have the mate -> process and remove mate from cache
            - Not have the mate -> add to cache and wait for mate
            - NB TODO: actually this will fail with bowtie2 because one read alignment can have multiple mate alignments.
                - should keep all reads in cache, ideally forever or at least per chr in case there are more alignments
                - at the moment we are just ignoring these and outputting these to the log as orphaned reads at the end
        - If not:
            - Count as alignment with unmapped mate [but maybe only if no better alignment]
    - If not:
        - If mate1, count unmapped pair, else ignore (should never see 2x)
- By the end of the file we should have resolved all mates. Anything left is suspicous!

But, question how we deal with multiple alignments per pair. Actually we are only interested in the top aligment, no need to record every single bad secondary alignment.

"""

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
log.addHandler(sh)
import time

# use pandas to read the CSV file and write output files
import pandas as pd
# for defaultdict + sorting
import collections
import operator

# use python's argumentparser to support command line parameters like --probe=62
import argparse

# for pileup
import pysam
# for getting ref base
import pyfaidx
# for debugging umi distance
import distance

# import reader script for reading probe design
from .reader import read_new_probe_design

def aggregate(folder):
    log.info('Aggregating tables in %s', folder)

    agg = None
    files = [file for file in os.listdir(folder) if file.endswith('.stats_alignment.csv')]
    files.sort()
    for file in files:
        sample = file
        sample = sample.replace('.stats_alignment.csv', '')
        sample = sample.replace('_L001_R1_001', '')
        sample = sample.replace('__MIP_TRIMMED_', '')

        try:
            df = pd.read_csv(os.path.join(folder, file))
            log.info('Found %d rows in %s -> %s', len(df), file, sample)

            df['sample'] = sample

            if agg is None:
                agg = df
            else:
                agg = agg.append(df)
        except pd.errors.EmptyDataError:
            log.exception("No data found: %s", file)

    assert agg is not None and len(agg) > 0, \
        '\n\nABORTED: Did not find any valid read alignments for any sample. Please check configuration, probe design table and samples!\n\n'

    agg.sort_values(['sample', 'probe'], inplace=True)
    agg = agg[ ['sample', 'probe'] + [c for c in agg.columns if not c in ['sample', 'probe']] ]

    log.info('Found %d rows total, %d columns', len(agg), agg.shape[1])

    outname = os.path.join(folder, 'stats_alignment.csv')
    agg.to_csv(outname, index=False)
    print(outname)

def process_file(
    probes_path: str,
    input_path: str,
    output_path: str,
    min_mapq: int,
    min_consensus_count: int,
    include_primers: bool,
    use_naive_groups: bool,
    ignore_groups: bool,
    ignore_group_mismatch: bool = False,
    debug: bool = False,
    targets_path: str = None,
):
    """
    Loop through coordinate-sorted BAM file and collect alignment statistics.

    Args:
        probes_path (str): Path to probes design CSV file
        input_path (str): Path to input BAM file
        output_path (str): Prefix of output CSV file ('.stats_alignment.csv' will be appended)

    Output file columns:

    - read_pairs_total: Total number of read pairs (one pair may have multiple alignments)
    - alignments_total: Total number of read pair alignments (one pair may have multiple alignments)
    - alignments_good: Total number of alignments that were on target (coordinates match target coordinates) and fully covered the target region
    - umis_good: Total number of UMI groups (read families) that were on target and fully covered the target region
    - umis_good_coverage_ge_NN: Total number of UMI groups (read families) that were on target, fully covered the target region and also had at least NN supporting read pairs
    - alignments_partial: Total number of alignments that were on target but did not fully cover the target region
    - alignments_off_target: Total number of alignments with coordinates different from the target coordinates
    - alignments_unmapped_single: Total number of alignments where one mate was unmapped
    - alignments_unmapped_both: Total number of alignments where both mates were unmapped
    - alignments_flagged: Total number of alignments where at least one mate was flagged as QC fail or supplementary alignment
    - alignments_invalid_pair: Total number of alignments where mates were on different chromosomes or not in the expected orientation (one forward, one reverse)
    - alignments_mapq_lt_NN: Total number of alignments with at least one mate with mapping quality under NN
    - multimapping_alignments: Total number of read pairs that had multiple possible alignments
    - offtarget_locations: Locations and counts for the three most common off-target alignment locations

    Will only process a read pair once both alignments have been found, which means that
    the first mate has to be kept in memory until the second mate has been reached.
    This may lead to high memory usage if mates are far away and there are many alignments
    but seems to work fine in practice.

    However, reads with multiple alternative mate alignments may not be handled
    properly since an alignment is only kept in memory until the first matching mate has been found.
    These will be reported as "orphaned mates" at the end of the run.
    """
    assert input_path is not None, 'input file path missing'
    assert output_path is not None, 'output file path missing'

    # if we did not get a design path we can only look at on/off target
    if probes_path is None or len(probes_path) == 0:
        assert target_path is not None, 'target file path is missing'
        # TODO

    design = read_new_probe_design(probes_path)

    if include_primers:
        probe_column_start = 'probe_start'
        probe_column_end = 'probe_end'
    else:
        probe_column_start = 'target_start'
        probe_column_end = 'target_end'
    log.info('Using columns %s / %s to determine target regions', probe_column_start, probe_column_end)

    log.info('Processing file %s', input_path)
    if not os.path.isfile(input_path):
        sys.exit('input does not exist: %s' % input_path)

    with pysam.AlignmentFile(input_path, "rb") as samfile:
        t_shown = t_start = time.time()

        reads = collections.defaultdict(dict)
        reads_processed = set()
        reads_umi_checked = set()

        n_rows = 0
        n_pairs = collections.defaultdict(int)
        n_alignments = collections.defaultdict(int)
        n_alignments_good = collections.defaultdict(int)
        n_alignments_partial = collections.defaultdict(int)
        n_alignments_off_target = collections.defaultdict(int)
        n_alignments_unmapped_single = collections.defaultdict(int)
        n_alignments_unmapped_both = collections.defaultdict(int)
        n_alignments_flagged = collections.defaultdict(int)
        n_alignments_invalid = collections.defaultdict(int)
        n_alignments_low_mapq = collections.defaultdict(int)
        n_alignments_multimappers = collections.defaultdict(int)

        probe_umi_alignments = collections.defaultdict(collections.Counter)
        n_umis_good = collections.defaultdict(int)
        n_umis_good_coverage_ge_MIN = collections.defaultdict(int)
        probe_offtarget_locations = collections.defaultdict(collections.Counter)

        # make sure we have zeros for each probe id
        for pname in design['id']:
            for counter in [
                n_pairs,
                n_alignments,
                n_alignments_good,
                n_alignments_partial,
                n_alignments_off_target,
                n_alignments_unmapped_single,
                n_alignments_unmapped_both,
                n_alignments_flagged,
                n_alignments_invalid,
                n_alignments_low_mapq,
                n_alignments_multimappers,

                n_umis_good,
                n_umis_good_coverage_ge_MIN,
            ]:
                counter[pname] = 0

        iter = samfile.fetch(until_eof = True)
        for x in iter:
            qname = x.query_name
            assert qname is not None
            rindex = 1 if x.is_read1 else 2 if x.is_read2 else None
            assert rindex is not None
            other_rindex = 3 - rindex

            # supplementary alignments mess up the assumption that each mate has one other mate
            # so we ignore them for now
            if x.is_supplementary:
                # log.warn('Ignorning supplementary alignment for %s', x.query_name)
                continue

            # only count rows after this to fulfil assumption later
            n_rows += 1

            # parse info from read name (should be bowtie2 compatible)
            read_name, read_probe, read_umi = parse_extended_read_name(qname)

            rdata = {}
            rdata['index'] = rindex
            rdata['mapq'] = x.mapping_quality
            rdata['has_bad_flags'] = (x.is_qcfail) or \
                (x.is_supplementary)
                # (x.is_secondary) or \ #TODO: removed this to allow secondary hits from bowtie2, not sure if good idea or not. this column is often (always?) zero for bwa runs though...
            rdata['unmapped'] = x.is_unmapped
            if x.is_unmapped:
                rdata['chr'] = None
                rdata['start'] = None
                rdata['end'] = None
                rdata['reverse'] = None
            else:
                rdata['chr'] = x.reference_name
                rdata['start'] = x.reference_start + 1
                rdata['end'] = x.reference_end
                rdata['reverse'] = x.is_reverse
            rdata['probe'] = read_probe
            rdata['umi'] = read_umi

            pname = read_probe
            if debug:
                print('>', qname, '-', rindex, '(', rdata['probe'], '/', rdata['chr'], ':', rdata['start'], '-', rdata['end'], ')', 'mate mapped =', not x.mate_is_unmapped)

            # keep track of whether we processed this alignment pair now
            unprocessed = False

            # are we unmapped
            if x.is_unmapped:
                if debug:
                    print('is unmapped -- mate_is_unmapped =', x.mate_is_unmapped, ', read1 =', x.is_read1)

                # is our mate also unmapped, and this is read1
                if x.mate_is_unmapped and x.is_read1:
                    n_alignments_unmapped_both[pname] += 1
                else:
                    # if this is read2 we already counted read1 (or will count)
                    # if mate is not unmapped we will count this when we read the mate
                    unprocessed = True
            # are we mapped but our mate is unmapped?
            elif x.mate_is_unmapped:
                if debug:
                    print('is mapped but mate_is_unmapped ==', x.mate_is_unmapped)

                n_alignments_unmapped_single[pname] += 1
            # mapped pair
            else:
                # figure out where the mate should be
                mate_reference_id = x.next_reference_id
                mate_start = x.next_reference_start + 1 #one-indexed to agree with above

                # did we already see the paired mate?
                paired_mate_number = None
                if other_rindex in reads[qname]:
                    for other_mate_number, other_mate_alignment in enumerate(reads[qname][other_rindex]):
                        # is this the correct alignment that we are mated with?
                        if (mate_reference_id == 0 or mate_reference_id == other_mate_alignment['reference_id']) and mate_start == other_mate_alignment['start']:
                            paired_mate_number = other_mate_number
                            break

                if debug:
                    print('mate @ #', mate_reference_id, ':', mate_start, ' -> paired_mate_number =', paired_mate_number)

                if paired_mate_number is None:
                    # for debugging
                    rdata['qname'] = qname
                    rdata['mate_reference_id'] = mate_reference_id
                    rdata['mate_start'] = mate_start

                    # we will be matching using the id, not the chr name
                    rdata['reference_id'] = x.reference_id

                    # keep this read around until we find the mate
                    if not rindex in reads[qname]:
                        reads[qname][rindex] = list()
                    reads[qname][rindex].append(rdata)
                    unprocessed = True
                else:
                    # process the pair
                    rd = {}
                    rd[rindex] = rdata
                    rd[other_rindex] = reads[qname][other_rindex][paired_mate_number]

                    if debug:
                        print('Processing pair:')
                        print(rd)

                    # make sure this is a good pair
                    is_valid_pair = True
                    is_valid_pair = is_valid_pair and rd[1]['reverse'] != rd[2]['reverse']
                    is_valid_pair = is_valid_pair and rd[1]['chr'] == rd[2]['chr']
                    if debug:
                        print('is_valid_pair =', is_valid_pair)

                    # sanity checks
                    if rd[1]['probe'] != rd[2]['probe']:
                        log.error('Probe mismatch for %s', qname)
                        print(rd)
                        assert False

                    # this should not happen anymore, since we catch these abvoe
                    assert not (rd[1]['unmapped'] or rd[2]['unmapped'])

                    # check probe data
                    probe_data = design.loc[pname, :]
                    covers_entire_probe = False
                    matches_probe = rd[1]['chr'] == probe_data['chr']
                    if matches_probe:
                        if rd[1]['reverse']:
                            matches_probe = matches_probe and rd[2]['start'] == probe_data[probe_column_start]
                            matches_probe = matches_probe and rd[1]['end'] == probe_data[probe_column_end]
                            covers_entire_probe = rd[2]['end'] >= rd[1]['start']
                        else:
                            matches_probe = matches_probe and rd[1]['start'] == probe_data[probe_column_start]
                            matches_probe = matches_probe and rd[2]['end'] == probe_data[probe_column_end]
                            covers_entire_probe = rd[1]['end'] >= rd[2]['start']
                    if debug:
                        print('matches_probe =', matches_probe)

                    # decide which bin to put this in
                    if not is_valid_pair:
                        n_alignments_invalid[pname] += 1
                    elif rd[1]['has_bad_flags'] or rd[2]['has_bad_flags']:
                        n_alignments_flagged[pname] += 1
                    elif rd[1]['mapq'] < min_mapq or rd[2]['mapq'] < min_mapq:
                        n_alignments_low_mapq[pname] += 1
                    elif matches_probe:
                        if not covers_entire_probe:
                            n_alignments_partial[pname] += 1
                        else:
                            n_alignments_good[pname] += 1

                            assert rd[1]['umi'] == rd[2]['umi']
                            probe_umi_alignments[pname][ rd[1]['umi'].encode('utf-8') ] += 1
                    else:
                        # otherwise we had a good mapping, but not in the expected location
                        n_alignments_off_target[pname] += 1

                        start_name = '%s:%d-%d' % (rd[1]['chr'], rd[1]['start'], rd[2]['end'])
                        if rd[1]['reverse']:
                            start_name = '%s:%d-%d' % (rd[1]['chr'], rd[2]['start'], rd[1]['end'])
                        probe_offtarget_locations[pname][start_name] += 1

                    if debug and len(reads_processed) > 5:
                        log.warn('Breaking for debug')
                        break

                    # remove paired mate from memory
                    del reads[qname][other_rindex][paired_mate_number]
                    # also remove the containing lists
                    if len(reads[qname][other_rindex]) == 0:
                        del reads[qname][other_rindex]
                    if len(reads[qname]) == 0:
                        del reads[qname]

            # keep track of processed read pairs
            if not unprocessed:
                n_alignments[pname] += 1

                # remember we processed this one
                if qname in reads_processed:
                    n_alignments_multimappers[pname] += 1
                else:
                    n_pairs[pname] += 1
                    reads_processed.add(qname)

            # print some stats
            t_now = time.time()
            if n_rows > 0 and (n_rows % 10000 == 0 or (t_now - t_shown > 60)):
                t_shown = t_now
                log.info("processed %d rows - %.f sec elapsed, %.4f sec/pair, %.1f row/sec",
                    n_rows,
                    t_now - t_start,
                    (t_now - t_start) / n_rows,
                    n_rows / (t_now - t_start))

        # summary stats
        log.info('processed %d alignment rows total.', n_rows)
        log.info('processed %d reads total.', len(reads_processed))

        # check if we have any reads left over?!
        orphan_reads = [r for rs in reads.values() for r in rs.values()]
        if len(orphan_reads) > 0:
            # TODO: this may happen if a single aligned read is paired with multiple alignments of its mate.
            # we will remember the first mate, process the pair once we get the second mate and then forget both
            # next time we see a new alignment for one of the mates we can never find the other alignment, since we have already deleted it
            # I don't think there is an easy way to resolve this unfortunately.
            # Options would be to keep _everything_ in memory for each contig (and lose a few), or to use name-sorted BAM files,
            # or to see if bowtie2 gives us some idea of how many mates there are for a given read, so we know which ones to keep around for how long (reference counting)?
            log.info('Found %d orphaned read mates. First 10:', len(orphan_reads))
            pprint.pprint(orphan_reads[0:10])

            # most of these will be on different references (in our test case all are)
            orphan_reads_same_ref = [orphan_read for orphan_read_list in orphan_reads for orphan_read in orphan_read_list if orphan_read['mate_reference_id'] == orphan_read['reference_id']]
            if len(orphan_reads_same_ref) > 0:
                log.warn('Of these, %d orphaned read mates with same ref id. First 10:', len(orphan_reads_same_ref))
                pprint.pprint(orphan_reads_same_ref[0:10])
            else:
                log.info('However, all of these involve pairs crossing chromosomes.')

        # at the end, process the umi stats
        for (pname, umi_raw_counts) in probe_umi_alignments.items():
            # use umi_tools to group raw UMIs together
            umi_to_group = find_umi_groups(umi_raw_counts)

            # import pprint
            # umi_group_debug = collections.defaultdict(list)
            # for umi, group in umi_to_group.items():
            #     umi_group_debug[group].append('%s=%d' % (umi, umi_raw_counts[umi]))
            # print(pname)
            # pprint.pprint(umi_group_debug)

            # now do a new count
            umi_group_counts = collections.Counter()
            for umi, group in umi_to_group.items():
                # we want the number of reads supporting the umi group here, which is
                # the sum of the numbers of reads supporting each grouped umi
                umi_group_counts[group] += umi_raw_counts[umi]

            # now use these counts
            for _, count in umi_group_counts.items():
                n_umis_good[pname] += 1
                if count >= min_consensus_count:
                    n_umis_good_coverage_ge_MIN[pname] += 1

        # and find the top offtarget locations for each probe
        top_offtarget_locations = dict(
            (
                ( pname, ';'.join(['%s=%d' % (location, count) for location, count in otl_counter.most_common(3)]) ) for pname, otl_counter in probe_offtarget_locations.items()
            )
        )

        # make a table
        df = pd.DataFrame(collections.OrderedDict((
            ('read_pairs_total', n_pairs),
            ('alignments_total', n_alignments),
            ('alignments_good', n_alignments_good),
            ('umis_good', n_umis_good),
            ('umis_good_coverage_ge_%d' % min_consensus_count, n_umis_good_coverage_ge_MIN),
            ('alignments_partial', n_alignments_partial),
            ('alignments_off_target', n_alignments_off_target),
            ('alignments_unmapped_single', n_alignments_unmapped_single),
            ('alignments_unmapped_both', n_alignments_unmapped_both),
            ('alignments_flagged', n_alignments_flagged),
            ('alignments_invalid_pair', n_alignments_invalid),
            ('alignments_mapq_lt_%d' % min_mapq, n_alignments_low_mapq),
            ('multimapping_alignments', n_alignments_multimappers),
            ('offtarget_locations', top_offtarget_locations)
            )))

        # fill with 0s and make ints
        df = pd.concat([df.loc[:, [c for c in df.columns if c != 'offtarget_locations']].fillna(0).astype(int), df['offtarget_locations']],
            axis=1)
        # set and sort by probe name (row index)
        df.index.rename('probe', inplace=True)
        df.sort_index(0, inplace=True)

        # ASSERTIONS
        # read_pairs_total should be == read_pairs_total from stats_overview
        assert df['read_pairs_total'].sum() == len(reads_processed) #each read pair should be in one probe row
        # assert df['alignments_total'].sum()*2 == n_rows #each alignment has two rows -- this falls apart if one mate is outside the target region!
        assert (df['read_pairs_total'] == df['alignments_total'] - df['multimapping_alignments']).all() #should have one non-multimapping alignment per pair
        assert (df['alignments_total'] == df[ [c for c in df.columns if c.startswith('alignments_') and c != 'alignments_total'] ].sum(axis=1)).all() #each alignment should be counted once in the alignment_XX columns

        df['fraction_pairs_good'] = 1.0 * df['alignments_good'] / df['read_pairs_total']

        if debug:
            print(df)
        else:
            outname = output_path
            df.to_csv(outname, index=True)
            print(outname)

def main():
    log.info('Called with arguments: "%s"', '" "'.join(sys.argv))

    # parse the arguments, which will be available as properties of args (e.g. args.probe)
    parser = argparse.ArgumentParser()
    # specify parameters
    parser.add_argument("--aggregate", help="folder with processed CSV files to aggregate")
    parser.add_argument("-d", "--design", help="CSV file with probes, arms sequences and locations")
    parser.add_argument("-i", "--input", help="input bam file name")
    parser.add_argument("-o", "--output", help="output file prefix + '.stats_alignment.csv'")
    parser.add_argument("--use-naive-groups", help="use UMIs from um tag instead of UG tags", action='store_true')
    parser.add_argument("--ignore-groups", help="ignore the UMI group entirely", action='store_true')
    parser.add_argument("--ignore-group-mismatch", help="ignore group mismatches (to deal with old alignments that have separate UMIs)", action='store_true')
    parser.add_argument("--min-consensus-count", help="minimum number of consistent read pairs supporting a UMI", default=2, type=int)
    parser.add_argument("--min-mapq", help="minimum mapping quality (remove read alignments with either read having mapq < X)", default=20, type=int)
    parser.add_argument("--include-primers", help="include primers when determining on/off target", action='store_true')

    parser.add_argument("--debug", help="enable debug output", action='store_true')
    args = parser.parse_args()

    if args.aggregate:
        aggregate(args.aggregate)
    else:
        process_file(args.design, args.input, args.output + '.stats_alignment.csv',
            args.min_mapq, args.min_consensus_count,
            args.include_primers, args.use_naive_groups, args.ignore_groups,
            args.ignore_group_mismatch,
            debug = args.debug)

    return 0

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    sys.exit(main())