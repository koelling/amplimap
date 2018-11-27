# -*- coding: utf-8 -*-
#python 3 compat
#http://python-future.org/compatible_idioms.html
from __future__ import print_function

import sys
import os
import re
import pprint

from .common import find_umi_groups, parse_extended_read_name

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
import collections
import operator

#use python's argumentparser to support command line parameters like --probe=62
import argparse

#for pileup
import pysam
#for getting ref base 
import pyfaidx
#for debugging umi distance
import distance

#import reader script for reading probe design
from .reader import read_new_probe_design

class ConsensusFilterException(Exception):
    def __init__(self, filter_column):
        self.filter_column = filter_column

def process_file(design, input, output,
        min_mapq, min_consensus_count,
        include_primers, use_naive_groups, ignore_groups,
        ignore_group_mismatch = False,
        debug = False):
    assert design is not None, 'design file missing'
    assert input is not None, 'input file missing'
    assert output is not None, 'output file missing'

    design = read_new_probe_design(design)

    if include_primers:
        probe_column_start = 'probe_start'
        probe_column_end = 'probe_end'
    else:
        probe_column_start = 'target_start'
        probe_column_end = 'target_end'
    log.info('Using columns %s / %s to determine target regions', probe_column_start, probe_column_end)

    log.info('Processing file %s', input)
    if not os.path.isfile(input):
        sys.exit('input does not exist: %s' % input)

    def process_aln_pile(outbam, alns_by_umi):
        #determine counts per umi
        umi_counts = dict( [(umi, len(umi_alns)) for umi, umi_alns in alns_by_umi.items() ]  )
        #group these raw umis into umi groups
        umi_to_group = find_umi_groups( umi_counts )

        group_firstread = {}
        group_qnames = collections.defaultdict(list)
        group_cigars = collections.defaultdict(collections.Counter())
        group_max_mapqual = collections.defaultdict(int)
        group_bases = {}
        group_max_basequal = {}

        for umi, umi_alns in alns_by_umi.items():
            umi_group = umi_to_group[umi]

            for aln in umi_alns:
                #we also need to handle different mate starts
                umi_group_mate = '%d_%d_%d' % (umi_group, aln.next_reference_id, aln.next_reference_start)
                if not umi_group_mate in group_firstread:
                    group_firstread[umi_group_mate] = aln
                    group_bases[umi_group_mate] = collections.defaultdict(dict)
                    group_max_basequal[umi_group_mate] = collections.defaultdict(dict)

                group_qnames[umi_group_mate].append(aln.query_name)
                group_cigars[umi_group_mate][aln.cigarstring] += 1
                group_max_mapqual[umi_group_mate] = max(group_max_mapqual[umi_group_mate], aln.mapping_quality)

                #go through read base by base to form consensus
                for index in range( len(aln.query_sequence) ):
                    group_bases[umi_group_mate][index][aln.query_sequence[index]] += 1
                    group_max_basequal[umi_group_mate][index][aln.query_sequence[index]] = max(
                        group_max_basequal[umi_group_mate][index][aln.query_sequence[index]],
                        aln.query_qualities[index]
                        )

        for umi_group_mate, aln in group_firstread.items():
            aln.query_name = 'group__%s' % umi_group_mate
            aln.set_tag('rc', len(group_qnames[umi_group_mate]))
            aln.set_tag('rn', ';'.join(group_qnames[umi_group_mate]))

            aln.cigarstring = group_cigars[umi_group_mate].most_common(1)[0][0]
            aln.mapping_quality = group_max_mapqual[umi_group_mate]

            #TODO: do we need to do something about the flags?

            #determine most common base at each position
            my_bases = [counts.most_common(1)[0][0] for counts in group_bases[umi_group_mate]]
            my_sequence = ''.join( my_bases )
            #determine highest quality score for most common base at each position
            my_qualities = [group_maxqual[umi_group_mate][base] for base in my_bases]

            #set sequence and qualities
            #note: need to set seq before qual
            aln.query_sequence = my_sequence
            aln.query_qualities = my_qualities

    with tempfile.TemporaryFile() as tmp:
        with pysam.AlignmentFile(input, 'rb') as inbam:
            with pysam.AlignmentFile(tmp, 'wb', template = inbam) as outbam:
                t_shown = t_start = time.time()

                counters = collections.Counter()
                mates_of_invalid_reads = set()
                current_location = (None, None)
                current_alns = collections.defaultdict(list)

                #doc: If until_eof is True, all reads from the current file position will be returned in order as they are within the file. Using this option will also fetch unmapped reads.
                for aln in inbam.fetch(until_eof = True):
                    rindex = 1 if aln.is_read1 else 2 if aln.is_read2 else None
                    assert rindex is not None
                    other_rindex = 3 - rindex

                    #only count rows after this to fulfil assumption later
                    counters['alignments_in'] += 1

                    try:
                        #supplementary alignments mess up the assumption that each mate has one other mate
                        #so we ignore them for now
                        if aln.is_supplementary:
                            raise ConsensusFilterException('is_supplementary')

                        #we also don't want unmapped reads here
                        if aln.is_unmapped:
                            raise ConsensusFilterException('is_unmapped')

                        #and no unmapped mates
                        if aln.mate_is_unmapped:
                            raise ConsensusFilterException('mate_is_unmapped')

                        if debug:
                            print(aln)

                        #parse info from read name (should be bowtie2 compatible)
                        read_name, read_probe, read_umi = parse_extended_read_name(aln.query_name)

                        #make sure we have a valid read pair
                        is_valid_pair = True
                        is_valid_pair = is_valid_pair and aln.is_reverse != aln.mate_is_reverse #should be opposite
                        is_valid_pair = is_valid_pair and aln.reference_id == aln.next_reference_id #should be on same chr
                        if debug:
                            print('is_valid_pair =', is_valid_pair)

                        #and no unmapped mates
                        if not is_valid_pair:
                            raise ConsensusFilterException('not_valid_pair')

                        #check probe data
                        probe_data = design.loc[pname, :]
                        matches_probe = aln.reference_name == probe_data['chr']
                        if matches_probe:
                            #TODO: this may be the wrong way round!
                            my_column = probe_column_start
                            if aln.is_reverse or (not aln.is_reverse and aln.is_read2):
                                matches_probe = matches_probe and aln.reference_end + 1 == probe_data[probe_column_end]
                            else:
                                matches_probe = matches_probe and aln.reference_start + 1 == probe_data[probe_column_start]
                        if debug:
                            print('matches_probe =', matches_probe)

                        #and no unmapped mates
                        if not matches_probe:
                            raise ConsensusFilterException('mismatches_probe')

                        #update flags if this is a mate of a read that has been filtered
                        my_tuple = (aln.query_name, aln.reference_id, aln.reference_start)
                        if my_tuple in mates_of_invalid_reads:
                            mates_of_invalid_reads.remove( my_tuple )
                            raise ConsensusFilterException('mate_filtered')

                        #is this a new pile of alignments?
                        if current_location[0] != aln.reference_id or current_location[1] != aln.reference_start:
                            process_aln_pile(outbam, current_alns)

                            current_alns = collections.defaultdict(list)
                            current_location = (aln.reference_id, aln.reference_start)

                        assert read_umi is not None and len(read_umi) > 0
                        read_umi_bytes = read_umi.encode('utf-8')
                        current_alns[read_umi_bytes].append(aln)
                    except ConsensusFilterException as ex:
                        #if the row is filtered, we raise an exception with message being the corresponding counter column
                        counters[ex.filter_column] += 1

                        #keep track of mates we still need to filter out
                        if ex.filter_column != 'mate_filtered':
                            mates_of_invalid_reads.add( (aln.query_name, aln.next_reference_id, aln.next_reference_start) )

                    #print some stats
                    t_now = time.time()
                    if n_rows > 0 and (n_rows % 10000 == 0 or (t_now - t_shown > 60)):
                        t_shown = t_now
                        log.info("processed %d alignment rows - %.f sec elapsed, %.4f sec/pair, %.1f row/sec",
                            counters['alignments_in'],
                            t_now - t_start,
                            (t_now - t_start) / counters['alignments_in'],
                            counters['alignments_in'] / (t_now - t_start))

                #run this one final time to catch the last pile
                process_aln_pile(outbam, current_alns)

            #summary stats
            log.info('processed %d alignment rows (reads) total.', counters['alignments_in'])

        log.info('Re-processing file to remove invalid pairs...', counters['alignments'])
        with pysam.AlignmentFile(tmp, 'rb') as inbam:
            with pysam.AlignmentFile(output, 'wb', template = inbam) as outbam:
                for aln in inbam.fetch(until_eof = True):
                    counters['alignments_first_pass'] += 1

                    my_tuple = (aln.query_name, aln.reference_id, aln.reference_start)
                    if my_tuple in mates_of_invalid_reads:
                        mates_of_invalid_reads.remove( my_tuple )
                        counters['mate_filtered_second_pass'] += 1
                    else:
                        counters['alignments_second_pass'] += 1
                        outbam.write(aln)

        log.info('read %d rows from tmp file, wrote %d rows.', counters['alignments_first_pass'], counters['alignments_second_pass'])
        print(output)

def main():
    log.info('Called with arguments: "%s"', '" "'.join(sys.argv))
    
    #parse the arguments, which will be available as properties of args (e.g. args.probe)
    parser = argparse.ArgumentParser()
    #specify parameters
    parser.add_argument("--aggregate", help="folder with processed CSV files to aggregate")
    parser.add_argument("-d", "--design", help="CSV file with probes, arms sequences and locations")
    parser.add_argument("-i", "--input", help="input bam file name")
    parser.add_argument("-o", "--output", help="output bam file name")
    parser.add_argument("--use-naive-groups", help="use UMIs from um tag instead of UG tags", action='store_true')
    parser.add_argument("--ignore-groups", help="ignore the UMI group entirely", action='store_true')
    parser.add_argument("--ignore-group-mismatch", help="ignore group mismatches (to deal with old alignments that have separate UMIs)", action='store_true')
    parser.add_argument("--min-consensus-count", help="minimum number of consistent read pairs supporting a UMI", default=2, type=int)
    parser.add_argument("--min-mapq", help="minimum mapping quality (remove read alignments with either read having mapq < X)", default=20, type=int)
    parser.add_argument("--include-primers", help="include primers when determining on/off target", action='store_true')

    parser.add_argument("--debug", help="enable debug output", action='store_true')
    args = parser.parse_args()

    process_file(args.design, args.input, args.output,
        args.min_mapq, args.min_consensus_count,
        args.include_primers, args.use_naive_groups, args.ignore_groups,
        args.ignore_group_mismatch,
        debug = args.debug)

    return 0

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    sys.exit(main())