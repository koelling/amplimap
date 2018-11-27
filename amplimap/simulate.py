# -*- coding: utf-8 -*-
"""
``amplimap.simulate`` contains methods for simulating variants in reads by replacing a specific DNA sequence with a different sequence,
some percentage of the time.
"""

import random
import gzip
import Bio.Seq
import collections
import re

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
log.addHandler(sh)

def check_parameters(wildcards):
    assert re.match(r"^[ACGT]{10,}", wildcards['search']), '\nplease provide at least 10 characters of ACGT to search for\n'
    assert re.match(r"^[ACGT]{1,}", wildcards['replace']), '\nplease provide at least 1 character of ACGT to replace\n'
    assert float(wildcards['percentage']) > 0.0, '\nreplace_percentage must be larger than 0\n'
    assert float(wildcards['percentage']) <= 100.0, '\nreplace_percentage must be smaller or equal to 100\n'

def make_simulated_reads(input, output, wildcards, config):
    #figure out search/replace strings both on fwd and reverse strands
    search_fwd = wildcards['search']
    search_rev = str(Bio.Seq.Seq(search_fwd).reverse_complement())
    replace_fwd = wildcards['replace']
    replace_rev = str(Bio.Seq.Seq(replace_fwd).reverse_complement())

    #we need a bytes-like object here, not strings
    search_fwd = search_fwd.encode("utf-8")
    search_rev = search_rev.encode("utf-8")
    replace_fwd = replace_fwd.encode("utf-8")
    replace_rev = replace_rev.encode("utf-8")

    #counters for statistics
    counts = collections.Counter()
    count_cols = ['occurrences_fwd_r1', 'occurrences_fwd_r2', 'occurrences_rev_r1', 'occurrences_rev_r2', 'replacements_fwd_r1', 'replacements_fwd_r2', 'replacements_rev_r1', 'replacements_rev_r2', 'total_umi_occurrences', 'total_umi_replacements']

    #remember for each umi whether to replace in the read or not
    replace_by_umi = {}

    with gzip.open(input[0], 'r') as file_in0, gzip.open(input[1], 'r') as file_in1, \
    gzip.open(output[0], 'w') as file_out0, gzip.open(output[1], 'w') as file_out1, \
    open(output['stats'], 'w') as stats:
        #write header line to stats file
        stats.write('sample,%s\n' % ','.join(count_cols))

        #keep track of umis
        occurrence_umis = set()
        replacement_umis = set()

        #loop through input lines (we need map(list, x) here to get lists we can assign values to instead of tuples)
        for ix, read_lines in enumerate(map(list, zip(file_in0, file_in1))):
            if ix % 4 == 1:
                #decide to replace for both mates -- but cache results for umi so we do the same thing for every member of the consensus group
                my_umi = None
                if config['parse_reads']['umi_one'] > 0 or config['parse_reads']['umi_two'] > 0:
                    #extract UMI - note that this does not do any of the smarter UMI grouping we do later on
                    my_umi = read_lines[0][:config['parse_reads']['umi_one']] + b'_' + read_lines[1][:config['parse_reads']['umi_two']] #not the b'_' here -- need to use bytes instead of str
                    assert len(my_umi) == config['parse_reads']['umi_one'] + 1 + config['parse_reads']['umi_two']

                    if not my_umi in replace_by_umi:
                        replace_by_umi[my_umi] = random.random() <= float(wildcards['percentage']) / 100.0
                    do_replace = replace_by_umi[my_umi]
                else:
                    do_replace = random.random() <= float(wildcards['percentage']) / 100.0

                #loop over both mates in the pair
                mates_found = 0
                for readix in range(2):
                    #first try to find fwd, then reverse -- this is slightly weird, but prevents double replacement
                    if search_fwd in read_lines[readix]:
                        if my_umi is not None:
                            occurrence_umis.add(my_umi)
                        mates_found += 1
                        counts['occurrences_fwd_r%d' % (readix+1)] += 1
                        if do_replace:
                            read_lines[readix] = read_lines[readix].replace(search_fwd, replace_fwd, 1)
                            counts['replacements_fwd_r%d' % (readix+1)] += 1
                            if my_umi is not None:
                                replacement_umis.add(my_umi)
                    elif search_rev in read_lines[readix]:
                        if my_umi is not None:
                            occurrence_umis.add(my_umi)
                        mates_found += 1
                        counts['occurrences_rev_r%d' % (readix+1)] += 1
                        if do_replace:
                            read_lines[readix] = read_lines[readix].replace(search_rev, replace_rev, 1)
                            counts['replacements_rev_r%d' % (readix+1)] += 1
                            if my_umi is not None:
                                replacement_umis.add(my_umi)

                if mates_found == 1:
                    counts['only_one_mate'] += 1                   
            
            file_out0.write(read_lines[0])
            file_out1.write(read_lines[1])

        counts['total_umi_occurrences'] = len(occurrence_umis)
        counts['total_umi_replacements'] = len(replacement_umis)

        stats.write('%s,%s\n' % (wildcards['sample_with_lane'], ','.join([str(counts[col]) for col in count_cols])))  

def stats_replacements_agg(input, output):
    import pandas as pd

    merged = None
    for file in input:
        print('Reading', file, '...')
        try:
            df = pd.read_csv(file, index_col = False)
            print('Data shape:', str(df.shape))

            if merged is None:
                merged = df
            else:
                merged = merged.append(df, ignore_index = True)
        except pd.io.common.EmptyDataError:
            print('No data for', file, ', skipping.')

    print('Merged data shape:', str(merged.shape))

    for c in ['occurrences', 'replacements']:
        merged['total_%s' % c] = (merged['%s_fwd_r1' % c] + merged['%s_rev_r1' % c] + merged['%s_fwd_r2' % c] + merged['%s_rev_r2' % c])

    merged['replaced_percentage_reads'] = 100.0 * merged['total_replacements'] / merged['total_occurrences']
    merged['replaced_percentage_umis'] = 100.0 * merged['total_umi_replacements'] / merged['total_umi_occurrences']

    if merged['total_replacements'].sum() <= 0:
        merged.to_csv(output[0]+'.error.csv', index = False)
        print('Wrote replacement stats to %s' % output[0]+'.error.csv')
    assert merged['total_replacements'].sum() > 0, 'ERROR: no replacements made!'

    merged.to_csv(output[0], index = False)

    