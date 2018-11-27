#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module provides the amplimap.merge_folders.main() function called by the ``amplimap_merge`` script.
This script merges coverage data and variant calls from different working directories together,
making it possible to merge samples sequenced in different runs into a single output file.
"""

#python 3 compat
#http://python-future.org/compatible_idioms.html
from __future__ import print_function

import sys
import os
import re

import logging
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

import time

#use pandas to read the CSV file and write output files
import pandas as pd
#for defaultdict + sorting
import collections, itertools
import operator

import argparse

files_to_merge = ['variants_raw/variants_summary.csv', 'variants_raw/variants_summary_filtered.csv', 'bams/coverages/coverage_full.csv']

def join_nonempty(values):
    return ';'.join(values[values.notnull() & (values.str.len() > 0)])

def merge_folders(output_folder, folders, force, unique_sample_id_column, additional_coverage, enforce_integer_ids):
    data = {}

    #read in
    for folder in folders:
        log.info('Reading files from %s...', folder)

        for file in files_to_merge:     
            fn = os.path.join(folder, file)
            try:
                folder_data = pd.read_csv(fn, index_col = False, dtype = 'object')
                log.info('%s: shape = %s', fn, folder_data.shape)

                folder_data['Folder'] = folder
                if not 'Notes' in folder_data.columns:
                    folder_data['Notes'] = ''

                if file in data:
                    #reorder if we have to
                    if list(data[file].columns) != list(folder_data.columns):
                        log.warn('Inconsistent column names! Will proceed but column order may be changed')
                        log.warn('%s =\t%s', 'previous', ','.join(data[file].columns))
                        log.warn('%s =\t%s', 'this file', ','.join(folder_data.columns))
                    data[file] = pd.concat([data[file], folder_data], ignore_index = True)
                else:
                    data[file] = folder_data
            except pd.io.common.EmptyDataError:
                log.info('Skipping empty file: %s', fn)

    #add additional coverage data from another file -- to override coverage numbers with sanger sequencing
    if additional_coverage:
        fn = os.path.join(additional_coverage)
        try:
            folder_data = pd.read_csv(fn, index_col = False, dtype = 'object')
            log.info('%s: shape = %s', fn, folder_data.shape)

            folder_data['Folder'] = additional_coverage
            if not 'Notes' in folder_data.columns:
                folder_data['Notes'] = ''

            data['bams/coverages/coverage_full.csv'] = pd.concat([data['bams/coverages/coverage_full.csv'], folder_data], ignore_index = True)
        except pd.io.common.EmptyDataError:
            log.info('Empty additional coverage file: %s', additional_coverage)
            raise

    #handle duplicates
    if unique_sample_id_column:
        log.info('Combining duplicate values of %s into single row...', unique_sample_id_column)

        variants_null_sample_column = data['variants_raw/variants_summary_filtered.csv'][unique_sample_id_column].isnull()
        if enforce_integer_ids:
            #just remove the .0 from end of strings rather than trying to make them ints, which causes problems with NAs
            data['variants_raw/variants_summary_filtered.csv'][unique_sample_id_column] = data['variants_raw/variants_summary_filtered.csv'][unique_sample_id_column].replace(
                re.compile(r'\.0$'), '')
            log.info('Removing .0 from end of ID column!')
        if variants_null_sample_column.any():
            log.warn('Dropping %d variant rows without %s to make unique table', variants_null_sample_column.sum(), unique_sample_id_column)
            variants_with_sample = data['variants_raw/variants_summary_filtered.csv'][~variants_null_sample_column]
        else:
            variants_with_sample = data['variants_raw/variants_summary_filtered.csv']

        #merge variant table - don't need to fix dtype here since we just check for equality
        #in this case we can drop the previous index, since it's not meaningful (we may not even have to reset it at all, not sure...)
        data['variants_summary_filtered.unique.csv'] = variants_with_sample.drop_duplicates([unique_sample_id_column, 'Chr', 'Start', 'Ref', 'Alt']).reset_index(drop = True)
        log.info('Variants: Combined %d/%d rows into %d rows (first row of each sample kept)',
            len(variants_with_sample), len(data['variants_raw/variants_summary_filtered.csv']), len(data['variants_summary_filtered.unique.csv']))

        #merge coverage table
        coverage_group_cols = [unique_sample_id_column, 'Target']
        coverage_ignore_cols = ['basepairs', 'Sample']
        coverage_aggregation_ops = {'min_coverage': max, 'cov_per_bp': max, 'sum_coverage': max, 'fraction_zero_coverage': min,
            'Notes': join_nonempty, 'Folder': join_nonempty}
        #detect additional columns to join
        for colname in data['bams/coverages/coverage_full.csv'].columns:
            if not colname in coverage_group_cols + coverage_ignore_cols + list(coverage_aggregation_ops.keys()):
                coverage_aggregation_ops[colname] = join_nonempty
                log.info('Coverage: Detected additional column to join: %s', colname)

        #fix dtypes for numeric columns
        for colname in ['min_coverage', 'sum_coverage']:
            data['bams/coverages/coverage_full.csv'][colname] = data['bams/coverages/coverage_full.csv'][colname].astype(int)
        for colname in ['cov_per_bp', 'fraction_zero_coverage']:
            data['bams/coverages/coverage_full.csv'][colname] = data['bams/coverages/coverage_full.csv'][colname].astype(float)
        log.info('Coverage: Fixed datatypes of numeric columns for min/max')

        coverage_null_sample_column = data['bams/coverages/coverage_full.csv'][unique_sample_id_column].isnull()
        if enforce_integer_ids:
            #just remove the .0 from end of strings rather than trying to make them ints, which causes problems with NAs
            data['bams/coverages/coverage_full.csv'][unique_sample_id_column] = data['bams/coverages/coverage_full.csv'][unique_sample_id_column].replace(
                re.compile(r'\.0$'), '')
            log.info('Removing .0 from end of ID column!')
        if coverage_null_sample_column.any():
            log.warn('Dropping %d coverage rows without %s to make unique table', coverage_null_sample_column.sum(), unique_sample_id_column)
            coverage_with_sample = data['bams/coverages/coverage_full.csv'][~coverage_null_sample_column]
        else:
            coverage_with_sample = data['bams/coverages/coverage_full.csv']

        #actually do the aggregation
        #in this case we want to keep the index, since this is coverage_group_cols, which is dnaid + target (I think?)
        data['coverage_full.unique.csv'] = coverage_with_sample.groupby(coverage_group_cols).aggregate(coverage_aggregation_ops).reset_index(drop = False)
        #reorder columns
        coverage_unique_cols = [colname for colname in coverage_with_sample.columns if not colname in coverage_ignore_cols]
        assert pd.Series(data['coverage_full.unique.csv'].columns).isin(coverage_unique_cols).all()
        data['coverage_full.unique.csv'] = data['coverage_full.unique.csv'][coverage_unique_cols]
        #log
        log.info('Coverage: Combined %d/%d rows (mean min_coverage = %f) into %d rows (mean min_coverage = %f) taking best values per sample',
            len(coverage_with_sample), len(data['bams/coverages/coverage_full.csv']), coverage_with_sample['min_coverage'].mean(),
            len(data['coverage_full.unique.csv']), data['coverage_full.unique.csv']['min_coverage'].mean())
    else:
        log.info('Not combining rows since unique-sample-id-column not provided.')

    #write out
    log.info('Writing output files...')
    for file in data.keys():
        log.info('%s: final shape = %s', file, data[file].shape)
        bnfile = os.path.basename(file)

        fn = os.path.join(output_folder, bnfile)

        #handle existing files
        if os.path.exists(fn):
            if force:
                log.info('Overwriting file: %s', fn)
            else:
                raise Exception('Cannot write to file, exists: %s. Aborting, delete file or set --force to overwrite.' % fn)
        else:
            log.info('Creating new file: %s', fn)

        #add hash for first column name
        data[file].columns = ['#' + c if ix == 0 else c for ix, c in enumerate(data[file].columns)]
        #write to file
        data[file].to_csv(fn, index = False)

def main():
    log.info('Called with arguments: "%s"', '" "'.join(sys.argv))
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--force", help="overwrite existing files", action='store_true')
    parser.add_argument("-c", "--unique-sample-id-column", help="column which contains the unique sample id (to remove duplicates)")
    parser.add_argument("-a", "--additional-coverage", help="file with additional coverage data, which will be added to the merged coverages")
    parser.add_argument("--enforce-integer-ids", help="force all sample id columns to be in integer format", action='store_true')
    parser.add_argument("OUTPUT_FOLDER", help="output folder")
    parser.add_argument("INPUT_FOLDER", help="input folders to merge (analysis directories)", nargs="+")
    args = parser.parse_args()

    merge_folders(args.OUTPUT_FOLDER, args.INPUT_FOLDER, args.force, args.unique_sample_id_column, args.additional_coverage, args.enforce_integer_ids)

    return 0

if __name__ == '__main__':
    sys.exit(main())
