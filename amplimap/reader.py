# -*- coding: utf-8 -*-
"""
This module contains methods related to reading input files, such as targets.csv, probes.csv etc.
"""

import re
import sys
import os, hashlib
import pandas as pd
import Bio.Seq

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
log.addHandler(sh)

mipgen_columns = ['chr','ext_probe_sequence','lig_probe_sequence','mip_scan_start_position','mip_scan_stop_position','probe_strand','mip_name']
heatseq_columns = ['chromosome','ext_sequence','lig_sequence','target_start','target_stop','probe_strand','probe_id']
probe_columns = ['id', 'first_primer_5to3', 'second_primer_5to3', 'chr', 'target_start', 'target_end', 'strand']

class AmplimapReaderException(Exception):
    """
    Will be raised if reading one of the standard input files has failed, with the aim of providing a more useful error message in the user-visible output.
    """
    def __init__(self, e: Exception, filename: str, should_have_header: bool = None):
        if isinstance(e, ValueError):
            err = 'Invalid value. This usually indicates that there were non-numeric characters (eg. letters, spaces, ...) in a column that should only contain numbers.\n\nMessage: {}'.format(
                str(e)
            )
        elif isinstance(e, AssertionError):
            err = 'Assertion error. This usually indicates that one of the values is invalid or in an unexpected range.\n\nMessage: {}'.format(
                str(e)
            )
        else:
            err = str(e)

        header_note = ''
        if should_have_header is not None:
            if should_have_header:
                header_note = 'Please note that this file SHOULD start with a header line specifying the column names.'
            else:
                header_note = 'Please note that this file SHOULD NOT contain any header line.'

        #add better message to exception
        message = '\n\nError while reading file {}:\n{}\n\n{}\n\n'.format(
            filename,
            err,
            header_note
        )

        #init base exception
        super().__init__(message)

def get_code_versions(path: str = '.') -> dict:
    """
    Get the file modification times of common code files for versioning.
    """
    VERSIONS = {}
    VERSIONS['_amplimap'] = str(__version__)
    VERSIONS['Snakefile'] = str(os.path.getmtime(os.path.join(path, 'Snakefile')))
    VERSIONS['parse_reads'] = str(os.path.getmtime(os.path.join(path, 'parse_reads.py')))
    VERSIONS['pileup'] = str(os.path.getmtime(os.path.join(path, 'pileup.py')))
    VERSIONS['stats_alignment'] = str(os.path.getmtime(os.path.join(path, 'stats_alignment.py')))
    VERSIONS['coverage'] = str(os.path.getmtime(os.path.join(path, 'coverage.py')))
    VERSIONS['simulate'] = str(os.path.getmtime(os.path.join(path, 'simulate.py')))
    VERSIONS['variants'] = str(os.path.getmtime(os.path.join(path, 'variants.py')))
    return VERSIONS

def get_file_hashes(path: str = '.') -> dict:
    """
    Get SHA256 hashes for common input files, so we can make sure these didn't change between runs.
    """
    hashes = {}
    for fn in ['probes.csv', 'targets.bed', 'targets.csv', 'snps.txt']:
        fpath = os.path.join(path, fn)
        if os.path.isfile(fpath):
            file_hash = hashlib.sha256()
            with open(fpath, 'rb', buffering = 0) as file:
                for chunk in iter(lambda : file.read(128*1024), b''):
                    file_hash.update(chunk)
            hashes[fn] = file_hash.hexdigest()
        else:
            hashes[fn] = ''

    return hashes

def merge_probes_by_id(rows: pd.DataFrame) -> pd.Series:
    """
    Merge multiple probes.csv rows with the same probe ID together, to handle cases where MIPGEN generated multiple version because of a SNP.
    """
    #grab first row as template (use .loc instead of .iloc?)
    first_row = rows.iloc[0].copy()

    first_row['n_merged'] = len(rows)
    first_row['n_wildcards'] = None

    if len(rows) > 1:
        first_row['n_wildcards'] = 0

        #make sure these are all the same
        assert (rows['chr'] == first_row['chr']).all(), 'found multiple rows with same probe name but different chr: %s' % first_row['id']
        assert (rows['target_start'] == first_row['target_start']).all(), 'found multiple rows with same probe name but different target_start: %s' % first_row['id']
        assert (rows['target_end'] == first_row['target_end']).all(), 'found multiple rows with same probe name but different target_end: %s' % first_row['id']
        assert (rows['strand'] == first_row['strand']).all(), 'found multiple rows with same probe name but different strand: %s' % first_row['id']

        #merge primers together
        for primer_col in ['first_primer_5to3', 'second_primer_5to3']:
            assert (rows[primer_col].str.len() == len(first_row[primer_col])).all(), 'found multiple rows with same probe id but primer sequences were different lengths'
            for ixrow in range(1, len(rows)):
                for ixchar in range(len(first_row[primer_col])):
                    #if this row doesn't match the first row's sequence
                    if first_row[primer_col][ixchar] != rows.iloc[ixrow][primer_col][ixchar]:
                        #replace character by dot, but only if it hasn't been replaced yet
                        if first_row[primer_col][ixchar] != '.':
                            first_row['n_wildcards'] += 1
                            first_row[primer_col] = first_row[primer_col][:ixchar] + '.' + first_row[primer_col][(ixchar+1):]

        #make sure we actually replaced something with a wildcard
        assert first_row['n_wildcards'] > 0, 'found multiple rows with same probe name but primer sequences were identical: %s' % first_row['id']
        assert first_row['n_wildcards'] < 10, 'found too many differences between primer sequences: %s' % first_row['id']

    del first_row['id'] #remove this, will be attached again later
    return first_row

def read_and_convert_mipgen_probes(path: str) -> pd.DataFrame:
    """
    Read probes file from MIPGEN in CSV format and generate an amplimap probes.csv from it.
    """
    try:
        design = pd.read_csv(path)

        log.info('Read mipgen table from %s -- found %d probes', path, len(design))
        original_len = len(design)
        design.dropna(axis=0, how='all', inplace=True)
        if original_len != len(design):
            log.info('%d probes after dropping empty rows', len(design))
        assert len(design) > 0, 'mipgen table is empty!'

        #make sure we have all the columns we expect
        for column in mipgen_columns:
            if not column in design.columns:
                raise Exception('invalid mipgen table: column %s is missing!' % column)

        #rename columns
        design.rename(columns = {'mip_scan_start_position': 'target_start', 'mip_scan_stop_position': 'target_end', 'probe_strand': 'strand', 'mip_name': 'id'}, inplace=True)

        #in smMIPs, the first arm (read1) seems to be the reverse primer, which would normally be the second arm
        #its sequence in the table seems to be on the forward strand, which needs to be RC'd to get the actual read sequence (5p->3p)
        design['first_primer_5to3'] = [str(Bio.Seq.Seq(s).reverse_complement()) for s in design['lig_probe_sequence'].str.upper()]
        design['second_primer_5to3'] = design['ext_probe_sequence'].str.upper()

        #detect SNPs and replace them by dot
        original_len = len(design)
        design = design.groupby('id').apply(merge_probes_by_id)
        if original_len != len(design):
            log.info('%d/%d probes left after merging probes by ID', len(design), original_len)
        #reset index to turn id into column again
        design.reset_index(inplace=True)  

        #reorder and extract columns
        design = design[probe_columns]
    except Exception as e:
        raise AmplimapReaderException(e, filename = path,  should_have_header = True).with_traceback(sys.exc_info()[2])

    return design

def read_and_convert_heatseq_probes(path: str) -> pd.DataFrame:
    """
    UNTESTED: Read probes file from Roche heatseq in CSV format and generate an amplimap probes.csv from it.
    """
    try:
        design = pd.read_csv(path, sep='\t', comment='#')

        log.info('Read heatseq table from %s -- found %d probes', path, len(design))
        original_len = len(design)
        design.dropna(axis=0, how='all', inplace=True)
        if original_len != len(design):
            log.info('%d probes after dropping empty rows', len(design))
        assert len(design) > 0, 'heatseq table is empty!'

        #make sure we have all the columns we expect
        for column in heatseq_columns:
            if not column in design.columns:
                raise Exception('invalid heatseq table: column %s is missing!' % column)

        #rename columns
        design.rename(columns = {'chromosome': 'chr', 'target_stop': 'target_end', 'probe_strand': 'strand', 'probe_id': 'id'}, inplace=True)

        #In this protocol, read1 should contain extension primer and read2 the ligation primer
        #Note that this is the opposite of smMIPs!
        design['second_primer_5to3'] = [str(Bio.Seq.Seq(s).reverse_complement()) for s in design['lig_sequence'].str.upper()]
        design['first_primer_5to3'] = design['ext_sequence'].str.upper()

        #detect SNPs and replace them by dot
        original_len = len(design)
        design = design.groupby('id').apply(merge_probes_by_id)
        if original_len != len(design):
            log.info('%d/%d probes left after merging probes by ID', len(design), original_len)
        #reset index to turn id into column again
        design.reset_index(inplace=True)  

        #reorder and extract columns
        design = design[probe_columns]
    except Exception as e:
        raise AmplimapReaderException(e, filename = path,  should_have_header = True).with_traceback(sys.exc_info()[2])

    return design

def read_new_probe_design(path: str, reference_type: str = 'genome') -> pd.DataFrame:
    """
    Read amplimap probes.csv file and return pandas dataframe.
    """
    try:
        design = pd.read_csv(path)

        log.info('Read probe design table from %s -- found %d probes', path, len(design))
        if list(design.columns) == mipgen_columns:
            #NB: smmip data seems to be in F2R1 orientation (second read = fwd in genomic coordinates) for fwd probes
            #but F1R2 orientation (second read = rev) for rev probes.
            #cs-tag data seems to be in F1R2 orientation for fwd targets. unclear for rev targets, but presumably F2R1?
            #in other words, CS-tag is in gene orientation, while smMIP is in opposite
            #so both are swapped for MIPs.
            #is this why sequences in probes.csv are currently so confusing?
            log.info('Detected old MIPGEN format, converting...')

            #read the probes file again in old mipgen format and convert
            design = read_and_convert_mipgen_probes(path)

        design = process_probe_design(design, reference_type)
    except Exception as e:
        raise AmplimapReaderException(e, filename = path,  should_have_header = True).with_traceback(sys.exc_info()[2])

    return design

def process_probe_design(design: pd.DataFrame, reference_type: str = 'genome') -> pd.DataFrame:
    """
    Read amplimap probes.csv file and return pandas dataframe.
    """
    #now drop empty rows
    original_len = len(design)
    design.dropna(axis=0, how='all', inplace=True)
    if original_len != len(design):
        log.info('%d probes after dropping empty rows', len(design))
    assert len(design) > 0, 'Probe design table is empty!'

    #make sure we have all the columns we expect
    for column in probe_columns:
        if not column in design.columns:
            raise Exception('invalid probe design table: column %s is missing!' % column)

    #fix data types
    design = design.astype(
        {
            'id': str,
            'chr': str,
            'strand': str,
            'first_primer_5to3': str,
            'second_primer_5to3': str,
            #these need to be int instead of uint, otherwise they'll get turned into floats when subtracting an int (as we do below)
            'target_start': int,
            'target_end': int,
        },
        errors = 'raise'
    )

    #calculate some additional columns
    design.loc[design['strand'] == '+', 'probe_start'] = design.loc[design['strand'] == '+', 'target_start'] - design.loc[design['strand'] == '+', 'first_primer_5to3'].str.len()
    design.loc[design['strand'] == '-', 'probe_start'] = design.loc[design['strand'] == '-', 'target_start'] - design.loc[design['strand'] == '-', 'second_primer_5to3'].str.len()
    design.loc[design['strand'] == '+', 'probe_end'] = design.loc[design['strand'] == '+', 'target_end'] + design.loc[design['strand'] == '+', 'second_primer_5to3'].str.len()
    design.loc[design['strand'] == '-', 'probe_end'] = design.loc[design['strand'] == '-', 'target_end'] + design.loc[design['strand'] == '-', 'first_primer_5to3'].str.len()
    design['target_start_0'] = design['target_start'] - 1 #we have 1-based coords now
    design['probe_start_0'] = design['probe_start'] - 1 #we have 1-based coords now
    design['capture_size'] = design['probe_end'] - design['probe_start_0']
    design['target_length'] = design['target_end'] - design['target_start_0']

    invalid_coords = ~(design['target_length'] > 0)
    if invalid_coords.any():
        print(design.loc[invalid_coords, ['id', 'chr', 'probe_start', 'probe_end', 'strand']])
        raise Exception('Encountered probe regions with invalid target coordinates. Please check the input file and the output above. Note that the probe_end position must always be greater than the probe_start.')

    # #chr should start with chr
    # if reference_type == 'genome':
    #     design.loc[~design['chr'].str.startswith('chr'), 'chr'] = ['chr' + c for c in design.loc[~design['chr'].str.startswith('chr'), 'chr']]

    #check that probe ids don't contain unexpected characters (space breaks bam files for example, since umi will be added as tag)
    assert design['id'].str.contains(r'^[a-zA-Z0-9._:+-]*$').all(), 'probe names may only contain alphanumeric characters and ._:+-'

    #check for dupes
    any_dup = design.duplicated(subset = ['chr', 'probe_start', 'probe_end', 'strand'], keep = False)
    if any_dup.any():
        print('These lines in the design are duplicated:')
        print(design.loc[any_dup, ['chr', 'probe_start', 'probe_end', 'id', 'strand']])

        dup = design.duplicated(subset = ['chr', 'probe_start', 'probe_end', 'strand'], keep = 'first')
        log.warn('Dropping %d/%d duplicated rows', dup.sum(), len(design))
        design = design[~dup]
        log.warn('Now have %d rows after dropping duplicates', len(design))

        raise Exception('Stopping because of duplicate rows in probes.csv!')

    #make sure probes are uppercase
    design['first_primer_5to3'] = design['first_primer_5to3'].str.upper()
    design['second_primer_5to3'] = design['second_primer_5to3'].str.upper()

    #check strand
    if 'strand' in design.columns:
        assert (design['strand'].isin(['+', '-'])).all(), 'strand must be + or -'

    #get reverse complement of sequence for smart trimming and matching to arms reverse read
    for c in ['first_primer_5to3', 'second_primer_5to3']:
        design['%s__rc' % c] = [str(Bio.Seq.Seq(s).reverse_complement()) for s in design[c]]

    #set the index to the probe id (useful for probe stats)
    design.set_index('id', drop = False, inplace = True, verify_integrity = True)

    return design

def read_sample_info(path: str) -> pd.DataFrame:
    """
    Read amplimap sample_info.csv file and return pandas dataframe.
    """
    try:
        sample_infos = pd.read_csv(path, dtype = 'object')

        assert not sample_infos.columns.str.contains('^Unnamed').any(), \
            'Missing column names in sample info table! Please provide names for all columns in the first row.'

        log.info('Read sample info table from %s -- found %d rows', path, len(sample_infos))
        sample_infos.dropna(axis=0, how='all', inplace=True)
        log.info('%d after dropping empty rows', len(sample_infos))

        assert sample_infos.columns[0] == 'Sample' and sample_infos.columns[1] == 'Targets', \
            'Missing Sample/Targets columns! Got %d columns: %s' % (len(sample_infos.columns), '/'.join(sample_infos.columns))

        #split up Targets column into list
        sample_infos['Target'] = sample_infos['Targets'].str.split(';')

        #make new data frame by stacking the Target lists
        #inspired by: https://stackoverflow.com/questions/38372016/split-nested-array-values-from-pandas-dataframe-cell-over-multiple-rows
        non_target_columns = [c for c in sample_infos.columns if c != 'Target']
        sample_infos.set_index(non_target_columns, inplace = True) #this needs need to be all unstacked columns
        #stack each column (only Target actually)
        sample_infos = sample_infos.apply(lambda col: col.apply(pd.Series).stack(), axis=0) 
        #turn multilevel index back into columns
        sample_infos.reset_index(inplace = True)
        #drop Targets column and target index column generated by reset_index
        sample_infos.drop(['Targets', 'level_%d' % len(non_target_columns)], axis=1, inplace = True)

        #check for dupes
        if sample_infos.duplicated(['Sample', 'Target']).any():
            print(targets[sample_infos.duplicated(['Sample', 'Target'], keep=False)])
            raise Exception('Invalid sample_infos table -- found duplicated Sample/Target IDs!')

        #set proper index
        sample_infos.set_index(['Sample', 'Target'], drop = True, inplace = True, verify_integrity = True)

        print(sample_infos.head(3))
    except Exception as e:
        raise AmplimapReaderException(e, filename = path,  should_have_header = True).with_traceback(sys.exc_info()[2])

    return sample_infos

def read_targets(path: str, check_overlaps: bool = False, reference_type: str = 'genome', file_type: str = 'bed') -> pd.DataFrame:
    """
    Read amplimap targets.csv or targets.bed file and return pandas dataframe.
    """
    try:
        if file_type == 'bed':
            skiprows = 0
            with open(path, 'r') as f:
                for line in f:
                    if line.startswith('track') or line.startswith('browser'):
                        skiprows += 1
                    else:
                        break

            if skiprows > 0:
                log.info('Skipping %d browser/track rows from target BED files', skiprows)

            targets = pd.read_table(path, header = None, skiprows = skiprows)
            log.info('Read targets from %s: %d rows', path, len(targets))

            if targets.shape[1] == 6:
                targets.columns = ['chr', 'start_0', 'end', 'id', 'score', 'strand']
                #targets = targets[['chr', 'start_0', 'end', 'id']] #remove score/strand
            elif targets.shape[1] == 4:
                targets.columns = ['chr', 'start_0', 'end', 'id']
            elif targets.shape[1] == 3:
                targets.columns = ['chr', 'start_0', 'end']
                targets['id'] = ['target_%d' % i for i in range(len(targets))]
            else:
                print(targets)
                raise Exception('Invalid target BED file: expected to find either 6, 4 or 3 columns separated by tabs.')
        elif file_type == 'csv':
            targets = pd.read_csv(path)

            #check columns
            if not ('chr' in targets.columns and 'start' in targets.columns and 'end' in targets.columns):
                print(targets)
                raise Exception('Invalid target CSV: first line should contain column names, which must include "chr", "start" and "end".')

            #add id if need be
            if not 'id' in targets.columns:
                targets['id'] = ['target_%d' % i for i in range(len(targets))]

            #fix start and calculate start_0
            targets['start'] = targets['start'].astype(int)
            targets['start_0'] = targets['start'] - 1
        else:
            raise Exception('invalid target file type')

        #fix data types
        targets = targets.astype(
            {
                'id': str,
                'chr': str,
                'start_0': int,
                'end': int,
            },
            errors = 'raise'
        )

        #make sure we have no duplicated ids
        if targets.duplicated('id').any():
            print(targets[targets.duplicated('id', keep=False)])
            raise Exception('Invalid target table: found duplicated target ID(s)!')

        #check for overlaps
        if check_overlaps:
            targets_dict = targets.to_dict()
            for ixrow1 in range(len(targets_dict['chr'])-1):
                for ixrow2 in range(ixrow1+1, len(targets_dict['chr'])):
                    if targets_dict['chr'][ixrow1] == targets_dict['chr'][ixrow2] \
                    and targets_dict['start_0'][ixrow1] <= targets_dict['end'][ixrow2] \
                    and targets_dict['end'][ixrow1] >= targets_dict['start_0'][ixrow2]:
                        print(targets.iloc[ixrow1])
                        print(targets.iloc[ixrow2])
                        raise Exception(
                            'Invalid target table: target regions %s and %s overlap! Please merge them or shorten one of the regions' % (
                                targets_dict['id'][ixrow1],
                                targets_dict['id'][ixrow2]
                            )
                        )

        #chr should always start with chr
        targets['chr_original'] = targets['chr']
        # if reference_type == 'genome':
        #     targets.loc[~targets['chr'].str.startswith('chr'), 'chr'] = ['chr' + c for c in targets.loc[~targets['chr'].str.startswith('chr'), 'chr']]

        targets['length'] = targets['end'] - targets['start_0']

        targets_wrong_length = targets['length'] <= 0
        if targets_wrong_length.any():
            print(targets[targets_wrong_length])
            raise Exception('Invalid target table: found targets with length <= 0! Please check start and end coordinates.')

        targets['type'] = 'target'
    except Exception as e:
        raise AmplimapReaderException(e, filename = path,  should_have_header = file_type == 'csv').with_traceback(sys.exc_info()[2])

    return targets

def write_targets_bed(path: str, targets: pd.DataFrame):
    """
    Write targets dataframe to bed file.
    """
    if not 'score' in targets.columns:
        targets['score'] = '.'
    if not 'strand' in targets.columns:
        targets['strand'] = '.'

    targets[['chr_original', 'start_0', 'end', 'id', 'score', 'strand']].to_csv(path, header = False, index = False, sep = '\t')

def read_snps_txt(path: str, reference_type: str = 'genome') -> pd.DataFrame:
    """
    Read amplimap snps.txt file and return pandas dataframe.
    """
    try:
        log.info('Loading SNPs from %s', path)
        snps = pd.read_table(path, header = None) #pos is 1-based!

        #detect potential issues with spaces instead of tabs
        if snps.shape[1] == 1:
            if ' ' in snps.iloc[0, 0]:
                log.info('SNP table does not seem to be tab-separated, trying again with spaces')
                snps = pd.read_table(path, header = None, sep='\s+')

        if snps.shape[1] == 5:
            snps.columns = ['chr', 'pos', 'id', 'snp_ref', 'snp_alt']
        elif snps.shape[1] == 4:
            snps.columns = ['chr', 'pos', 'id', 'snp_ref']
        elif snps.shape[1] == 3:
            snps.columns = ['chr', 'pos', 'id']
        elif snps.shape[1] == 2:
            snps.columns = ['chr', 'pos']
            snps['id'] = ['snp_%d' % i for i in range(len(snps))]
        else:
            print(snps)
            raise Exception('Invalid SNP table -- expected to find between 2 and 5 columns.')

        #fix data types
        snps = snps.astype(
            {
                'id': str,
                'chr': str,
                'pos': int,
            },
            errors = 'raise'
        )

        # #chr should always start with chr    
        # if reference_type == 'genome':
        #     snps.loc[~snps['chr'].str.startswith('chr'), 'chr'] = ['chr' + c for c in snps.loc[~snps['chr'].str.startswith('chr'), 'chr']]

        #do we have genotypes for a snp?
        snps['snp_has_genotypes'] = False
        if 'snp_ref' in snps and 'snp_alt' in snps:
            snps.loc[snps['snp_ref'].notnull() & snps['snp_alt'].notnull(), 'snp_has_genotypes'] = True

        #stop if we have duplicate chr/pos
        if snps.duplicated(['chr', 'pos']).any():
            print(snps[snps.duplicated(['chr', 'pos'], keep=False)])
            raise Exception('Invalid SNPs table -- found duplicated chromosome/position!')

        #set the index to the snp id (stops if duplicate IDs)
        snps.set_index('id', drop = False, inplace = True, verify_integrity = True)

        log.info('Found SNP info for %d SNPs -- shape = %s', len(snps), str(snps.shape))
    except Exception as e:
        raise AmplimapReaderException(e, filename = path,  should_have_header = True).with_traceback(sys.exc_info()[2])

    return snps
