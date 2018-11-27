# -*- coding: utf-8 -*-
"""
This module contains methods to parse, merge, annotate and filter variant calls.
"""

import pandas as pd
import numpy as np
import re
import os

#for overlap with targest bed
import collections
import interlap

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
log.addHandler(sh)

#functions for reading input files
from .reader import read_targets, read_sample_info

#datatypes for variant summary tables
variants_dtypes = {
    'Ref': str,
    'Alt': str,
    'Start': int,
    'End': int,
    'GeneDetail.refGene': str,
    'Func.refGene': str,
    'ExonicFunc.refGene': str,
    'genomicSuperDups': str,
    'SIFT_pred': str,
    'Polyphen2_HDIV_pred': str,
    'LRT_pred': str,
    'MutationTaster_pred': str,
    'MutationAssessor_pred': str,
    'FATHMM_pred': str,
    'GERP++_RS': float,
    'phyloP100way_vertebrate': float,
}
variants_na_values = ['', 'NA', '.']

def load_gene_exons(file: str, genes: list) -> dict:
    """Load list of exon chr, strand, start and end locations for each gene in genes."""
    genes = genes.unique()

    #load gene exon information from flat table
    print('Loading exon table from ', file)
    gene_exon_df = pd.read_table(file, header = None)
    gene_exon_df.rename(columns={2: 'chr', 3: 'strand', 8: 'n_exons', 9: 'exon_starts', 10: 'exon_ends', 12: 'gene'}, inplace=True)
    print('Loaded ', len(gene_exon_df), ' rows.')

    #filter down to only the genes we actually need
    gene_exon_df = gene_exon_df[ gene_exon_df['gene'].isin(genes) ]
    print('Filtered to ', len(gene_exon_df), ' rows for ', len(genes), ' genes.')

    #identify standard chrs
    regex_standard_chr = re.compile(r'^chr([0-9MXY]+)$')

    #make dict of genes to exon locations
    gene_exons = {}
    for gex in gene_exon_df.itertuples():
        try:
            #skip chrY FIXME: do we want this
            if gex.chr == 'chrY':
                continue

            exon_data = {
                'chr': gex.chr,
                'strand': gex.strand,
                'starts': [int(x) for x in gex.exon_starts.split(',') if len(x) > 0],
                'ends': [int(x) for x in gex.exon_ends.split(',') if len(x) > 0],
            }

            #do we already have exons for this gene?
            if gex.gene in gene_exons:
                #is this on the same chr (and not an Alt chr)?
                if gene_exons[gex.gene]['chr'] == exon_data['chr']:
                    #if we are on the same strand let's merge them
                    if gene_exons[gex.gene]['strand'] == exon_data['strand']:
                        n_total = 0
                        n_added = 0
                        for start, end in zip(exon_data['starts'], exon_data['ends']):
                            n_total += 1

                            #try to find existing pair that matches                        
                            for existing_start, existing_end in zip(gene_exons[gex.gene]['starts'], gene_exons[gex.gene]['ends']):
                                if existing_start == start and existing_end == end:
                                    break
                            #otherwise add it
                            else:
                                n_added += 1
                                gene_exons[gex.gene]['starts'].append(start)
                                gene_exons[gex.gene]['ends'].append(end)

                        #log.debug('Merging exon data from duplicate rows for gene %s - added %d/%d', gex.gene, n_added, n_total)
                    else:
                        # print(gex.gene)
                        # print(gene_exons[gex.gene])
                        # print(exon_data)
                        raise Exception('Exon strand mismatch for {}'.format(gex.gene))
                #if the existing copy is on a standard chr we want to keep it
                elif regex_standard_chr.match(gene_exons[gex.gene]['chr']):
                    #normally this should only happen if we are on an alt chr now
                    if not regex_standard_chr.match(gex.chr):
                        # log.debug('Skipping copy of gene on nonstandard chr for gene exon table: {} @ {}/{}'.format(
                        #     gex.gene, gex.chr, gene_exons[gex.gene]['chr']
                        # ))
                        continue
                    else:
                        raise Exception('Two copies of gene on different standard chrs encounted for gene exon table: {} @ {}/{}'.format(
                            gex.gene, gex.chr, gene_exons[gex.gene]['chr']
                        ))
                #if the existing copy was not on a standard chr we will just overwrite it with this copy
                else:
                    pass
                    # log.debug('Overwriting copy of gene on nonstandard chr for gene exon table: {} @ {}/{}'.format(
                    #     gex.gene, gex.chr, gene_exons[gex.gene]['chr']
                    # ))

            #error checking
            assert len(exon_data['starts']) == len(exon_data['ends'])
            assert len(exon_data['starts']) == gex.n_exons
            for start, end in zip(exon_data['starts'], exon_data['ends']):
                assert end > start, 'Exon end is not after start for {}'.format(gex.gene)

            #finally save the exon data
            gene_exons[gex.gene] = exon_data
        except Exception as e:
            print('Igoring error while processing exon data: {}'.format(
                str(e)))

    #make sure we found exons for every gene
    missing_genes = set(genes) - set(gene_exons.keys())
    if len(missing_genes) > 0:
        print('Failed to find exon location information for: %s' % (','.join(missing_genes)))

    return gene_exons

def find_closest_exon(row: pd.Series, gexs: dict) -> int:
    """Get distance of variant to the closest exon of the gene it has been annotated with."""

    if row['Gene.refGene'] is None:
        return None
    if not row['Func.refGene'] in ['intronic', 'upstream', 'downstream']:
        return None
    if not row['Gene.refGene'] in gexs:
        return None

    gex = gexs[row['Gene.refGene']]
    assert row['Chr'] == gex['chr']
    
    #upstream of start = negative, downstream of end = positive
    start_dist = [row['End'] - boundary_pos for boundary_pos in gex['starts'] if boundary_pos >= row['End']]
    end_dist = [row['Start'] - boundary_pos for boundary_pos in gex['ends'] if boundary_pos <= row['Start']]
    
    #find smallest abs(value) for each
    closest_start = max(start_dist) if len(start_dist) > 0 else None
    closest_end = min(end_dist) if len(end_dist) > 0 else None
    
    #figure out whether we are closer to start or end
    if closest_start is None and closest_end is None:
        closest_distance = None
    elif closest_start is None:
        closest_distance = closest_end
    elif closest_end is None:
        closest_distance = closest_start
    else:
        closest_distance = closest_start if abs(closest_start) < closest_end else closest_end
    
    #if gene strand is negative it's actually the opposite
    if closest_distance is None:
        return None
    elif gex['strand'] == '+':
        return (closest_distance)
    elif gex['strand'] == '-':
        return (closest_distance * -1)
    else:
        raise Exception('Invalid strand')

#see: /hts/data6/smcgowan/hts_scripts/TableAnnovar_to_spreadsheet.pl
def calculate_del_score(merged: pd.DataFrame):
    """
    Add a column DeleteriousScore to dataframe which contains a count of how many tools have assigned this variant a deletious scores.

    Score ranges from 0-6, corresponding to the tools SIFT, Polyphen2, LRT, MutationTaster, GERP++ and phyloP100way_vertebrate.

    Additionally, any stopgain, frameshift or splicing variants are always set to 6.
    """
    assert merged['GERP++_RS'].dtype == float
    assert merged['phyloP100way_vertebrate'].dtype == float
        
    merged['DeleteriousScore'] = 0
    merged.loc[merged['SIFT_pred'] == 'D', 'DeleteriousScore'] += 1
    merged.loc[merged['Polyphen2_HDIV_pred'] == 'D', 'DeleteriousScore'] += 1
    merged.loc[merged['LRT_pred'] == 'D', 'DeleteriousScore'] += 1
    merged.loc[merged['MutationTaster_pred'].isin(['A', 'D']), 'DeleteriousScore'] += 1
    merged.loc[merged['GERP++_RS'] >= 5, 'DeleteriousScore'] += 1
    merged.loc[merged['phyloP100way_vertebrate'] > 1, 'DeleteriousScore'] += 1
    #set score = 6 for some special cases
    merged.loc[merged['ExonicFunc.refGene'].notnull() & merged['ExonicFunc.refGene'].str.contains('^stopgain|^frameshift'), 'DeleteriousScore'] = 6
    merged.loc[merged['Func.refGene'].notnull() & merged['Func.refGene'].str.contains('splic'), 'DeleteriousScore'] = 6
    merged.loc[merged['ExonicFunc.refGene'].notnull() & merged['ExonicFunc.refGene'].str.contains('splic'), 'DeleteriousScore'] = 6

def merge_variants(input, output):
    """Merge individual Annovar CSV files together."""

    merged = None
    for file in input:
        sname = os.path.basename(file).split('.')[0]
        print('Reading', file, ' ( Sample =', sname, ')...')
        try:
            df = pd.read_csv(file, index_col = False, dtype = variants_dtypes, na_values = variants_na_values)

            print('Data shape:', str(df.shape))
            if len(df) > 0:
                df['Sample'] = sname

                if merged is None:
                    merged = df
                else:
                    merged = merged.append(df, ignore_index = True)
            else:
                print ('Ignoring - empty!')
        except pd.io.common.EmptyDataError:
            print('No data for', file, ', skipping.')

    if merged is None:
        print('No variants found in any of the samples, creating empty file!')
        open(output[0], 'a').close()
    else:
        print('Merged data shape:', str(merged.shape))
        merged.to_csv(output[0], index = False)

def make_summary(input: list, output: list, config: dict, exon_table_path: str = None):
    """Load merged Annovar CSV file (plus targets and sample info), process them and output a new CSV file."""

    #load targets (to add a targets column)
    targets = None
    if len(input['targets']) > 0:
        targets = read_targets(input['targets'], file_type = 'bed', reference_type = 'genome')

    #load sample information table
    sample_info = None
    if len(input['sample_info']) > 0:
        sample_info = read_sample_info(input['sample_info'][0])

    #load merged variant table, process it and write to CSV
    try:
        merged = pd.read_csv(input['merged'][0], index_col = False, dtype = variants_dtypes, na_values = variants_na_values)
        merged = make_summary_dataframe(
            merged,
            targets,
            sample_info,
            genome_name = config['general']['genome_name'],
            include_gbrowse_links = 'include_gbrowse_links' in config['annotate'] and config['annotate']['include_gbrowse_links'],
            include_exon_distance = 'include_exon_distance' in config['annotate'] and config['annotate']['include_exon_distance'],
            include_score = 'include_score' in config['annotate'] and config['annotate']['include_score'],
            exon_table_path = exon_table_path,
        )
        merged.to_csv(output[0], index = False)
    except pd.io.common.EmptyDataError:
        print('No variants found, creating empty file!')
        open(output[0], 'a').close()

def make_summary_dataframe(
        merged: pd.DataFrame,
        targets: pd.DataFrame = None,
        sample_info: pd.DataFrame = None,
        genome_name: str = None,
        include_gbrowse_links: bool = False,
        include_exon_distance: bool = False,
        include_score: bool = False,
        exon_table_path: str = None,
    ) -> pd.DataFrame:
    """Process merged Annovar dataframe (with optional targets and sample_info data frames) into a large summary table."""

    #process targets (to add a targets column)
    target_intervals = None
    if targets is not None:
        target_intervals = collections.defaultdict(interlap.InterLap)
        for target in targets.itertuples():
            target_intervals[target.chr].add( (int(target.start_0), int(target.end), target) ) #note the double parentheses!

    #process sample information table
    sample_info_columns = []
    if sample_info is not None:
        sample_info_columns = list(sample_info.columns)

    vcf = merged['Otherinfo'].apply(lambda x: pd.Series(x.split('\t')))
    vcf.replace('.', np.nan, inplace=True) #replace dots with NAs again
    vcf.columns = ['Chr', 'Pos', 'ID', 'Ref', 'Alt', 'Qual', 'Filter', 'Info', 'Fields', 'SampleData']
    assert (merged['Chr'] == vcf['Chr']).all(), 'chr mismatch'

    #We are using the approach from /hts/data6/smcgowan/hts_scripts/TableAnnovar_to_spreadsheet.pl, which only looks for exons
    #for the given gene (as annotated by Annovar). This should work fine, but is a bit inefficient.
    #See here for better, more general solutions: https://www.biostars.org/p/53561/
    #or actually just use the same approach as for targest (interlap.InterLap)
    if include_exon_distance:
        assert 'Gene.refGene' in merged.columns, 'The Annovar refGene database is required to calculate the distance to the nearest exon.'

        merged['DistanceNearestExon'] = merged.apply(find_closest_exon, axis = 1, gexs = load_gene_exons(
            exon_table_path,
            merged.loc[merged['Func.refGene'].isin(['intronic', 'upstream', 'downstream']), 'Gene.refGene']))

    #get gt info
    vcf_sample = pd.DataFrame([ dict( zip( row.Fields.split(':'), row.SampleData.split(':') ) ) for row in vcf.itertuples() ])
    vcf_sample.replace('.', np.nan, inplace=True) #replace dots with NAs again

    #parse out some fields
    vcf_infos = vcf['Info'].apply(lambda x: pd.Series(dict( [ tuple(pair.split('=')) if '=' in pair else (pair, True) for pair in x.split(';') ] )))
    #handle multiple alt alleles (not usually expected)
    if 'TR' in vcf_infos:
        vcf_infos.loc[vcf_infos['TR'].str.contains(','), 'TR'] = '-1' 
        vcf_infos.loc[vcf_infos['NF'].str.contains(','), 'NF'] = '-1' 
        vcf_infos.loc[vcf_infos['NR'].str.contains(','), 'NR'] = '-1' 
    if 'NR' in vcf_sample:
        vcf_sample.loc[vcf_sample['NR'].str.contains(','), 'NR'] = '-1' #also need to fix this, for assertion below
        vcf_sample.loc[vcf_sample['NV'].str.contains(','), 'NV'] = '-1' #also need to fix this, for assertion below
    #now rename coverage fields and convert to ints
    have_fwd_ref = False
    vcf_infos['Total_Coverage'] = -1
    vcf_infos['Total_Coverage_fwd'] = -1
    vcf_infos['Total_Coverage_rev'] = -1
    vcf_infos['Total_Reads_Alt'] = -1
    vcf_infos['Total_Reads_Alt_fwd'] = -1
    vcf_infos['Total_Reads_Alt_rev'] = -1

    #PLATYPUS will give us TC/TR
    if 'TC' in vcf_infos:
        have_fwd_ref = True
        vcf_infos['Total_Coverage'] = vcf_infos['TC'].astype(int)
        vcf_infos['Total_Coverage_fwd'] = vcf_infos['TCF'].astype(int)
        vcf_infos['Total_Coverage_rev'] = vcf_infos['TCR'].astype(int)
    if 'TR' in vcf_infos:
        vcf_infos['Total_Reads_Alt'] = vcf_infos['TR'].astype(int)
        vcf_infos['Total_Reads_Alt_fwd'] = vcf_infos['NF'].astype(int)
        vcf_infos['Total_Reads_Alt_rev'] = vcf_infos['NR'].astype(int)

    #GATK will give us an AD value
    if 'AD' in vcf_sample:
        have_fwd_ref = False
        vcf_infos['Total_Reads_Ref'] = [int(x.split(',')[0]) for x in vcf_sample['AD']]
        vcf_infos['Total_Reads_Alt'] = [int(x.split(',')[1]) for x in vcf_sample['AD']]
        if 'DP' in vcf_sample:
            #gatk will give us DP
            vcf_infos['Total_Coverage'] = vcf_sample['DP'].astype(int)
        else:
            #mutect2 won't
            vcf_infos['Total_Coverage'] = [int(x.split(',')[0]) + int(x.split(',')[1]) for x in vcf_sample['AD']]

    #some error checking
    if ((vcf_infos['Total_Coverage_fwd'] != -1) & (vcf_infos['Total_Coverage_rev'] != -1)).any():
        assert (vcf_infos['Total_Coverage'] == vcf_infos['Total_Coverage_fwd'] + vcf_infos['Total_Coverage_rev']).all(), 'total coverage mismatch'
    if ((vcf_infos['Total_Reads_Alt_fwd'] != -1) & (vcf_infos['Total_Reads_Alt_rev'] != -1)).any():
        assert ((vcf_infos['Total_Reads_Alt'] == -1) | (vcf_infos['Total_Reads_Alt'] == vcf_infos['Total_Reads_Alt_fwd'] + vcf_infos['Total_Reads_Alt_rev'])).all(), 'variant coverage mismatch'
    if 'NR' in vcf_sample:
        assert ((vcf_infos['Total_Reads_Alt'] == -1) | (vcf_infos['Total_Coverage'] == vcf_sample['NR'].astype(int))).all(), 'coverage match'
    if 'NV' in vcf_sample:
        assert ((vcf_infos['Total_Reads_Alt'] == -1) | (vcf_infos['Total_Reads_Alt'] == vcf_sample['NV'].astype(int))).all(), 'variant coverage match'

    #add info from VCF to table
    merged['Var_Zygosity'] = ['HOM' if gt == '1/1' else 'Het' if gt in ['0/1', '1/0'] else '???' for gt in vcf_sample['GT']]
    merged['Var_FailedFilters'] = ['' if (filt == 'PASS' or filt == '.') else filt for filt in vcf['Filter']]
    #mutect2 gives us AF
    if 'AF' in vcf_sample:
        merged['Var_AltFraction'] = [float(x.split(',')[0]) for x in vcf_sample['AF']]
    else:
        merged['Var_AltFraction'] = 1.0 * vcf_infos['Total_Reads_Alt'] / vcf_infos['Total_Coverage']
        merged.loc[vcf_infos['Total_Reads_Alt'] == -1, 'Var_AltFraction'] = None
    merged['Var_TotalCoverage'] = vcf_infos['Total_Coverage']
    if have_fwd_ref:
        merged['Var_Ref_fwd'] = (vcf_infos['Total_Coverage_fwd'] - vcf_infos['Total_Reads_Alt_fwd'])
        merged['Var_Ref_rev'] = (vcf_infos['Total_Coverage_rev'] - vcf_infos['Total_Reads_Alt_rev'])
    merged['Var_Alt'] = vcf_infos['Total_Reads_Alt']
    if have_fwd_ref:
        merged['Var_Alt_fwd'] = vcf_infos['Total_Reads_Alt_fwd']
        merged['Var_Alt_rev'] = vcf_infos['Total_Reads_Alt_rev']
    merged['Var_QualVariant'] = vcf['Qual'].astype(float)
    if 'GQ' in vcf_sample:
        merged['Var_QualSample'] = vcf_sample['GQ'].astype(float)
    else:
        merged['Var_QualSample'] = None
    merged['Comments'] = ''

    #manually add filter status for low coverage
    merged.loc[merged['Var_TotalCoverage'] < 10, 'Var_FailedFilters'] = [filt + ';CovLt10' if len(filt) > 0 else 'CovLt10' for filt in merged.loc[merged['Var_TotalCoverage'] < 10, 'Var_FailedFilters']]

    #calculate the deleterious score
    if include_score:
        calculate_del_score(merged)

    #add gbrowse and regionseq links
    if include_gbrowse_links:
        merged['GBrowse'] = ""
        merged['RegionSeq'] = ""
        if genome_name == 'hg19':
            merged['GBrowse'] = [ '=HYPERLINK("https://gbrowse2.molbiol.ox.ac.uk/fgb2/gbrowse/CRANIOFACIAL_GRCh37/?name=%s:%d-%d", "GBrowse")' % (
                row.Chr, row.Start - 25, row.End + 25
            ) for row in merged[['Chr', 'Start', 'End']].itertuples()]
            merged['RegionSeq'] = [ '=HYPERLINK("https://gbrowse2.molbiol.ox.ac.uk/cgi-bin/varSeqRegion.cgi?var_chr=%s&var_posn=%d&ref=%s&var=%s&rm=mode_2", "RegionSeq")' % (
                row.Chr, row.Start, row.Ref, row.Alt
            ) for row in merged[['Chr', 'Start', 'Ref', 'Alt']].itertuples()]
        elif genome_name == 'hg38':
            merged['GBrowse'] = [ '=HYPERLINK("https://gbrowse2.molbiol.ox.ac.uk/fgb2/gbrowse/CRANIOFACIAL_GRCh38/?name=%s:%d-%d", "GBrowse")' % (
                 row.Chr, row.Start - 25, row.End + 25
            ) for row in merged[['Chr', 'Start', 'End']].itertuples()]
        else:
            print('Not generating GBrowse/RegionSeq columns since reference build is neither hg19 nor hg38.')

    #NB: always take the first one
    if target_intervals is not None:
        merged['Target'] = [ (
            next(target_range[2].id for target_range in target_intervals[row.Chr].find( (row.Start, row.End) )) #note the double parenthesis for tuple
        ) for row in merged[['Chr', 'Start', 'End']].itertuples()]

    #join sample_info
    if sample_info is not None:
        merged = merged.join(sample_info, on = ['Sample', 'Target'], how = 'left')

    #drop useless columns
    merged.drop(['Otherinfo'], axis=1, inplace=True)
    #drop extra populations from ExAC/gnomAD (all columns that start with ExAC/gnomAD and don't end with _ALL)
    for db_prefix in ['ExAC', 'gnomAD_genome', 'gnomAD_exome']:
        merged.drop([column for column in merged.columns if column.startswith(db_prefix) and not column.endswith('_ALL')], axis=1, inplace=True)

    #first columns
    first_cols = []
    if include_score:
        first_cols.append('DeleteriousScore')
    if 'Gene.refGene' in merged.columns:
        first_cols += ['Gene.refGene', 'Func.refGene', 'ExonicFunc.refGene', 'AAChange.refGene']
    if include_exon_distance:
        first_cols.append('DistanceNearestExon')
    first_cols = first_cols + [c for c in merged.columns if c.startswith('Var_') and c != 'Var_Qual']

    #last columns
    last_cols = ['Chr', 'Start', 'End', 'Ref', 'Alt']
    if include_gbrowse_links:
        last_cols.append('GBrowse')
    last_cols += ['Sample'] + sample_info_columns + ['Comments']

    #columns to remove
    ignored_cols = ['GeneDetail.refGene']

    #select and reorder columns
    merged.sort_values(['Sample', 'Chr', 'Start'], inplace=True)     
    merged = merged[first_cols + [c for c in merged.columns if not c in first_cols + last_cols + ignored_cols] + last_cols]

    #output
    return merged

def make_summary_condensed(input, output):
    """Make condensed summary table that only contains a subset of columns."""

    #load sample information table
    sample_info_columns = []
    if len(input['sample_info']) > 0:
        sample_info = read_sample_info(input['sample_info'][0])
        sample_info_columns = list(sample_info.columns)

    for do_filter in ['filtered', 'unfiltered']:
        try:
            #load
            merged = pd.read_csv(input['summary'], index_col = False, dtype = variants_dtypes, na_values = variants_na_values)

            #filter
            if do_filter == 'filtered':
                print('Got %d rows with %d passing filters and %d coverage >= 10' % (len(merged), (merged['Var_FailedFilters'].isnull()).sum(), (merged['Var_TotalCoverage'] >= 10).sum()))
                merged = merged[ (merged['Var_FailedFilters'].isnull()) & (merged['Var_TotalCoverage'] >= 10) ]
                print('Got %d rows after filtering' % len(merged))
            else:
                print('Not filtering, got %d rows' % len(merged))

            #find column with dbsnp-id (taking first that starts with avsnp)
            avsnp_columns = [c for c in merged.columns if c.startswith('avsnp')]
            if len(avsnp_columns) > 0:
                merged['ID_dbsnp'] = merged[avsnp_columns[0]]
            #rename column
            merged['Genotype'] = merged['Var_Zygosity']
            #get desired columns
            output_columns = ([] if do_filter == 'filtered' else ['Var_FailedFilters']) + [
                    'Sample', 'Target',
                    'Chr', 'Start', 'Ref', 'Alt', 'Genotype',
                    'ID_dbsnp', 'DeleteriousScore',
                    'Gene.refGene', 'Func.refGene', 'ExonicFunc.refGene', 'AAChange.refGene'
                ] + sample_info_columns
            #make sure the columns actually exist
            output_columns = [c for c in output_columns if c in merged.columns]

            merged = merged[output_columns]

            #output
            merged.to_csv(output[do_filter], index = False)
        except pd.io.common.EmptyDataError:
            print('No variants found, creating empty file!')
            open(output[do_filter], 'a').close()

def make_summary_excel(input, output):
    """UNTESTED: make Excel table for merged table"""
    import openpyxl

    try:
        merged = pd.read_csv(input[0], index_col = False)

        #reset header format
        pd.formats.format.header_style = None

        #requires openpyxl
        with pd.ExcelWriter(output[0]) as xlsx:
            merged.to_excel(xlsx, sheet_name='All (%d)' % len(merged))

            merged['_sheet'] = 'Unknown'
            for sheet, regex in {
                'ncRNA': 'ncRNA',
                'Intronic intergenic': 'intronic|intergenic',
                'Splicing': 'splic',
                'Up-downstream': 'stream',
                'UTR': 'utr'
            }.items():
                merged.loc[merged['Func.refGene'].notnull() & merged['Func.refGene'].str.contains(regex, case=False), '_sheet'] = sheet

            for sheet, regex in {
                'Stops': 'stop',
                'Nonsyn': 'nonsynonymous',
                'Synon': 'synonymous',
                'Fshifts': 'frameshift'
            }.items():
                merged.loc[merged['ExonicFunc.refGene'].notnull() & merged['ExonicFunc.refGene'].str.contains(regex, case=False), '_sheet'] = sheet

            sheets = ['Nonsyn', 'Stops', 'Fshifts', 'Splicing', 'Synon', 'Intronic intergenic', 'Up-downstream', 'UTR', 'ncRNA', 'Unknown']
            assert merged['_sheet'].isin(sheets).all()
            for sheet in sheets:
                my_merged = merged[merged['_sheet'] == sheet]
                if len(my_merged) > 0:
                    my_merged.drop('_sheet', axis=1, inplace=True)
                    my_merged.to_excel(xlsx, sheet_name='%s (%d)' % (sheet, len(my_merged)))

            #add formatting
            workbook = xlsx.book
            for worksheet in workbook:
                #http://openpyxl.readthedocs.io/en/default/api/openpyxl.worksheet.worksheet.html
                #http://openpyxl.readthedocs.io/en/default/formatting.html
                #add colour scale on deleteriousness score
                worksheet.conditional_formatting.add('A1:A%d'%worksheet.max_row, openpyxl.formatting.rule.ColorScaleRule(start_color = 'ffffff', end_color = 'aa0000'))
                #causes error, not clear why...
                pass

            xlsx.save()
    except pd.io.common.EmptyDataError:
        print('No variants found, creating empty file!')
        open(output[0], 'a').close()
