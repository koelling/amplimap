# -*- coding: utf-8 -*-
"""
This is a custom read mapper for amplicon data.
It simply places the reads at the location they were expected to come from,
based on the target coordinates of their associated probe.
It will perform an aligment within the target region to handle small differences
but will not correctly handle chimeric off-target reads or errors.
"""

import sys, os, re, time, collections

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
log.addHandler(sh)

#for bam creation
import pysam
#for getting ref chrs 
import pyfaidx
#for fastq reading
import Bio.SeqIO
import Bio.pairwise2
import gzip

#amplimap imports
from .common import parse_extended_read_name
from .reader import read_new_probe_design
from .version import __title__, __version__

class AmplimapNoAlignment(Exception):
    pass

def align_and_find_cigar(read, ref, debug = False, reverse = False) -> tuple:
    """
    Align read and reference sequence using parwise global alignment and return start offset and CIGAR ops.
    """ 
    assert not reverse, 'not in use'

    #we always want to be left-aligned with the proper seq (although it's not clear whether this actually makes a diff)
    if reverse:
        read = read[::-1]
        ref = ref[::-1]
    alignments = Bio.pairwise2.align.globalms(read, ref, 2, -1, -3, -.5, one_alignment_only = True, penalize_end_gaps = False)
    if len(alignments) == 0:
        raise AmplimapNoAlignment()

    alignment = alignments[0]
    #now we need to reverse these again, so that the cigar string (which is in genomic orientation) makes sense
    if reverse:
        alignment = (alignment[0][::-1], alignment[1][::-1])
    return find_cigar_for_alignment(len(read), alignment, debug)

def find_cigar_for_alignment(read_len, alignment, debug) -> tuple:
    """
    Generate tuple of start offset and CIGAR operations for given alignment.
    """ 
    cigartuples = []
    cigar_ref_consumed = 0
    cigar_read_consumed = 0
    cigar_read_start_offset = 0

    if debug:
        print(alignment)

    for char_read, char_ref in zip(alignment[0], alignment[1]):
        char_cigar = None
        if char_read == '-' and char_ref == '-':
            raise Exception('both gaps?!')
        elif char_read == '-':
            #if we already aligned some bases then this is a deletion
            if cigar_read_consumed > 0:
                char_cigar = 2 #DEL
            #otherwise we have a start offset
            else:
                cigar_read_start_offset += 1
            cigar_ref_consumed += 1
        elif char_ref == '-':
            #if we already started the read this is an insertion, otherwise we soft-clip
            if cigar_ref_consumed > 0:
                char_cigar = 1 #INS
            else:
                char_cigar = 4 #SOFTCLIP
            cigar_read_consumed += 1
        else:
            char_cigar = 0 #MATCH
            cigar_read_consumed += 1
            cigar_ref_consumed += 1

        #this may be None if we have start offset
        if char_cigar is not None:
            #start new tuple or increase counter for current tuple
            if len(cigartuples) == 0 or cigartuples[-1][0] != char_cigar:
                cigartuples.append( [char_cigar, 1] )
            else:
                cigartuples[-1][1] += 1

        # if debug:
        #     print(char_read, cigar_read_consumed, char_ref, cigar_ref_consumed, char_cigar, cigar_read_start_offset)
        
        #cigar should not get longer than read
        if cigar_read_consumed == read_len:
            break

    if debug:
        print(cigar_read_start_offset, cigartuples, sum( [x[1] for x in cigartuples if x[0] in [0, 1, 4]] ))
        print(read_len)

    #ensure that the last operation isn't an insertion (it should just be soft-clipped)
    if cigartuples[-1][0] == 1:
        cigartuples[-1][0] = 4

    #make sure we consumed the entire read
    assert cigar_read_consumed == read_len
    assert sum( [x[1] for x in cigartuples if x[0] in [0, 1, 4]] ) == read_len

    #make these into actual tuples
    cigartuples = [ tuple(x) for x in cigartuples ]
            
    return cigar_read_start_offset, cigartuples

def create_bam(
    sample,
    files_in1,
    files_in2,
    ref_fasta,
    probes_dict,
    output,
    has_trimmed_primers = True,
    debug = False
):
    """
    Create a BAM file with reads placed at their expected locations, adjusted through pairwise alignment to the target sequences.

    This will give reasonable results as long as probes capture the exact target sequences, but will generate
    alignments with many mismatches if there are any discrepancies.
    """
    assert len(files_in1) == len(files_in2)

    tStart = time.time()
    counters = collections.Counter()

    ref_idx = pyfaidx.Faidx(ref_fasta, rebuild=False)
    bam_header = {
        'HD': { 'VN': '1.0' },
        'SQ': [ {'LN': record.rlen, 'SN': name} for name, record in ref_idx.index.items() ],
        'RG': [ { 'ID': sample, 'SM': sample } ],
        'PG': [ { 'ID': __title__, 'PN': __title__, 'VN': __version__ } ],
    }

    chr_indices = { chrom: index for index, chrom in enumerate(ref_idx.index.keys()) }

    with pysam.AlignmentFile(output, "wb", header = bam_header) as pairedreads:
        for ixfile in range(len(files_in1)):
            file1 = files_in1[ixfile]
            file2 = files_in2[ixfile]
            log.info('Processing %s and %s (#%d)', file1, file2, ixfile)

            counters['files'] += 1
            opener = gzip.open if file1.endswith('.gz') else open
            with opener(file1, 'rt') as hdl1, opener(file2, 'rt') as hdl2:
                for read_pair in zip( Bio.SeqIO.QualityIO.FastqGeneralIterator(hdl1),
                        Bio.SeqIO.QualityIO.FastqGeneralIterator(hdl2) ):
                    counters['pairs_total'] += 1
                    if counters['pairs_total'] % 50000 == 0:
                        log.info("processed %d pairs - %.f sec elapsed, %.4f sec/pair, %.1f pairs/sec",
                            counters['pairs_total'],
                            time.time() - tStart,
                            (time.time() - tStart) / counters['pairs_total'],
                            counters['pairs_total'] / (time.time() - tStart))

                    if debug and counters['pairs_total'] > 10:
                        print('DEBUG - stopping after ', counters['pairs_total'])
                        break

                    #extract and parse read name
                    read_names_original = [
                        read_pair[0][0].split('\t')[0],
                        read_pair[1][0].split('\t')[0],
                    ]
                    assert len(read_names_original[0]) > 0
                    assert read_names_original[0] == read_names_original[1]
                    read_name, read_probe, read_umi = parse_extended_read_name(read_names_original[0])

                    probe_chr = probes_dict['chr'][read_probe]
                    if not probe_chr in chr_indices:
                        raise Exception('Probe {} is associated with chromosome {}, but this entry does not exist in the reference fasta file!'.format(
                            read_probe, probe_chr
                        ))
                    probe_chr_index = chr_indices[probe_chr]

                    read_lens = [ len(read_pair[read_number][1]) for read_number in range(2) ]

                    #untested: if we haven't trimmed off the primers then we need to start aligning from the primer start location!
                    if has_trimmed_primers:
                        probe_start = probes_dict['target_start_0'][read_probe] + 1
                        probe_end = probes_dict['target_end'][read_probe]
                    else:
                        probe_start = probes_dict['probe_start_0'][read_probe] + 1
                        probe_end = probes_dict['probe_end'][read_probe]

                    if probes_dict['strand'][read_probe] == '+':
                        read_starts = [ probe_start, probe_end - read_lens[1] + 1 ]
                        read_reverse = [ False, True ]
                    elif probes_dict['strand'][read_probe] == '-':
                        read_starts = [ probe_end - read_lens[0] + 1, probe_start ]
                        read_reverse = [ True, False ]
                    else: raise Exception('Unexpected strand for probe {}'.format(read_probe))

                    #NOTE: this SHOULD BE one-based based on documentation
                    #but actually seem to be ZERO-based -- at least the sequence we get for PRRX1-Ex1
                    #starts CGGA but should start GGA; ends TTC but should end TTCT if we just use probe_start and probe_end
                    #they are always in genomic sense
                    probe_target_sequence = str(ref_idx.fetch(probe_chr, probe_start, probe_end)).upper()

                    #sanity check that we got the right sequence
                    if has_trimmed_primers:
                        assert len(probe_target_sequence) == probes_dict['target_length'][read_probe]
                    else:
                        assert len(probe_target_sequence) == probes_dict['capture_size'][read_probe]                        

                    if debug:
                        print(read_name, read_probe, read_umi)
                        print(probe_chr, probe_chr_index, probe_start, probe_end)
                        print(read_starts)
                        print(read_reverse)                       
                        print(probe_target_sequence)               
                        print(read_pair)

                    try:
                        #pre-process alignments to make sure the mate starts are actually correct
                        read_cigars = []
                        read_sequences = []
                        read_tags_for_pysam = []
                        for read_number in range(2):
                            #copy over our custom tags from FASTQ file
                            read_tags = [ tag.split(':') for tag in read_pair[read_number][0].split('\t')[1:] ]
                            read_tags_for_pysam.append (
                                [("RG", sample, "Z")] + [
                                    (tag_name, int(tag_value) if tag_type == 'i' else tag_value, tag_type) for tag_name, tag_type, tag_value in read_tags
                                ]
                            )
                            if debug:
                                print(read_tags)
                                print(read_tags_for_pysam[read_number])

                            #figure out sequence
                            read_sequence = str(Bio.Seq.Seq(read_pair[read_number][1]).reverse_complement()) if read_reverse[read_number] else read_pair[read_number][1]
                            read_sequences.append(read_sequence)

                            #align read to target sequence -- note both of these are in genomic sense!
                            try:
                                cigar_read_start_offset, cigartuples = align_and_find_cigar(read_sequence, probe_target_sequence)
                                read_tags_for_pysam[read_number].append( ("so", cigar_read_start_offset, 'i') )

                                #remember cigar
                                read_cigars.append(cigartuples)
                                #adjust start -- need to use zero-based coords here but probe_start is 1-based
                                read_starts[read_number] = probe_start - 1 + cigar_read_start_offset
                            except AssertionError:
                                cigar_read_start_offset, cigartuples = align_and_find_cigar(read_sequence, probe_target_sequence,
                                    debug = True)
                                raise

                        if debug:
                            print(read_cigars)
                            print(read_starts)

                        for read_number in range(2):

                            #create aligned segment
                            a = pysam.AlignedSegment()
                            a.mapping_quality = 255 #always best quality
                            a.query_name = read_names_original[0]
                            a.query_sequence = read_sequences[read_number]
                            a.query_qualities = pysam.qualitystring_to_array(read_pair[read_number][2][::-1] if read_reverse[read_number] else read_pair[read_number][2])
                            a.set_tags( read_tags_for_pysam[read_number] )
                            a.cigartuples = read_cigars[read_number]

                            a.reference_id = probe_chr_index
                            a.reference_start = read_starts[read_number]
                            a.next_reference_id = probe_chr_index
                            a.next_reference_start = read_starts[1 - read_number]
                            #a.template_length = read_lens[read_number]

                            a.is_paired = True
                            a.is_proper_pair = True
                            a.is_read1 = read_number == 0
                            a.is_read2 = read_number == 1
                            a.is_reverse = read_reverse[read_number]
                            a.mate_is_reverse = read_reverse[1-read_number]

                            if debug:
                                print(a)
                            pairedreads.write(a)
                            if debug:
                                break
                    #normally we always get an alignment, but apparently sometimes we don't?
                    except AmplimapNoAlignment:
                        counters['no_alignment'] += 1
                        pass


    log.info('%s done - %d pairs in total, %d without alignment', sample, counters['pairs_total'], counters['no_alignment'] )

    log.info("BAM file created: %s", output)