.. image:: https://raw.githubusercontent.com/koelling/amplimap/master/amplimap_logo_400px.png
	:width: 400px

==========================================================
amplimap: amplicon mapping and analysis pipeline
==========================================================


amplimap is a pipeline to process and analyse short-read sequencing
data from **smMIP-panels** or **PCR-based amplicons**. It was designed
to support both **germline variant calling** as well as **quantification of
variant allele frequencies** using ultra-high-coverage pileups with and without
**UMIs**.

amplimap takes fastq or bam files generated from an Illumina
sequencer (tested with MiSeq, HiSeq, and NextSeq), assigns each read pair to a
probe/amplicon, trims and aligns the reads and then generates a set of
basic statistics. After that, two different types of analyses can be
performed:

1. Germline variant calling using Platypus/GATK and annotation with Annovar,
   to generate an annotated table of potentially
   pathogenic variants.

2. Per-basepair “pileup” of reads to generate nucleotide counts of
   each base at each target basepair and each sample, for analysis of
   allele frequencies.

Built on top of Snakemake and Python 3, amplimap is entirely
automated and can be run on a single machine as well as on a HPC cluster
(eg. LSF, SGE).

Links
--------

- Package: https://pypi.org/project/amplimap/
- Code: https://github.com/koelling/amplimap/
- Documentation: https://amplimap.readthedocs.io/


Basic installation
-------------------
::

	pip install amplimap

Requires Python 3.5+