==========================================================
amplimap documentation
==========================================================

A mapping and analysis pipeline for targeted NGS data (version |version|)
--------------------------------------------------------------------------

amplimap is a command-line tool to automate the processing and analysis of data from targeted next-generation sequencing (NGS) experiments with PCR-based amplicons or capture-based enrichment systems.

From raw sequencing reads, amplimap generates a variety of output files including read alignments, per-basepair nucleotide counts, target coverage data and annotated variant calls.

In addition to its focus on user-friendliness and reproducibility, amplimap supports advanced features such as the generation of consensus base calls for read families based on molecular identifiers/barcodes (UMIs) and the detection of chimeric reads caused by amplification of off-target loci.

Overview
----------
To run amplimap you create a directory containing a small set of input files (see :ref:`usage`):

- A subdirectory with FASTQ.GZ or BAM files representing your different samples (tested with Illumina MiSeq, HiSeq and NextSeq)

- Optionally: Files describing the targeted genomic regions, the primers you used or other custom configuration parameters

Then you can run ``amplimap`` to generate a variety of different output files, depending on your experiment.
These include, for example:

1. A target coverage table, showing you how well-covered each target region was in each sample.

2. A table of germline variants in your samples, annotated with gene, impact, population frequencies, deleteriousness scores, etc.

3. A per-basepair “pileup” table telling you how often each nucleotide was seen in each sample at each position.

Built on top of `Snakemake <https://snakemake.readthedocs.io/>`_ and Python 3, amplimap is entirely
automated and can be run on a single machine as well as on an HPC cluster
(e.g. LSF, SGE).

Pipeline diagram & cheat sheet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: amplimap-diagram.png
   :scale: 20%

:download:`Download PDF version <amplimap-diagram.pdf>`

Supported experimental protocols
---------------------------------
amplimap is compatible with most targeted sequencing protocols that generate paired-end short read data.

For protocols utilising PCR or smMIPs each read should start with a known primer (or targeting arm) sequence, followed by the amplified target DNA.
Reads can optionally contain a unique molecular identifier (UMI) sequence in front of the primer, which can be used to group reads into families.
Data should be available as demultiplexed FASTQ.GZ files, with each pair of files representing a different sample.

For capture-based protocols data can be provided in FASTQ.GZ or unmapped/mapped BAM format, which may contain UMIs as BAM tags.
See :ref:`running-capture` for details.

Some of the protocols we have analyzed with amplimap include:

- PCR-based targeted resequencing (single/multiplex)
- `smMIPs with and without UMIs <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3638140/>`_
- Probe based target enrichment, for example:

  - `IDT xGen Lockdown probes <https://www.idtdna.com/pages/products/next-generation-sequencing/hybridization-capture/custom-probes-panels/xgen-lockdown-probes>`_
  - `Twist Bioscience Custom Panels <https://twistbioscience.com/products/ngs#product-featured-2911>`_

Documentation contents
------------------------

.. toctree::
   :maxdepth: 2

   extended_installation
   quickstart
   usage
   configuration
   tutorials
   advanced
   code
   references

Links
--------

- Package: https://pypi.org/project/amplimap/
- Code: https://github.com/koelling/amplimap/
- Documentation: https://amplimap.readthedocs.io/

Citation and License
--------------------
Licensed under the Apache License, version 2.0.
Copyright 2020 Nils Koelling.
When you use amplimap,
please cite the `amplimap paper <https://academic.oup.com/bioinformatics/article/35/24/5349/5539690>`_
in your work:

   Nils Koelling, Marie Bernkopf, Eduardo Calpena, Geoffrey J Maher, Kerry A Miller, Hannah K Ralph, Anne Goriely, Andrew O M Wilkie, amplimap: a versatile tool to process and analyze targeted NGS data, Bioinformatics, Volume 35, Issue 24, 15 December 2019, Pages 5349–5350, https://doi.org/10.1093/bioinformatics/btz582
