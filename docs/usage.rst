Input and the working directory
-------------------------------
For each experiment (for example each sequencing run), you need to create a working directory
inside which you run amplimap. All the input and output files will be kept in this directory.

Input overview
~~~~~~~~~~~~~~

Required input
^^^^^^^^^^^^^^^^

-  :ref:`reads-in`: subdirectory containing pairs of gzipped FASTQ files

   -  Alternative: :ref:`unmapped-bams`

-  :ref:`probes-csv`: probe design table in CSV format

   -  Alternative: :ref:`probes-mipgen-csv`

Optional input
^^^^^^^^^^^^^^^^

-  :ref:`targets-csv`: target regions in CSV format 
  
  - Alternative: :ref:`targets-bed`

-  ``config.yaml``: configuration file (see :doc:`configuration`)
-  :ref:`snps-txt`: list of known SNPs to genotype
-  :ref:`sample-info-csv`: sample annotation for variant and coverage
   tables

Input files
~~~~~~~~~~~~~~

.. _reads-in:

reads_in/
^^^^^^^^^^^^^^^^^^^^^^^^

This subdirectory should contain pairs of gzipped FASTQ files from
paired-end sequencing.
Filenames should follow standard Illumina naming conventions
(ending in ``_L001_R1_001.fastq.gz`` and
``_L001_R2_001.fastq.gz``, where ``L001`` is the lane).
Data from multiple lanes (eg. ``L001`` and ``L002``) will be
merged, as long as all samples have a pair of files for each lane.

.. _unmapped-bams:

Unmapped BAM files (unmapped_bams_in/)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Instead of the reads_in directory you can also provide your sequencing data
as unmapped BAM files inside a 
directory called ``unmapped_bams_in/``.
If your reads contained UMIs, these should have been trimmed off already and be
provided in the bam file using the ``RX`` tag.

This is useful for running amplimap on data from capture-based protocols.
See also :ref:`running-capture`.

If you run amplimap with unmapped BAM files, probes can no longer be identified
based on the primer arms. Thus, some of the output files can not be generated.

.. _probes-csv:

probes.csv
^^^^^^^^^^^^^^^^^^^^^^^^

This file contains information about the sequences of the primers and
their expected location on the reference genome. It can be generated
directly from the output table created by *MIPGEN* - see below!

The columns are:

1. ``id``: a unique name for the probe
2. ``first_primer_5to3``: sequence of the forward primer, in 5’ to 3’
   orientation (needs to match start of read 1)
3. ``second_primer_5to3``: sequence of the reverse primer, in 5’ to 3’
   orientation (needs to match start of read 2)
4. ``chr``: chromosome (should start with *chr*)
5. ``target_start``: start of the target region (after the primer,
   1-based coordinate system)
6. ``target_end``: end of the target region
7. ``strand``: ``+`` or ``-``

The first row of the file should be the header listing these column names. Any
additional columns are ignored.

No two probe arms should be so similar to each other that a read could
match both of them. Otherwise, the read will be ambigous and may end up
being counted multiple times, once for each probe. This will skew any
downstream analysis.

If you have two probes for the same region to account for a SNP in the
arm sequence, you need to provide a single merged entry for this. To
create such a merged probe, simply replace all ambiguous (SNP)
nucleotides in the arm sequences with a dot (``.``). This way, any
nucleotide will be allowed in this location and reads from either
version of the probe will be counted together.

MIP names cannot contain characters other than alphanumeric characters
(``A-Z``, ``0-9``), or ``_:+-``. Avoid using multiple colons in a row
(eg. ``::``) since this is used as a field separator internally.

The file needs to be in plain CSV format with UNIX/Windows (not Mac)
style line endings.

.. _probes-mipgen-csv:

MIPGEN probe design table (probes_mipgen.csv)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Instead of the standard ``probes.csv`` file you can also provide a
*MIPGEN* probe design table. Simply save it in CSV format using the
filename ``probes_mipgen.csv``. When you run amplimap, this file will
automatically be converted into a ``probes.csv`` file with the right
format.

If *MIPGEN* generated two versions of the same probe to account for a
SNP, ``amplimap`` will detect this based on the identical value in the
``mip_name`` column and merge them into a single line, replacing any
differences in the primer sequences by a dot (see above). Any duplicate
probe names that differ in their location, or by more than 10 characters
will cause an error.

.. _targets-csv:

targets.csv
^^^^^^^^^^^^^^^^^^^^^^^^

List of target regions (eg. exons, not the MIPs themselves) in CSV format.
This file should contain the following columns:

1. ``chr`` (chromosome, should start with *chr*)
2. ``start`` (start position, 1-based coordinate system)
3. ``end`` (end position)
4. [optional] ``id`` (name of the target)

The first row of the file should be the header listing these column names. Any
additional columns are ignored.

Variants will only be called inside these target regions! If any of the
target regions overlap, they will be merged for variant calling and
cause an error when trying to calculate pileups.

.. _targets-bed:

targets.bed
''''''''''''''''''''''''''''''''''''''''''''''''''''''''

You can also provide this data in BED format. In that case, the file should be called
``targets.bed`` and use the standard BED columns (chromosome, 0-based start position,
end position, id). The score and strand columns may be included, but do not
have any effect on the pipeline.
Note that BED files do *not* contain column headers!

.. _snps-txt:

snps.txt (for allele counting pileup)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have certain SNPs that you want to generate pileups for, you can
provide a list in tab-separated text format here. The columns are:

1. ``chr`` (should start with *chr*)
2. ``pos`` (1-based coordinate system)
3. ``id``
4. ``snp_ref``
5. ``snp_alt`` (only a single alt allele is supported)

This will generate the ``pileup_snps`` directory with
reference/alternate allele counts for each SNP. The filter column in the
pileup tables will reflect whether the observed alleles matched the SNP
alleles, or whether additional alleles were found.

Note: This file do *not* contain column headers!


.. _sample-info-csv:

sample_info.csv (for variant calling)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file can be provided to add additional sample information columns
to the coverage and variant tables. If provided, it always needs to
start with these two columns:

1. ``Sample``: the sample id, including the ``_S123`` part

   -  needs to match the sample identifiers of the input fastq files
   -  example: ``Barcode-001_S1``

2. ``Targets``: a semicolon-separated list of target ids

   -  needs to match the ids provided in ``targets.csv``
   -  example: ``GENE1-Ex1.1;GENE1-Ex1.2;GENE1-Ex2;GENE3``

These should then be followed by one or more annotation columns, which
can contain information like the id of the corresponding individual or
other information about the samples. All of these columns will be copied
into the coverage and variant tables.

A single sample id (= barcode) can have multiple rows with different
annotation columns, as long as none of the targets are the same. In
other words, any combination of sample/target id may only occur once.

If there are two overlapping target regions and a variant call is made
in the overlapping part, it can get assigned to either of them. To avoid
errors due to this, overlapping target regions must always be listed in
pairs and never be split up. For example, if the targets ``GENE1-Ex1a``
and ``GENE1-Ex1b`` overlap, you should never have a row where you only
list ``GENE1-Ex1a`` or only list ``GENE1-Ex1b``. They should always be
listed together (``GENE1-Ex1a;GENE1-Ex1b``) or not at all.

This file needs to be in plain CSV format with UNIX/Windows (not Mac)
style line endings.


Running amplimap
----------------

The pipeline is based on Snakemake, which uses predefined rules to
figure out what it needs to do to generate a certain output file. For
example, if you wanted a file called
``analysis/variants_umi/variants_summary.csv``, the pipeline will
automatically go through the rules to figure out what commands it needs
to run to generate this file from your fastq files.

To run the pipeline, just enter ``amplimap`` in your terminal:

::

    amplimap

By default, this will only start a so-called dryrun. This will not
actually run any of the code yet, but will make sure that the expected
input files are present and tell you which jobs it would be running.

If the output of this dryrun looks as expected, you can start the actual
pipeline by adding the ``--run`` parameter:

::

    amplimap --run

This will go through the first few steps of the pipeline but will not
run the more advanced analysis-specific parts.

To run these additional steps, you need to specify so-called *target
rules* at the end of the command line. Some of these are:

``bams`` (just do the alignment and create the BAM files)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    amplimap bams

``variants`` and ``variants_umi`` (germline variant calling/annotation analysis)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To call variants from raw reads:

::

    amplimap variants

To call variants from UMI-deduplicated (but not consensus) reads:

::

    amplimap variants_umi

``coverages`` (make target coverage tables)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    amplimap coverages

``pileups`` (do low frequency analysis)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    amplimap pileups

``variants_low_frequency`` (low-frequency/somatic variant calling/annotation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[EXPERIMENTAL!] To call low-frequency variants using Mutect2 use this command:

::

    amplimap variants_low_frequency

This function is still experimental and has not been thoroughly tested.


Multiple targets
~~~~~~~~~~~~~~~~

You can also group together multiple *target rules*:

::

    amplimap variants coverages

Running on a cluster
~~~~~~~~~~~~~~~~~~~~~~

You can specify the additional parameter ``--cluster=qsub`` to run jobs
in parallel on a SGE cluster:

::

    amplimap --cluster=qsub

    amplimap --cluster=qsub variants

    amplimap --cluster=qsub pileups

This can speed up the processing by an order of magnitude, as commands
will be run in parallel instead of sequentially. However, this process
is a bit more complex and may lead to unexpected errors. If you get an
error message, try running the standard command without the
``--cluster`` parameter instead.

You can set the number of jobs to submit by setting the ``--njobs``
parameter:

::

    amplimap --cluster=qsub --njobs=5

To use other cluster environments (such as LSF), add an entry with the submission command
to the ``clusters:`` section of the config file.

Output: the ``analysis`` directory
----------------------------------

All analysis results will be written to the subdirectory ``analysis`` inside
the working directory. These include:

-  ``reads_parsed/stats_samples.csv``: Sample statistics - number of
   read pairs matching expected arms per sample, etc.
-  ``reads_parsed/stats_reads.csv``: Read statistics - number of reads
   per probe per sample, number of UMIs per probe per sample, etc.
-  ``bam/``: BAM files with aligned reads
-  ``stats_alignment/stats_alignment.csv``: Alignment statistics for
   each sample - number of read pairs and unique UMI groups aligning in
   the expected location, etc.
-  ``reads_parsed/``: unknown arm files - sequences from the start of reads that didn’t
   match any of the expected primer sequences. Will only include data for the
   first 10,000 read pairs with unknown arms.

In addition, the ``analysis`` directory will contain
``config_used.yaml``, which is a copy of the configuration that was used
at the time the pipeline was first run. Note that this will not be
updated if you run the pipeline a second time, unless you delete the old
copy first.

Target-specific output
~~~~~~~~~~~~~~~~~~~~~~

The ``analysis`` directory will contain further subdirectories for the
different analyses that were performed by the pipeline:

Germline variant calling and annotation analysis: ``coverages`` and ``variants``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  ``umi_dedup/``: deduplicated BAM files, with one read pair per UMI
   group chosen at *random* (no consensus calling, no minimum coverage
   per UMI)
-  ``bam/coverage/``: min/average/zero coverage of target regions based
   on raw reads (ignoring mapping quality)
-  ``umi_dedup/coverage/``: min/average/zero coverage of target regions
   after UMI deduplication (ignoring mapping quality)
-  ``variants_raw/``: variant calls from Platypus/GATK with ANNOVAR
   annotation based on raw reads (ignoring UMIs)

   -  full summary table: ``variants_raw/variants_summary.csv`` includes
      summary of all variants in all samples, with deleteriousness
      score, etc.
   -  filtered summary table:
      ``variants_raw/variants_summary_filtered.csv`` all variants from
      summary table that pass Platypus quality filters and have a
      coverage of at least 10

-  ``variants_umi/``: variant calls from Platypus with ANNOVAR
   annotation based on UMI deduplicated reads (requires UMIs - otherwise
   as above)

Low-frequency variation analysis: ``pileups``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  ``pileup/``: target region pileup tables

   -  per-basepair pileups based on UMI groups:
      ``pileup/pileups_long.csv``
   -  coverage of UMI groups over target regions:
      ``pileup/target_coverage.csv``

-  ``pileup_snps/``: SNP pileup tables (optional, requires ``snps.txt``)

   -  per-SNP pileups based on UMI groups:
      ``pileup_snps/target_snps_pileups_long.csv``

Additional output
~~~~~~~~~~~~~~~~~~~~~~

In addition to the ``analysis`` directory, these folders may be created:

-  ``cluster_logs/``: directory with log files for each job submitted to
   the cluster (contain error messages if cluster submission fails)