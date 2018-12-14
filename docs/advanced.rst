Advanced usage
---------------

.. _running-capture:

Ignoring read families (UMI groups)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Even when your reads contain unique molecular identifiers (UMIs) you
may want to ignore them to perform a pileup on the raw reads.
To do this, set ``ignore_umis: true`` under ``general:``

::

    general:
      ignore_umis: true


Running on capture-based data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Capture-based sequencing results in reads of varying length from the
target regions, which do not contain specific primer sequences. 
This means that some features of amplimap, such as the trimming of primer sequences
or low quality bases and the detection of off-target capture events will no longer be
applicable.

However, amplimap can still produce pileups and variant calls
from the raw FASTQ or unmapped/mapped BAM files.

To do this, create a working directory containing a
``targets.bed`` or ``targets.csv`` file as usual, as well as a
``config.yaml`` file containing at least this setting:

::

    general:
      use_raw_reads: true

If you have gzipped FASTQ files without UMIs, put them into the usual :ref:`reads-in` directory.

If you have BAM files with unmapped reads, do not create a ``reads_in`` directory
but instead create a subdirectory called ``unmapped_bams_in`` inside your
working directory. Place the files
there (see also :ref:`unmapped-bams`).

If you have BAM files with reads that have already been mapped to the genome,
you can put these into a directory called ``mapped_bams_in`` instead.

Then you can amplimap as usual to generate coverage data, pileups or variants calls.

If your reads have UMIs, these should be provided in the BAM file
as a tag, the name of which should be given in the config file:

::

    general:
      umi_tag_name: "RX"
      use_raw_reads: true

UMI-deduplicated variant calling (``variants_umi``) is currently
not supported for capture-based data.

Simulating variation
~~~~~~~~~~~~~~~~~~~~~~~~

In order to test that the pipeline works as expected it may be helpful
to simulate variants and check whether these are correctly reported in
the output files. This can be particularly helpful in noisy regions,
where an abundance of errors may mask real variants. To make this sort
of analysis easier, the pipeline comes with a simple testing framework,
which can generate variant reads based on an existing data set.

Simulation of variant reads is based on a very straightforward approach:
we search for a given **search sequence** in each read pair and replace
it with a **replacement sequence**. So to simulate a A>G SNP in the
sequence “AAAAAAAAAAA”, we would simply ask the testing framework to
search for “AAAAAAAAAAA” and replace it with “AAAAAGAAAAA”. The search
sequence should be chosen such that it is unique in the targeted
regions, but not too long to avoid problems with sequencing errors (see
Caveats).

The replacement is done on both the forward and reverse strand and in
both mates of the read pairs. The pipeline is then run using these
edited reads and the usual output files are generated.

To simulate low-frequency or heterozygous variants, the **replacement
percentage** of reads that should be edited can be specified as any
number between 0 to 100. When choosing the reads to edit, the simulation
algorithm takes into account UMIs such that read pairs with the same UMI
will be treated the same way.

Running a simulation
^^^^^^^^^^^^^^^^^^^^^

To run a simulation of variant calls, use the following command, where
SEARCH is the **search sequence**, REPLACE is the **replacement
sequence** and PERCENTAGE is the **replacement percentage**:

::

    amplimap test__SEARCH_REPLACE_PERCENTAGE/test_variants.done

For example:

::

    amplimap test__TTTACCTCTATTGTTGG_TTTACCTATATTGTTGG_50/test_variants.done

Similarly, pileups can be simulated like this:

::

    amplimap test__SEARCH_REPLACE_PERCENTAGE/test_pileups.done

For example:

::

    amplimap test__TTTACCTCTATTGTTGG_TTTACCTATATTGTTGG_0.5/test_pileups.done

Simulation statistics
^^^^^^^^^^^^^^^^^^^^^

Some statistics on the number of found occurrences of the search string
and the number of replacements that have actually been carried out can
be found in the file
``test__SEARCH_REPLACE_PERCENTAGE/stats_replacements/stats_replacements.csv``.
Note that the percentage of replacements may be significantly different
from the percentage that had been specified, especially if the read
count or UMI count is low.

Caveats
^^^^^^^^^^^^^^^^^^^^^

This testing framework is meant to be a useful first step towards
simulating data, but it is not as sophisticated as other alternatives.
In particular, there are several caveats that should be kept in mind
when using this:

-  Only exact matches of the search sequence will be replaced, so any
   reads that contain mismatches in the region will be ignored. This
   means that random sequencing errors in the read will prevent
   replacement, potentially leading to a lower variant read percentage
   than specified.
-  Only full matches of the search sequence will be replaced, so search
   sequences that are only partially covered by a read will never be
   edited. Thus, the locations for simulations should be chosen such
   that they are fully contained in all overlapping reads.
-  Matches of the sequence will be replaced regardless of the genomic
   location. Consequently, if the chosen sequence is not unique,
   multiple variants may be introduced.
-  Matches inside the sequences will be replaced as well. This may cause
   problems with matching primer sequences to expected probe arms.


Merging runs
~~~~~~~~~~~~

To merge data from multiple runs together, use the ``amplimap_merge``
script. You can run ``amplimap_merge --help`` to see the parameters.
Here is an example:

::

    amplimap_merge /data/OUTPUT_FOLDER /data/working_directory1/analysis /data/working_directory2/analysis /data/working_directory3/analysis

This will merge the variant summary and coverage files from
``/data/working_directory1``, ``2`` and ``3`` together and save them in
a folder called ``/data/OUTPUT_FOLDER``. If you only want to get one row
per sample, you can use the ``--unique-sample-id-column`` to specify the
column name containing the sample ID (eg. ``DNAId``). This will generate
an additional file called ``variants_summary_filtered.unique.csv``,
which contains all unique filtered variants, and another file called
``coverage_full.unique.csv``, which contains the highest coverage numbers
observed for each sample.

For example:

::

    amplimap_merge --unique-sample-id-column=DNAId /data/OUTPUT_FOLDER /data/working_directory1/analysis /data/working_directory2/analysis /data/working_directory3/analysis



Additional Notes
~~~~~~~~~~~~~~~~~~

Platypus variant filters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The filters that a variant may have failed are described here:
http://www.well.ox.ac.uk/~gerton/Platypus/ng.3036-S1.pdf

Using ``screen``
^^^^^^^^^^^^^^^^^^^^^

While the pipeline is running, you normally need to keep your SSH
terminal connected. When the connection is lost, the pipeline run will
be aborted.

However, you can use the ``screen`` tool to make it sure it keeps
running even when you are not connected. To do this, run the command
``screen`` in the terminal. This will start a ``screen`` session, inside
which you can now run any normal commands. Even if you now disconnect
your SSH session, any commands that are running inside ``screen`` will
continue to run. To reconnect to the ``screen`` session later and check
the status of the pipeline, connect to the same server and type
``screen -r`` (r = reattach).

To scroll up and down in ``screen`` you need to use a special key
combination: Press ``Ctrl``-``A``, and then the ``ESC`` key to activate
copy mode. In copy mode, you can use the arrow keys or ``Ctrl``-``U`` to
go up and ``Ctrl``-``D`` to go down, as well as ``?`` and ``/`` to
search backwards/forwards. Press ``ESC`` again to get back to normal
typing mode.

Linking files
^^^^^^^^^^^^^^

Instead of copying large amounts of data into the working directory you
can also just create a link from the working directory to the actual
location of the files. This way, only one copy of the files is kept on
the file system.

This is particularly useful if you make multiple working directories for
the same set of samples, to analyse them with different parameters.

To create a link, use the ``ln -s`` command in the terminal, like this:

::

    ln -s /path/to/source/location name_of_link

So for example, to link the ``probes.csv`` file from another directory
into the current directory with the same name, you can run:

::

    ln -s /other/directory/probes.csv probes.csv

You can also link multiple files using wildcards - for example, to link
all fastq.gz files from your data directory into the ``reads_in``
folder:

::

    ln -s /path/to/data/directory/*.fastq.gz reads_in/