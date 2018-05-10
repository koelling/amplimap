Configuration
-------------

.. _default-config:

Default configuration
~~~~~~~~~~~~~~~~~~~~~

The default configuration file is called :file:`config_default.yaml`. It is
located in the amplimap basedir, which is usually the
:file:`amplimap.{VERSION}.egg/amplimap` directory located in Python’s
:file:`site-packages`
(where *VERSION* is the amplimap version number, eg. |version|).
You can run ``amplimap --basedir`` to get the path to the basedir.

Any settings in this file will be applied every time you run amplimap.
This is particularly helpful for setting up correct paths for the reference
files (genome build, aligner index, reference genome fasta).

You can also save this file under :file:`/etc/amplimap/{VERSION}/config.yaml`
(where *VERSION* is the amplimap version number, eg. |version|) or provide a different
path in the  :envvar:`AMPLIMAP_CONFIG` environment variable.

Local configuration
~~~~~~~~~~~~~~~~~~~~~

To specify experiment-specific settings, you can place a file called :file:`config.yaml` in your working
directory. Any setting that is specified in this local configuration
file will override the default configuration. This is useful for setting
analysis-specific parameters, such as the quality filters, UMI lengths,
etc.

To see the configuration that amplimap will use, based on your global and local
configuration files, run ``amplimap --print-config``.

Common configuration changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reference genome paths
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

amplimap requires a reference genome and associated indices (such as a FASTA index
or a bwa index) to run. For example, to specify paths for ``hg19`` and ``hg38``
of the human genome and set the default to ``hg38``:

::

    paths:
      hg38:
        bwa: "/PATH/TO/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        fasta: "/PATH/TO/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        annovar: '/PATH/TO/annovar/humandb'
      hg19:
        bwa: "/PATH/TO/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
        fasta: "/PATH/TO/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
        annovar: '/PATH/TO/annovar/humandb'
    general:
      genome_name: "hg38"


It is recommended that you specify these paths in the :ref:`default-config` file.

Once you have set up these paths, you can then choose the genome to use for each
experiment by specifying the ``genome_name`` in your local configuration file:

::

    general:
      genome_name: "hg38"

However, you can also specify the paths and genome in the same file.
Add a new section with a name of your choice under ``paths:`` and then
set your ``genome_name`` to the same name. For example:

::

    paths:
      mm10:
        bwa: "mm10/bwa"
        fasta: "mm10/genome.fa"
    general:
      genome_name: "mm10"

Note that when you are doing variant annotation with Annovar your
``genome_name`` has to match the name that Annovar uses.

Running with UMIs (eg. for smMIPs)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If one or both of your reads start with UMIs, you have to specify their lengths
in the configuration file using the ``umi_one:`` and ``umi_two:`` settings
under ``parse_reads:``. In addition, you probably want to set ``ignore_umis: false``
under ``general:`` to enable UMI grouping for pileups and alignment stats.

For example, to process an experiment with 5bp UMIs on each read, your
:file:`config.yaml` might look like this:


::

    general:
      ignore_umis: false
    parse_reads:
      umi_one: 5
      umi_two: 5

Note that it is very important to specify the correct lengths here, since
these UMIs will be trimmed off before amplimap tries to match the start of the read
to the expected primer sequence. If the length is incorrect, the primer sequences
will never match the reads and all of the reads will be discarded.

Primer trimming
^^^^^^^^^^^^^^^^^^^^^^

By default, primer (extension/ligation) arms are removed from the
beginnings and, if applicable, ends of reads before alignment. This is
particularly important when using overlapping (tiled) probes, since the
primers would otherwise skew the observed allele frequencies or even
prevent a variant from being called in the first place. They can also
lead to misalignment of off-target sequences that were inadvertendly
captured, introducing false positives. However, removing them also means
that only the targeted region in-between the arms will be aligned to the
genome. This can be problematic if its sequence is not unique, leading
to off-target alignment and reads with mapping quality 0. To turn off
primer trimming, specify ``trim_primers: false`` under ``general:``.

Quality trimming of reads
^^^^^^^^^^^^^^^^^^^^^^^^^^

Reads can optionally be trimmed at their beginnings/ends to remove
low-quality bases. This may be helpful to remove potentially noisy base
calls during variant calling, although most variant callers should be
able to account for this independently. To enable this, set a quality
trimming threshold, which is the highest probability of an errorneous
call that you would like to allow. The default (which results in quality
trimming being turned off) is ``false``, a suggested value to enable
quality trimming would be 0.01 (1%): ``quality_trim_threshold: 0.01``.

Minimum mapping quality (for pileups only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, no mapping quality filter is applied for the pileup and
alignment stats tables. If you think that filtering out low-quality
mappings may improve your results, you can change this by setting a
minimum mapping quality in the ``pileup:`` section using something like
``min_mapq: 20``. Note that this setting has no effect on coverage and
standard variant calling!

Support for modules
~~~~~~~~~~~~~~~~~~~~~

amplimap has some basic support for loading and unloading optional software packages
through the modules system. To use this feature, specify the modules that should be loaded
for each of the software packages listed under ``modules:``.
If you leave a setting empty, no module will be loaded and the software will have to be
available without loading a module.