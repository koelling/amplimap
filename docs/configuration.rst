.. _configuration:

Configuration
-------------

.. _default-config:

Default configuration
~~~~~~~~~~~~~~~~~~~~~

The default configuration file is called :file:`config_default.yaml`. It is
located in the amplimap basedir, which is usually the
:file:`amplimap.{VERSION}.egg/amplimap` directory located in Python’s
:file:`site-packages`
(where *VERSION* is the amplimap version number, e.g. |version|).
You can run ``amplimap --basedir`` to get the path to the basedir.

Any settings in this file will be applied every time you run amplimap.
This is particularly helpful for setting up correct paths for the reference
files (genome build, aligner index, reference genome fasta).

You can also save this file under :file:`/etc/amplimap/{VERSION}/config.yaml`
(where *VERSION* is the amplimap version number, e.g. |version|) or provide a different
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

.. _config-tools:

Selecting the aligner and variant caller
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

amplimap can work with different aligners and variant callers.

Supported aligners, specified through ``align: aligner:``, are:

  - BWA (``bwa``)
  - Bowtie2 (``bowtie2``)
  - STAR (``star``)

Supported variant callers, specified through ``variants: caller:``, are:

  - Platypus (``platypus``)
  - GATK 4 (``gatk``)
  - weCall (``wecall``, experimental)
  - Octopus (``octopus``, experimental)

For example:

::

    align:
      aligner: "bowtie2"
    variants:
      caller: "octopus"

Additional aligners and variant callers can also be added by specifying the relevant commands under `tools:`.
See the comments in the config file for details.

.. _config-reference:

Reference genome paths
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

amplimap requires a reference genome and associated indices (such as a FASTA index
or a bwa index) to run.
It is recommended that you specify these paths in the :ref:`default-config` file.
For example, to specify paths for ``hg19`` and ``hg38``
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

For suggestions on where to obtain these references and how to create indices see
:ref:`installation-setup`.

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

.. _config-annovar:

Setting up Annovar
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Annovar is the software amplimap uses to annotate variant calls.
For licensing reasons it needs to be downloaded and installed manually.
Please see the `Annovar website <http://annovar.openbioinformatics.org/en/latest/user-guide/download/>`_
for details.

Once you have Annovar installed we recommend that you download the following indices:

- refGene, esp6500siv2_all, 1000g2014oct_all, avsnp147, cosmic82, dbnsfp33a, clinvar_20150629, dbscsnv11, exac03, gnomad_genome, gnomad_exome

Finally, you need to specify the path to your Annovar index directory in your :ref:`default-config` file (see :ref:`config-reference`).
If you downloaded a different set of indices you also need to adjust the annovar protocols and operations parameters.

Running with UMIs (e.g. for smMIPs)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If one or both of your reads start with UMIs, you have to specify their lengths
in the configuration file using the ``umi_one:`` and ``umi_two:`` settings
under ``parse_reads:``.

For example, to process an experiment with 5bp UMIs on each read, your
:file:`config.yaml` might look like this:


::

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
primer trimming, specify ``trim_primers: false`` under ``parse_reads:``.

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



All configuration options
~~~~~~~~~~~~~~~~~~~~~~~~~
A commented list of all configuration options supported by amplimap and their default values
is available in config_default.yaml:

.. include:: ../amplimap/config_default.yaml
  :start-line: 4
  :code: yaml
