.. _installation:

Installation through pip
------------------------
If you already have all of the required external software available (see below)
you can install amplimap directly through pip.
Please note that this **requires Python 3.5 or 3.6** and does not currently
work with Python 3.7 due to problems with the pysam package. It also
does not work with any Python version lower than 3.5.

If you do not have the dependencies and the right version of Python available
please see our :doc:`extended_installation` to install amplimap through Miniconda.

::

  # you may need to use `pip` instead of `pip3`
  pip3 install amplimap

If this does not work, you can try to install it manually:

::

  # install required python3 packages
  # you may need to use `pip` instead of `pip3`
  pip3 install setuptools Cython numpy

  # download and install amplimap
  # you may need to use `python` instead of `python3`
  git clone --depth=1 https://github.com/koelling/amplimap.git
  cd amplimap
  python3 setup.py install


You can also :download:`download our requirements.txt file <../requirements.txt>`,
which contains a full list of all Python packages used by amplimap, and a known
working version.

Setup
~~~~~~~~~

To finish setting up amplimap you probably want to add the paths to the
reference genome files you will be using
(e.g. bwa index and reference genome fasta) to the :ref:`default-config`.

Requirements
~~~~~~~~~~~~~~~

- Linux environment (should also work on MacOS, Windows 10 Linux Subsystem)
- Python 3.5+ with setuptools, Cython and numpy

  - Further Python dependencies are listed in ``requirements.txt`` but can also be installed automatically by ``setup.py``

- Reference genome FASTA file, with indices

- Required software:

  - At least one read aligner: BWA (tested with v0.7.12), Bowtie2 (tested with v2.2.5), STAR (tested with v2.5.1b)
  - bedtools (tested with v2.27.1)
  - samtools (tested with v1.5)

- Additional software for germline variant calling (optional):

  - At least one variant caller: Platypus 0.8.1+, GATK 4+, Octopus
  - Annovar (tested with v2015-06-17)
  - bcftools (tested with v1.5)

- Additional software for low-frequency variant calling (optional):

  - Mutect2 (from GATK 4, tested with v4.0)

- Additional software for capture probe processing (optional):

  - Picard Tools 2+ (tested with v2.3.0)

If you do not have these required tools available yet, see :doc:`extended_installation`
for details on how to install them.
