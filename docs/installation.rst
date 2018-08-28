.. _installation:

Installation
------------
amplimap is available through pip, so this command should install
the pipeline as well as all required python packages.

::

  #you may need to use `pip` instead of `pip3`
  pip3 install amplimap

If this does not work, you can try to install it manually:

::

  #install required python3 packages
  #you may need to use `pip` instead of `pip3`
  pip3 install setuptools Cython numpy

  #download and install amplimap
  #you may need to use `python` instead of `python3`
  git clone --depth=1 https://github.com/koelling/amplimap.git
  cd amplimap
  python3 setup.py install

Setup
~~~~~~~~~

To finish setting up amplimap you probably want to add the paths to the
reference genome files you will be using
(eg. bwa index and reference genome fasta) to the :ref:`default-config`.

Requirements
~~~~~~~~~~~~~~~

- Linux environment (should also work on MacOS, Windows 10 Linux Subsystem)
- Python 3.5+ with setuptools, Cython and numpy

  - Further Python dependencies will be installed automatically

- Reference genome FASTA with indices

- Required software:

  - At least one read aligner: bwa (tested with v0.7.12), Bowtie2 (tested with v2.2.5), STAR (tested with v2.5.1b)
  - bedtools (tested with v2.27.1)
  - samtools (tested with v1.5)

- Additional software for germline variant analysis (optional):

  - At least one variant caller: Platypus 0.8.1+, GATK 4+
  - Annovar (tested with v2015-06-17)
  - bcftools (tested with v1.5)

- Additional software for somatic variant analysis (optional):

  - Mutect2 (from GATK 4, tested with v4.0)

- Additional software for capture probe processing (optional):

  - Picard Tools 2+ (tested with v2.3.0)

If you do not have these required tools available yet, see :doc:`extended_installation`
for details on how to install them.

