=========================================
Installation guide
=========================================

We recommend that you use `Miniconda <https://conda.io/en/latest/miniconda.html>`_
with `bioconda <https://bioconda.github.io/>`_ to install amplimap and its requirements
(such as Python 3.6, read aligners, etc). If you have `Docker <https://www.docker.com/>`_ you can
also use our Dockerfile instead: :ref:`installation-docker`.

If your machine already has all of the required
software installed you can also install amplimap through pip.
For more details, see :ref:`installation-pip`.

.. _installation-miniconda:

Installing amplimap through Conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Install Miniconda 3
-----------------------
Download and install `Miniconda with Python 3 <https://conda.io/en/latest/miniconda.html>`_
from the conda website.


2. Install amplimap environment
--------------------------------
Download amplimap's :download:`environment.yml file <../environment.yml>`
and use it to install amplimap and its requirements
into a new conda environment:

::

    conda env create --file environment.yml

.. conda create --name amplimap 'python>=3.4' pip setuptools numpy cython bwa bowtie2 star bedtools samtools bcftools gatk4 picard
.. conda activate amplimap
.. #conda env export > environment.yml

If you want to run germline variant calling and annotation you also need to `download and install
Annovar <http://annovar.openbioinformatics.org/en/latest/user-guide/download/>`_ manually. Make sure you also download
the relevant indices for the reference genome you want to use.


3. Activate amplimap environment
------------------------------------------------
Load the amplimap environnent by running this command:

::

    conda activate amplimap

You only need to run this command once per session,
e.g. when you open a new terminal window.

Your command line prompt should now start with ``(amplimap)``.
Run ``amplimap --version`` to confirm that the correct version of
amplimap has been installed and activated.


.. _installation-setup:

4. Set up your reference genome and indices
-------------------------------------------
Download the DNA (FASTA) file for the reference genome that you want to use, for example from the `Ensembl
FTP <https://www.ensembl.org/info/data/ftp/index.html>`_
or `iGenomes <https://support.illumina.com/sequencing/sequencing_software/igenome.html>`_.
When in doubt we recommend using the
primary_assembly file from Ensembl, for example ``Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz``.

Once you have downloaded this file you need to prepare it for use in amplimap:

::

    # decompress the file (if it ends in .gz)
    gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    # index FASTA file
    samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
    # create dictionary for GATK
    picard CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna.primary_assembly.fa
    # build bwa index, if you want to use bwa:
    bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
    # build bowtie2 index, if you want to use bowtie2:
    bowtie2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa


5. Update amplimap configuration
------------------------------------------

Finally, we recommend that you add the paths of the reference genome files to your ``config_default.yaml``.
This way, you don't need to specify these paths in every single directory-specific ``config.yaml``.
To find out where this file is located run:

::

    amplimap --basedir

Open the file ``config_default.yaml`` at this location and look for the settings under ``paths:``
corresponding to the indices you created.

Replace these with the full paths to your files. If you haven't generated one of the
files leave the corresponding setting empty.
For example, if you generated indices for bwa and bowtie2 (but not STAR or Annovar)
and always used the same FASTA filename as the prefix:

::

    paths:
      hg38:
        bwa: "/home/user/amplimap/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        bowtie2: "/home/user/amplimap/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        fasta: "/home/user/amplimap/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

If you are working with a different reference genome change ``hg38:`` to the appropriate abbreviation (e.g. ``mm10:``)
and also update the line ``genome_name: "hg38"`` below.

If you are using Annovar make sure you also provide the path to its indices directory under ``paths:``
and adjust the protocols/operations under ``annotate: annovar: protocols:`` to match the indices you
have downloaded.

Save the file and confirm that the settings are being read correctly by looking at the output of ``amplimap --print-config``.

6. Run amplimap!
-------------------
Now you are ready to run amplimap! Prepare a working directory
(see :ref:`usage`), change into it using ``cd`` and then run
``amplimap`` to get started.

If you get a message about the command not being found
please make sure you activated the conda environment as described above.

.. _installation-docker:

Installing amplimap through Docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We also have a `Docker image <https://hub.docker.com/r/koelling/amplimap>`_
available.
To use this, install `Docker <https://www.docker.com/>`_ and then
prefix your amplimap commands with ``docker run koelling/amplimap``,
forwarding directories from your host into the docker container using
Docker's ``-v`` parameter.

For example, here are some commands you could use to prepare
indices for an *E. coli* :download:`reference genome FASTA <../sample_data/ecoli.fasta>`
located under ``~/references/ecoli.fasta`` and then run amplimap
on some :download:`example data <../sample_data/example_wd.tar>`
located in ``~/data/example_wd``:

::

    # download the docker image (only need to run this once)
    docker pull koelling/amplimap

    # check version
    docker run koelling/amplimap amplimap --version

    # build indices for ~/references/ecoli.fasta
    docker run -v ~/references:/references koelling/amplimap samtools faidx /references/ecoli.fasta
    docker run -v ~/references:/references koelling/amplimap picard CreateSequenceDictionary R=/references/ecoli.fasta
    docker run -v ~/references:/references koelling/amplimap bwa index /references/ecoli.fasta

    # run amplimap with working directory ~/data/example_wd
    docker run -v ~/references:/references -v ~/data:/data koelling/amplimap amplimap --working-directory=/data/example_wd coverages pileups variants

Note that in this example you would have to provide the paths to your reference genome
in the ``~/data/example_wd/config.yaml`` file:

::

    paths:
      ecoli:
        bwa: "/references/ucsc.ecoli.fasta"
        fasta: "/references/ucsc.ecoli.fasta"
    general:
      genome_name: "ecoli"

You can avoid having to specify these paths every time by running a shell inside the Docker container
and adding your reference genome to your ``config_default.yaml`` as described here: :ref:`installation-setup`.

::

    docker run -t -i koelling/amplimap /bin/bash

To annotate variant calls you would also have to install Annovar inside the Docker container
and add the path to the Annovar indices to your config.


.. _installation-pip:

Installing amplimap through pip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you already have all of the required external software available
(see :ref:`installation-requirements`)
you can install amplimap directly through pip.
Please note that this **requires Python 3.5 or 3.6** and does not currently
work with Python 3.7 due to problems with the pysam package. It also
does not work with any Python version lower than 3.5.

If you do not have the dependencies and the right version of Python available
please see :ref:`installation-miniconda`.

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

To finish setting up amplimap you probably want to add the paths to the
reference genome files you will be using
(e.g. bwa index and reference genome fasta) to the :ref:`default-config`.
See :ref:`installation-setup` for more details.

.. _installation-requirements:

Requirements
~~~~~~~~~~~~~~~
Please note that, other than the Linux environment and the reference genome files,
all requirements **will be installed automatically** when you install amplimap
through conda.

- Linux environment (should also work on MacOS, Windows 10 Linux Subsystem)
- Python 3.5 or 3.6 with setuptools, Cython and numpy

  - Further Python dependencies are listed in ``requirements.txt``
    but can also be installed automatically by ``setup.py``.

- Required software:

  - At least one read aligner: BWA (tested with v0.7.12),
    Bowtie2 (tested with v2.2.5), STAR (tested with v2.5.1b)
  - bedtools (tested with v2.27.1)
  - samtools (tested with v1.5)

- Additional software for germline variant calling (optional):

  - At least one variant caller: Platypus 0.8.1+, GATK 4+
  - Annovar (tested with v2015-06-17)
  - bcftools (tested with v1.5)

- Additional software for low-frequency variant calling (optional):

  - Mutect2 (from GATK 4, tested with v4.0)

- Additional software for capture probe processing (optional):

  - Picard Tools 2+ (tested with v2.3.0)

- Reference genome FASTA file, with indices
