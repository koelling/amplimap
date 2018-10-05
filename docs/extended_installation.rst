=========================================
Extended installation guide for amplimap
=========================================

If you need to install amplimap and its requirements (such as Python, read aligners, etc) from scratch
we recommend that you use `bioconda <https://bioconda.github.io/>`_.

1. Install bioconda
--------------------
Follow the official `installation guide <https://conda.io/docs/user-guide/install/index.html>`_ on the
conda website to install Miniconda.

Then, set up the bioconda channels as described in the `bioconda documentation <https://bioconda.github.io/>`_.

2. Setup amplimap enviroment
--------------------------------
Download :download:`amplimap's environment file <../environment.yml>`, which contains a list of all the software used by amplimap.

Use the ``conda`` tool to install this software into a conda environment:

::
    
    conda env create -f environment.yml

.. conda create --name amplimap 'python>=3.4' pip setuptools numpy cython bwa bowtie2 star bedtools samtools bcftools gatk4 picard
.. source activate amplimap
.. #conda env export > environment.yml

If you want to run germline variant calling and annotation, you also need to `download and install
Annovar <http://annovar.openbioinformatics.org/en/latest/user-guide/download/>`_ manually. Make sure you also download
the relevant indices for the reference genome you want to use, 


3. Activate amplimap
------------------------------------------------
Load amplimap by running this command:

::

    source activate amplimap

You only need to run this command once per session, eg. when you open a new terminal window.

Run ``amplimap --version`` to confirm that amplimap has been installed and activated.


4. Set up your reference genome and indices
-------------------------------------------
Download the reference genome FASTA file that you want to use, for example from the `Ensembl
FTP <https://www.ensembl.org/info/data/ftp/index.html>`_. When in doubt, we recommend using the
primary_assembly version, for example ``Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz``.

Once you have downloaded this file, you need to prepare it for use in amplimap:

::

    #decompress the file (if it ends in .gz)
    gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

    #index FASTA file
    samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

    #build bwa index, if you want to use bwa:
    bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa

    #build bowtie2 index, if you want to use bowtie2:
    bowtie2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa


5. Update amplimap configuration
------------------------------------------

Finally, you need to add the paths of the reference genome files to your ``config_default.yaml``.
To find out where this file is located, run:

::

    amplimap --basedir

Open the file ``config_default.yaml`` at this location and look for the settings under ``paths:``
corresponding to the indices you created.
Replace these with the full paths to your files. If you haven't generated the corresponding
file, leave the setting empty. For example, if you generated indices for bwa and bowtie2 (but not STAR or Annovar)
and always used the same FASTA filename as the prefix:

::

    paths:
      hg38:
        bwa: "/home/user/amplimap/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        bowtie2: "/home/user/amplimap/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        fasta: "/home/user/amplimap/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

If you used a different reference genome, change the "hg38" to the appropriate abbreviation (eg. "mm10")
and also update the line ``genome_name: "hg38"`` below.

If you are using Annovar, make sure you also provide the path to its indices directory under ``paths:``
and adjust the protocols/operations under ``annotate: annovar: protocols:`` to match the indices you
have downloaded.

Save the file and confirm that the settings are being read correctly by looking at the output of ``amplimap --print-config``.

6. Run amplimap!
-------------------
Now you are ready to run amplimap! When you start a new session, activate the conda environment
as described above and then run the amplimap commands as usual:

::

    source activate amplimap

    amplimap --version
