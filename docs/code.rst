Code documentation
==================

Overview
------------------------
When ``amplimap`` is launched, the function amplimap.run.main() is called. This function
will perform various pre-processing steps (see below) and then launch Snakemake.
Snakemake, using the Snakefile provided in the amplimap package, will then determine
the jobs that need to be executed.

Testing amplimap
------------------------
To test amplimap, install the ``pytest`` package through pip and then
run this command in the amplimap source directory:

::

  pytest .


amplimap components
------------------------

Snakefile
^^^^^^^^^^^^^^^^

amplimap's ``Snakefile`` is used by the ``amplimap`` command-line executable
to determine the shell commands and library functions to execute to generate a given file.
Each "rule" in this file either executes a shell command (e.g. running ``bwa``),
executes a small set of Python commands or calls a function from the amplimap Python package.

Running ``snakemake`` directly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In theory, this file can also be used directly with the ``snakemake`` executable.
However, this will omit various steps performed by the ``amplimap`` executable,
including checking and loading the configuration files.
Thus, a configuration file needs to be provided with the ``--configfile`` parameter that contains all settings
included in ``config_default.yaml``.

In addition, you may want to add the following setting:

::

  general:
    amplimap_parent_dir: /path/to/parent/of/amplimap/package/directory

This directory will be inserted into the Python path when running Snakemake to make
sure ``import amplimap.xxx`` imports the correct files. If it is not provided
Snakemake will load the version of amplimap installed in the default Python path.


amplimap.run
^^^^^^^^^^^^^^^^

.. automodule:: amplimap.run
   :members:

amplimap.reader
^^^^^^^^^^^^^^^^

.. automodule:: amplimap.reader
   :members:

amplimap.parse_reads
^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: amplimap.parse_reads
   :members:

amplimap.naive_mapper
^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: amplimap.naive_mapper
   :members:

amplimap.pileup
^^^^^^^^^^^^^^^^

.. automodule:: amplimap.pileup
   :members:

amplimap.stats_alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: amplimap.stats_alignment
   :members:

amplimap.variants
^^^^^^^^^^^^^^^^^^

.. automodule:: amplimap.variants
   :members:

amplimap.coverage
^^^^^^^^^^^^^^^^^

.. automodule:: amplimap.coverage
   :members:

amplimap.merge_folders
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: amplimap.merge_folders
   :members:



Others
^^^^^^^^^^^^^^^^

.. automodule:: amplimap.common
   :members:

.. automodule:: amplimap.simulate
   :members:
