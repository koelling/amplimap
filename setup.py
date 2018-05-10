#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
if sys.version_info[0] < 3:
    sys.exit("Sorry, amplimap requires at least Python 3")

from setuptools import setup, find_packages
#"Note that when using setuptools, you should import it before Cython"
try:
    from Cython.Build import cythonize #this will fail if cython is not installed
except ImportError:
    raise ImportError('Could not load Cython! Please make sure the cython package is installed.')

#load version number (can't just import the package since we might miss requirements)
filename = 'amplimap/version.py'
exec(compile(open(filename, "rb").read(), filename, 'exec')) #python2/3 compatible replacement for execfile()

#load long description from readme
with open('README.rst') as f:
    long_description = f.read()

setup(
    name = __title__,
    version = __version__,

    packages = find_packages(),
    #these files will be added to the package directory
    package_data={ '': ['parse_reads_cy.pyx'] },
    #these files should be one level above the package directory (amplimap)
    data_files = [ ('', ['Snakefile', 'config_default.yaml']) ],

    # metadata for upload to PyPI
    author = "Nils Koelling",
    author_email = "git@nk.gl",
    description = "amplicon/smMIP mapping and analysis pipeline",
    long_description = long_description,
    #license = "MIT",
    keywords = "amplimap amplicon smmip mapping analysis pipeline",
    url = "https://github.com/koelling/amplimap/",
    download_url="https://github.com/koelling/amplimap/archive/v%s.tar.gz" % __version__,
    platforms=["any"],
    entry_points={
        'console_scripts': [
            'amplimap = amplimap.run:main',
            'amplimap_merge = amplimap.merge_folders:main',
            'amplimap_pileup = amplimap.pileup:main',
        ]
    },

    classifiers=[
        'Development Status :: 4 - Beta',
        #'License :: OSI Approved :: MIT License',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3',
    ],

    python_requires='>=3',

    install_requires=[
        'snakemake',
        'pyyaml',
        'six',
        'numpy',
        'biopython',
        'pandas',
        'argparse',
        'interlap',
        'pysam',
        'pyfaidx',
        'distance',
        'umi_tools',
    ],

    ext_modules = cythonize("amplimap/parse_reads_cy.pyx")
)
