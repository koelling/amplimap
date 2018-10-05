#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
if sys.version_info[0] < 3:
    sys.exit("Sorry, amplimap requires at least Python 3")

from setuptools import setup, find_packages, Extension

#load version number (can't just import the package since we might miss requirements)
filename = 'amplimap/version.py'
exec(compile(open(filename, "rb").read(), filename, 'exec')) #python2/3 compatible replacement for execfile()

#load long description from readme
import codecs
with codecs.open('README.rst', 'r', 'utf-8') as f:
    long_description = f.read()

setup(
    name = __title__,
    version = __version__,

    packages = find_packages(),
    #these files will be added to the package directory
    package_data={ '': ['parse_reads_cy.pyx', 'Snakefile', 'config_default.yaml'] },

    # metadata for upload to PyPI
    author = "Nils Koelling",
    author_email = "git@nk.gl",
    description = "amplicon/smMIP mapping and analysis pipeline",
    long_description = long_description,
    license = "Apache License, Version 2.0",
    keywords = "amplimap amplicon smmip mapping analysis pipeline",
    url = "https://github.com/koelling/amplimap/",
    download_url="https://github.com/koelling/amplimap/archive/v%s.tar.gz" % __version__,
    platforms=["any"],
    entry_points={
        'console_scripts': [
            'amplimap = amplimap.run:main',
            'amplimap_merge = amplimap.merge_folders:main',
            'amplimap_setup = amplimap.run_setup:main',
            #'amplimap_pileup = amplimap.pileup:main',
        ]
    },

    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3',
    ],

    python_requires='>=3',

    install_requires=[
        'snakemake>=3.11.2',
        'pyyaml>=3.12',
        'numpy>=1.13.1',
        'biopython>=1.69',
        'pandas>=0.20.3',
        'interlap>=0.2.5',
        'pysam>=0.11.1,<0.14', #pileups seem to have some issues in 0.14, so we force something in-between here
        'pyfaidx>=0.4.8.4',
        'distance>=0.1.3',
        'umi_tools>=0.5.0',
    ],

    setup_requires=[
        'setuptools>=18.0',
        'cython',
    ],

    ext_modules = [Extension("amplimap.parse_reads_cy", ["amplimap/parse_reads_cy.pyx"])]
)
