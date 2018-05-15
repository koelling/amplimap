"""
This does some basic testing of the whole pipeline. Note this needs to be run through pytest and doesn't use unittest.
"""

#make sure we can import from package directory
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) 

import shutil
import subprocess
import pandas as pd

import amplimap.run

#use this to load parse_reads_cy properly
import pyximport; pyximport.install()

#set config
os.environ['AMPLIMAP_CONFIG'] = "sample_data/config_default.yaml"

def init_wd(path):
    assert os.path.isdir(path)

    #remove previous output
    shutil.rmtree(os.path.join(path, 'analysis'), ignore_errors=True)

    #turn .fastq files into .fastq.gz
    os.system("gzip -c {}/reads_in/S1_L001_R1_001.fastq > {}/reads_in/S1_L001_R1_001.fastq.gz".format(path, path))
    os.system("gzip -c {}/reads_in/S1_L001_R2_001.fastq > {}/reads_in/S1_L001_R2_001.fastq.gz".format(path, path))

def test_version(capsys):
    amplimap.run.main(['--version'])
    captured = capsys.readouterr()
    assert captured.out.strip() == '{} {}'.format(amplimap.run.__title__, amplimap.run.__version__)

def test_config(capsys):
    amplimap.run.main(['--print-config'])
    captured = capsys.readouterr()
    assert 'Reading additional configuration from: sample_data/config_default.yaml' in captured.err
    # with capsys.disabled():
    #     sys.stdout.write(captured.err)
    #     sys.stdout.write(captured.out)

def test_normal_pileups(capsys):
    init_wd("sample_data/wd1/")

    #dry-run
    amplimap.run.main(['--working-directory={}'.format("sample_data/wd1/"), 'pileups'])
    captured = capsys.readouterr()
    assert '{} {} dry run successful.'.format(amplimap.run.__title__, amplimap.run.__version__) in captured.err.strip()

    #full run
    amplimap.run.main(['--working-directory={}'.format("sample_data/wd1/"), 'pileups', '--run'])
    captured = capsys.readouterr()
    assert '{} {} finished!'.format(amplimap.run.__title__, amplimap.run.__version__) in captured.err.strip()

    #check sample stats
    samples = pd.read_csv('sample_data/wd1/analysis/reads_parsed/stats_samples.csv')
    assert len(samples) == 1

    assert samples.loc[0, 'sample'] == 'S1'
    assert samples.loc[0, 'files'] == 1
    assert samples.loc[0, 'pairs_total'] == 6
    assert samples.loc[0, 'pairs_good_arms'] == 5

    #check pileups
    pileups = pd.read_csv('sample_data/wd1/analysis/pileup/pileups_long.csv')
    assert len(pileups) == 11

    assert pileups.loc[~pileups.pos.isin([35,37]), 'alts'].isnull().all()
    assert (pileups.loc[~pileups.pos.isin([35,37]), 'ref_hq_count'] == 4).all()
    assert (pileups.loc[~pileups.pos.isin([35,37]), 'nonref_hq_count'] == 0).all()

    #only these should be nonref (36 is low-quality in one read)
    assert pileups.loc[pileups.pos.isin([35,37]), 'alts'].notnull().all()
    assert (pileups.loc[pileups.pos.isin([35,37]), 'ref_hq_count'] == 3).all()
    assert (pileups.loc[pileups.pos.isin([35,37]), 'nonref_hq_count'] == 1).all()