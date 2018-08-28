"""
This does some basic testing of the whole pipeline. Note this needs to be run through pytest and doesn't use unittest.
"""

#make sure we can import from package directory
import sys, os
packagedir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, packagedir) 

import shutil
import subprocess
import pandas as pd

import amplimap.run

#we need to build the cython file here, since snakemake will call an external python that wouldn't
#inherit pyximport
os.system("cythonize -i {}".format(os.path.join(packagedir, "amplimap", "parse_reads_cy.pyx")))

#set config
os.environ['AMPLIMAP_CONFIG'] = os.path.join(packagedir, "sample_data", "config_default.yaml")

def init_wd(path, reads_in_path):
    assert os.path.isdir(path)

    #remove previous output
    shutil.rmtree(os.path.join(path, 'analysis'), ignore_errors=True)

    #remove previous reads_in and then prepare a new one
    shutil.rmtree(os.path.join(path, 'reads_in'), ignore_errors=True)
    os.mkdir(os.path.join(path, 'reads_in'))

    #turn .fastq files into .fastq.gz
    for file in os.listdir(reads_in_path):
        os.system("gzip -c {} > {}".format(os.path.join(reads_in_path, file), os.path.join(path, 'reads_in', '{}.gz'.format(file))))

def test_version(capsys):
    amplimap.run.main(['--version'])
    captured = capsys.readouterr()
    assert captured.out.strip() == '{} {}'.format(amplimap.run.__title__, amplimap.run.__version__)

def test_config(capsys):
    amplimap.run.main(['--print-config'])
    captured = capsys.readouterr()
    assert 'Reading additional configuration from: {}'.format(os.path.join(packagedir, "sample_data", "config_default.yaml")) in captured.err
    # with capsys.disabled():
    #     sys.stdout.write(captured.err)
    #     sys.stdout.write(captured.out)

def test_naive_pileups(capsys):
    wd_path = os.path.join(packagedir, "sample_data", "wd_naive")
    init_wd(wd_path, os.path.join(packagedir, "sample_data", "sample_reads_in"))

    #dry-run
    amplimap.run.main(['--working-directory={}'.format(wd_path), 'pileups'])
    captured = capsys.readouterr()
    assert '{} {} dry run successful.'.format(amplimap.run.__title__, amplimap.run.__version__) in captured.err.strip()

    #full run
    amplimap.run.main(['--working-directory={}'.format(wd_path), 'pileups', '--run'])
    captured = capsys.readouterr()
    assert '{} {} finished!'.format(amplimap.run.__title__, amplimap.run.__version__) in captured.err.strip()

    #check sample stats
    samples = pd.read_csv(os.path.join(wd_path, 'analysis', 'reads_parsed', 'stats_samples.csv'))
    assert len(samples) == 1

    assert samples.loc[0, 'sample'] == 'S1'
    assert samples.loc[0, 'files'] == 1
    assert samples.loc[0, 'pairs_total'] == 6
    assert samples.loc[0, 'pairs_good_arms'] == 5

    #check pileups
    pileups = pd.read_csv(os.path.join(wd_path, 'analysis', 'pileup', 'pileups_long.csv'))
    assert len(pileups) == 11

    assert pileups.loc[~pileups.pos.isin([35,37]), 'alts'].isnull().all()
    assert (pileups.loc[~pileups.pos.isin([35,37]), 'ref_hq_count'] == 4).all()
    assert (pileups.loc[~pileups.pos.isin([35,37]), 'nonref_hq_count'] == 0).all()

    #only these should be nonref (36 is low-quality in one read)
    assert pileups.loc[pileups.pos.isin([35,37]), 'alts'].notnull().all()
    assert (pileups.loc[pileups.pos.isin([35,37]), 'ref_hq_count'] == 3).all()
    assert (pileups.loc[pileups.pos.isin([35,37]), 'nonref_hq_count'] == 1).all()

def test_bwa_pileups(capsys):
    wd_path = os.path.join(packagedir, "sample_data", "wd_bwa")
    init_wd(wd_path, os.path.join(packagedir, "sample_data", "sample_reads_in"))

    #dry-run
    amplimap.run.main(['--working-directory={}'.format(wd_path), 'pileups'])
    captured = capsys.readouterr()
    assert '{} {} dry run successful.'.format(amplimap.run.__title__, amplimap.run.__version__) in captured.err.strip()

    #full run
    amplimap.run.main(['--working-directory={}'.format(wd_path), 'pileups', '--run'])
    captured = capsys.readouterr()
    assert '{} {} finished!'.format(amplimap.run.__title__, amplimap.run.__version__) in captured.err.strip()

    #check sample stats
    samples = pd.read_csv(os.path.join(wd_path, 'analysis', 'reads_parsed', 'stats_samples.csv'))
    assert len(samples) == 1

    assert samples.loc[0, 'sample'] == 'S1'
    assert samples.loc[0, 'files'] == 1
    assert samples.loc[0, 'pairs_total'] == 6
    assert samples.loc[0, 'pairs_good_arms'] == 5

    #check pileups
    pileups = pd.read_csv(os.path.join(wd_path, 'analysis', 'pileup', 'pileups_long.csv'))
    assert len(pileups) == 11

    assert pileups.loc[~pileups.pos.isin([35,37]), 'alts'].isnull().all()
    assert (pileups.loc[~pileups.pos.isin([35,37]), 'ref_hq_count'] == 4).all()
    assert (pileups.loc[~pileups.pos.isin([35,37]), 'nonref_hq_count'] == 0).all()

    #only these should be nonref (36 is low-quality in one read)
    assert pileups.loc[pileups.pos.isin([35,37]), 'alts'].notnull().all()
    assert (pileups.loc[pileups.pos.isin([35,37]), 'ref_hq_count'] == 3).all()
    assert (pileups.loc[pileups.pos.isin([35,37]), 'nonref_hq_count'] == 1).all()

def test_raw_read_pileups(capsys):
    wd_path = os.path.join(packagedir, "sample_data", "wd_bwa_raw")
    init_wd(wd_path, os.path.join(packagedir, "sample_data", "sample_reads_in"))

    #dry-run
    amplimap.run.main(['--working-directory={}'.format(wd_path), 'pileups'])
    captured = capsys.readouterr()
    assert '{} {} dry run successful.'.format(amplimap.run.__title__, amplimap.run.__version__) in captured.err.strip()

    #full run
    amplimap.run.main(['--working-directory={}'.format(wd_path), 'pileups', '--run'])
    captured = capsys.readouterr()
    assert '{} {} finished!'.format(amplimap.run.__title__, amplimap.run.__version__) in captured.err.strip()

    #check sample stats
    samples = pd.read_csv(os.path.join(wd_path, 'analysis', 'reads_parsed', 'stats_samples.csv'))
    assert len(samples) == 1

    assert samples.loc[0, 'sample'] == 'S1'
    assert samples.loc[0, 'files'] == 1
    assert samples.loc[0, 'pairs_total'] == 6
    assert samples.loc[0, 'pairs_good_arms'] == 6

    #check pileups
    pileups = pd.read_csv(os.path.join(wd_path, 'analysis', 'pileup', 'pileups_long.csv'))
    assert len(pileups) == 11

    assert pileups.loc[~pileups.pos.isin([35,37]), 'alts'].isnull().all()
    assert (pileups.loc[~pileups.pos.isin([35,37]), 'ref_hq_count'] == 4).all()
    assert (pileups.loc[~pileups.pos.isin([35,37]), 'nonref_hq_count'] == 0).all()

    #only these should be nonref (36 is low-quality in one read)
    assert pileups.loc[pileups.pos.isin([35,37]), 'alts'].notnull().all()
    assert (pileups.loc[pileups.pos.isin([35,37]), 'ref_hq_count'] == 3).all()
    assert (pileups.loc[pileups.pos.isin([35,37]), 'nonref_hq_count'] == 1).all()