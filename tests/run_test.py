"""
This does some basic testing of the whole pipeline. Note this needs to be run through pytest and doesn't use unittest.
"""

#make sure we can import from package directory
import sys, os, re, gzip
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
test_config_path = os.path.join(packagedir, "sample_data", "config_default.yaml")
os.environ['AMPLIMAP_CONFIG'] = test_config_path

def init_wd(path, reads_in_path, umi_one = 0, umi_two = 0):
    assert os.path.isdir(path)

    #remove previous output
    shutil.rmtree(os.path.join(path, 'analysis'), ignore_errors=True)

    #remove previous reads_in and then prepare a new one
    shutil.rmtree(os.path.join(path, 'reads_in'), ignore_errors=True)
    os.mkdir(os.path.join(path, 'reads_in'))

    #turn .fastq files into .fastq.gz
    for file in os.listdir(reads_in_path):        
        if umi_one > 0 and '_R1' in file:
            umi_len = umi_one
        elif umi_two > 0 and '_R2' in file:
            umi_len = umi_two
        else:
            umi_len = 0

        next_umi = None
        with open(os.path.join(reads_in_path, file), 'rt') as fin, gzip.open(os.path.join(path, 'reads_in', '{}.gz'.format(file)), 'wt') as fout:
            for ix, line in enumerate(fin):
                assert line.startswith('@') or ix % 4 != 0
                assert line.startswith('+') or ix % 4 != 2

                if umi_len > 0:
                    if ix % 4 == 0: #name
                        match = re.search(r'_UMI-([^_]+)', line)
                        assert match, 'UMI missing: %s' % line
                        next_umi = match.group(1)
                    elif ix % 4 == 1: #seq
                        line = '{}{}'.format(next_umi[0:umi_len], line)
                    elif ix % 4 == 3: #qual
                        line = '{}{}'.format('A' * umi_len, line)

                fout.write(line)

def test_version(capsys):
    amplimap.run.main(['--version'])
    captured = capsys.readouterr()
    assert captured.out.strip() == '{} {}'.format(amplimap.run.__title__, amplimap.run.__version__)

def test_config(capsys):
    amplimap.run.main(['--print-config'])
    captured = capsys.readouterr()
    assert 'Reading additional configuration file: {}'.format(os.path.join(packagedir, "sample_data", "config_default.yaml")) in captured.err

def test_config_env(capsys):
    extra_config_path = os.path.join(packagedir, "sample_data", "extra_config.yaml")
    os.environ['AMPLIMAP_CONFIG'] = extra_config_path
    amplimap.run.main(['--print-config'])
    captured = capsys.readouterr()
    os.environ['AMPLIMAP_CONFIG'] = test_config_path #reset, so we don't affect later tests

    # #for debugging:
    # with capsys.disabled():
    #     sys.stdout.write(captured.err)
    #     sys.stdout.write(captured.out)

    assert 'Reading additional configuration file: {}'.format(extra_config_path) in captured.err
    assert 'aligner: star' in captured.out

def test_config_env_invalid(capsys):
    extra_config_path = os.path.join(packagedir, "sample_data", "extra_config_invalid.yaml")
    os.environ['AMPLIMAP_CONFIG'] = extra_config_path
    amplimap.run.main(['--print-config'])
    captured = capsys.readouterr()
    os.environ['AMPLIMAP_CONFIG'] = test_config_path #reset, so we don't affect later tests

    assert 'Reading additional configuration file: {}'.format(extra_config_path) in captured.err
    assert 'Your configuration file(s) contain unknown or invalid settings:' in captured.err

def check_run(capsys, wd_path, rules = 'pileups'):
    #dry-run
    amplimap.run.main(['--working-directory={}'.format(wd_path), rules])
    captured = capsys.readouterr()
    assert '{} {} dry run successful.'.format(amplimap.run.__title__, amplimap.run.__version__) in captured.err.strip()

    #full run
    amplimap.run.main(['--working-directory={}'.format(wd_path), rules, '--run'])
    captured = capsys.readouterr()
    assert '{} {} finished!'.format(amplimap.run.__title__, amplimap.run.__version__) in captured.err.strip()    

def check_default_stats(wd_path, is_trimmed = True):
    samples = pd.read_csv(os.path.join(wd_path, 'analysis', 'reads_parsed', 'stats_samples.csv'))
    assert len(samples) == 1

    assert samples.loc[0, 'sample'] == 'S1'
    assert samples.loc[0, 'files'] == 2
    assert samples.loc[0, 'pairs_total'] == 7
    assert samples.loc[0, 'pairs_good_arms'] == 6
    assert samples.loc[0, 'pairs_r1_too_short'] == (1 if is_trimmed else 0)

def check_default_pileups(wd_path, expected_coverage = 5, include_too_short = False):
    pileups = pd.read_csv(os.path.join(wd_path, 'analysis', 'pileups', 'pileups_long.csv'))

    #we covered 11bp
    assert len(pileups) == 11

    #we should have 5 reads, except for the raw alignments where we include the pair with short r1/r2
    pileups['expected_coverage'] = expected_coverage
    if include_too_short:
        pileups.loc[(pileups.pos <= 32) | (pileups.pos >= 39), 'expected_coverage'] += 1

    #everything but 35/37 should be ref
    assert pileups.loc[~pileups.pos.isin([35,37]), 'alts'].isnull().all()
    assert (pileups.loc[~pileups.pos.isin([35,37]), 'ref_hq_count'] == pileups.loc[~pileups.pos.isin([35,37]), 'expected_coverage']).all()
    assert (pileups.loc[~pileups.pos.isin([35,37]), 'nonref_hq_count'] == 0).all()

    #only these should be nonref (36 is low-quality in one read)
    assert pileups.loc[pileups.pos.isin([35,37]), 'alts'].notnull().all()
    #one read from L001, one from L002
    assert (pileups.loc[pileups.pos == 35, 'nonref_hq_count'] == 2).all()
    assert (pileups.loc[pileups.pos == 35, 'ref_hq_count'] == pileups.loc[pileups.pos == 35, 'expected_coverage'] - 2).all()
    assert (set(pileups.loc[pileups.pos == 35, 'alts'].iloc[0].split(';')) == set(['A', 'G'])) #explicitly use iloc and no .all() here
    #just one in L001
    assert (pileups.loc[pileups.pos == 37, 'nonref_hq_count'] == 1).all()
    assert (pileups.loc[pileups.pos == 37, 'ref_hq_count'] == pileups.loc[pileups.pos == 37, 'expected_coverage'] - 1).all()
    assert (set(pileups.loc[pileups.pos == 37, 'alts'].iloc[0].split(';')) == set(['G'])) #explicitly use iloc and no .all() here

def test_naive_pileups_simulation(capsys):
    wd_path = os.path.join(packagedir, "sample_data", "wd_naive")
    init_wd(wd_path, os.path.join(packagedir, "sample_data", "sample_reads_in"))
    shutil.rmtree(os.path.join(wd_path, 'test__GGCAATATGT_GGCAATCTGT_100'), ignore_errors=True) #make sure this doesn't exist

    #run simulation, replacing A@30 to C
    amplimap.run.main(['--working-directory={}'.format(wd_path), 'test__GGCAATATGT_GGCAATCTGT_100/test_pileups.done', '--run'])
    captured = capsys.readouterr()
    assert '{} {} finished!'.format(amplimap.run.__title__, amplimap.run.__version__) in captured.err.strip()

    pileups = pd.read_csv(os.path.join(wd_path, 'test__GGCAATATGT_GGCAATCTGT_100', 'pileups', 'pileups_long.csv'))    
    assert len(pileups) == 11

    #we should have an A>C SNP at pos 30, in addition to the others
    assert pileups.loc[~pileups.pos.isin([30, 35, 37]), 'alts'].isnull().all()
    assert (pileups.loc[pileups.pos == 30, 'nonref_hq_count'] == 5).all()
    assert (pileups.loc[pileups.pos == 30, 'ref_hq_count'] == 0).all()
    assert (set(pileups.loc[pileups.pos == 30, 'alts'].iloc[0].split(';')) == set(['C'])) #explicitly use iloc and no .all() here


def test_naive_pileups(capsys):
    wd_path = os.path.join(packagedir, "sample_data", "wd_naive")
    init_wd(wd_path, os.path.join(packagedir, "sample_data", "sample_reads_in"))

    check_run(capsys, wd_path)
    check_default_stats(wd_path)
    check_default_pileups(wd_path)

def test_naive_pileups_notrim(capsys):
    wd_path = os.path.join(packagedir, "sample_data", "wd_naive_notrim")
    init_wd(wd_path, os.path.join(packagedir, "sample_data", "sample_reads_in"))

    check_run(capsys, wd_path)
    check_default_stats(wd_path, is_trimmed = False)
    #check_default_pileups(wd_path, expected_coverage = 6)

def test_bwa_pileups(capsys):
    wd_path = os.path.join(packagedir, "sample_data", "wd_bwa")
    init_wd(wd_path, os.path.join(packagedir, "sample_data", "sample_reads_in"))

    check_run(capsys, wd_path)
    check_default_stats(wd_path)
    check_default_pileups(wd_path)

def test_bwa_pileups_notrim(capsys):
    wd_path = os.path.join(packagedir, "sample_data", "wd_bwa_notrim")
    init_wd(wd_path, os.path.join(packagedir, "sample_data", "sample_reads_in"))

    check_run(capsys, wd_path)
    check_default_stats(wd_path, is_trimmed = False)
    #check_default_pileups(wd_path, expected_coverage = 6)

def test_raw_read_pileups(capsys):
    wd_path = os.path.join(packagedir, "sample_data", "wd_bwa_raw")
    init_wd(wd_path, os.path.join(packagedir, "sample_data", "sample_reads_in"))

    check_run(capsys, wd_path)

    #check custom stats
    samples = pd.read_csv(os.path.join(wd_path, 'analysis', 'reads_parsed', 'stats_samples.csv'))
    assert len(samples) == 1

    assert samples.loc[0, 'sample'] == 'S1'
    assert samples.loc[0, 'files'] == 2
    assert samples.loc[0, 'pairs_total'] == 7
    assert samples.loc[0, 'pairs_good_arms'] == 7
    assert samples.loc[0, 'pairs_r1_too_short'] == 0

    check_default_pileups(wd_path, include_too_short = True)

def test_umi_pileups(capsys):
    wd_path = os.path.join(packagedir, "sample_data", "wd_umis")
    init_wd(wd_path, os.path.join(packagedir, "sample_data", "sample_reads_in"), umi_one = 3, umi_two = 4)

    check_run(capsys, wd_path)
    check_default_stats(wd_path)
    check_default_pileups(wd_path, expected_coverage=4) #one less, because two read pairs have same umi!

    #check umi-specific stats
    stats_reads = pd.read_csv(os.path.join(wd_path, 'analysis', 'reads_parsed', 'stats_reads.csv'))
    assert len(stats_reads) == 1
    assert stats_reads.loc[0, 'sample'] == 'S1'
    assert stats_reads.loc[0, 'probe'] == 'Probe1'
    assert stats_reads.loc[0, 'read_pairs'] == 6
    assert stats_reads.loc[0, 'umis_total'] == 5
    assert stats_reads.loc[0, 'umis_coverage_max'] == 2

def test_umi_dedup(capsys):
    wd_path = os.path.join(packagedir, "sample_data", "wd_umis")
    init_wd(wd_path, os.path.join(packagedir, "sample_data", "sample_reads_in"), umi_one = 3, umi_two = 4)

    check_run(capsys, wd_path, rules = 'dedup_bams')
    
    import pysam
    #before dedup we had five aligned read pairs
    n_in = pysam.AlignmentFile(os.path.join(wd_path, 'analysis', 'bams', 'S1.bam')).count(until_eof=True)
    assert n_in == 5 * 2

    #after dedup we have four (two read pairs have same UMI)
    n_dedup = pysam.AlignmentFile(os.path.join(wd_path, 'analysis', 'bams_umi_dedup', 'S1.bam')).count(until_eof=True)
    assert n_dedup == 4 * 2

