#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module provides the ``amplimap.run.main()`` method called from the `amplimap` command-line executable,
as well as some helper functions for reading and checking the config files.
"""

import os
import sys

import snakemake
import argparse
import yaml

from .version import __title__, __version__
from .reader import AmplimapReaderException, read_new_probe_design, process_probe_design, read_and_convert_mipgen_probes, read_and_convert_heatseq_probes, read_targets, read_snps_txt

def check_config_keys(default_config, my_config, path = []):
    """Recursively check that config keys provided in my_config also exist in default_config (ignoring 'paths' and 'clusters')."""
    differences = []
    for key, value in my_config.items():
        if key in default_config:
            if isinstance(default_config[key], dict) == isinstance(value, dict):
                if isinstance(default_config[key], dict):
                    #both are dicts, keep checking
                    #but don't check paths/cluster because they can contain custom values
                    if not (len(path) == 0 and key in ['paths', 'clusters']):
                        differences += check_config_keys(default_config[key], value, path + [key])
                else:
                    #we're done here
                    pass
            else:
                differences.append(path + [key])
                #raise Exception('Config setting {} is invalid\n'.format(':'.join(path+[key])))
        else:
            differences.append(path + [key])
    return differences

def compare_config_dicts(my_config, used_config, path = []):
    """Recursively search for differences in values between two dicts."""
    differences = []
    for key, value in my_config.items():
        if key in used_config:
            if isinstance(value, dict):
                differences += compare_config_dicts(value, used_config[key], path + [key])
            else:
                if value != used_config[key]:
                    differences.append(path + [key])
        else:
            #key does not exist in used config
            sys.stderr.write('Warning - config key {} not set in used config\n'.format(':'.join(path+[key])))
    return differences

def read_config_file(print_config, path):
    """
    Helper function to read a config file, if it exists.
    Always returns a dict, but it may be empty.
    """
    this_config = {}
    if os.path.isfile(path):
        if print_config:
            sys.stderr.write('Reading additional configuration file: {}\n'.format(path))
        with open(path, 'r') as config_file:
            this_config = yaml.safe_load(config_file.read())
        
        #make sure we always get an empty dict, yaml.safe_load may give us None for an empty file
        if this_config is None:
            this_config = {}

    return this_config

def main(argv = None):
    """
    Main function for the ``amplimap`` executable. This function:

    - parses command line arguments
    - reads, merges and checks each of these config files, if they exist:
        + ``config_default.yaml`` in the amplimap package
        + ``/etc/amplimap/VERSION/config.yaml`` (where VERSION is the amplimap version)
        + ``$AMPLIMAP_CONFIG``
        + ``config.yaml`` in the working directory
    - checks for an existing analysis directory (and compares the amplimap version used to create it)
    - adds its own parent directory to the config file (to be inserted back into the python path inside Snakemake)
    - creates an analysis directory
    - writes ``config_used.yaml`` to the new analysis directory
    - creates a ``cluster_log`` directory (if running in cluster mode)
    - launches Snakemake, using the amplimap Snakefile, ``config_used.yaml`` as the config file and cluster parameters as specified in the command line arguments and config.
    """
    try:
        basedir = os.path.dirname(os.path.realpath(__file__))
        #sys.stderr.write('Called with arguments: "{}"\n'.format('" "'.join(sys.argv)))
        
        #parse the arguments, which will be available as properties of args (e.g. args.probe)
        parser = argparse.ArgumentParser(
            description = "amplimap v{} - amplicon mapping and analysis pipeline".format(__version__),
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
        #specify parameters
        parser.add_argument("-v", "--version", help="print version and exit", action="store_true")
        parser.add_argument("--basedir", help="print basedir and exit", action="store_true")
        parser.add_argument("--print-config", help="print configuration (including global and local settings) and exit", action="store_true")
        parser.add_argument("-r", "--run", help="actually run (will perform a dry run otherwise)", action="store_true")
        parser.add_argument("--resume", help="resume analysis in existing analysis directory", action="store_true")
        parser.add_argument("--cluster", help="specify a cluster type defined in your configuration files to run jobs on cluster.")
        parser.add_argument("--skip-file-check", help="skip check for changes in input files when resuming (not recommended)", action="store_true")
        parser.add_argument("--unlock", help="unlock working directory (Snakemake parameter)", action="store_true")
        parser.add_argument("--working-directory", help="path to the working directory", default=".")
        parser.add_argument("--ncores", help="number of local cores to run in parallel (only applies if --cluster is NOT set)", default=1, type=int)
        parser.add_argument("--njobs", help="number of cluster jobs to run in parallel (only applies if --cluster is set)", default=10, type=int)
        parser.add_argument("--latency-wait", help="How long to wait for output files to appear after job completes. Increase this if you get errors about missing output files. (Snakemake parameter)", default=5, type=int)
        parser.add_argument("--debug", help="debug mode", action="store_true")
        #parser.add_argument("--debug-dag", help="debug DAG", action="store_true")
        parser.add_argument("TARGET", help="targets to run (eg. pileups variants coverages)", nargs="*")
        if argv is None:
            args = parser.parse_args()
        else:
            args = parser.parse_args(argv)

        if args.debug:
            print('Incoming argv: {}'.format(str(argv)))
            print('Targets: {}'.format(str(args.TARGET)))

        if args.version:
            print('{} {}'.format(__title__, __version__))
            return 0

        if args.basedir:
            print(basedir)
            return 0

        #read base config to know which parameters etc are allowed
        default_config = read_config_file(args.print_config, os.path.join(basedir, 'config_default.yaml'))
        if not default_config:
            raise Exception('config_default.yaml file missing!')

        #add undocumented config keys to make sure these don't raise an error
        for key in ['include_gbrowse_links', 'include_exon_distance', 'include_score']:
            if not key in default_config['annotate']:
                default_config['annotate'][key] = False

        #override with data from /etc/amplimap, if exists
        etc_config = read_config_file(args.print_config, '/etc/amplimap/%s/config.yaml' % __version__)

        #override with data from $AMPLIMAP_CONFIG, if exists
        env_config = {}
        try:
            env_config = read_config_file(args.print_config, os.environ['AMPLIMAP_CONFIG'])
        except KeyError:
            pass

        #read local config
        local_config = read_config_file(args.print_config, os.path.join(args.working_directory, 'config.yaml'))
        if not local_config:
            if args.print_config:
                sys.stderr.write('No local config.yaml found, using default configuration.\n')

        #merge configs together
        config = default_config
        for my_config in [etc_config, env_config, local_config]:
            #check that all settings actually exist
            differences = check_config_keys(default_config, my_config)
            if len(differences) > 0:
                sys.stderr.write('Your configuration file(s) contain unknown or invalid settings:\n')
                for diff in differences:
                    sys.stderr.write('\t- {}\n'.format(':'.join(diff)))
                sys.stderr.write('Please check their spelling and location and try again.\n')
                return 1

            snakemake.utils.update_config(config, my_config)

        if args.print_config:
            yaml.dump(config, sys.stdout, default_flow_style=False)
            return 0

        #do some basic checks
        assert os.path.isdir(args.working_directory), 'working directory does not exist'
        assert os.path.isdir(os.path.join(args.working_directory, 'reads_in')) \
            or os.path.isdir(os.path.join(args.working_directory, 'bams_in')) \
            or os.path.isdir(os.path.join(args.working_directory, 'tagged_bams_in')) \
            or os.path.isdir(os.path.join(args.working_directory, 'unmapped_bams_in')), 'reads_in/, bams_in/ or unmapped_bams_in/ directory missing'
        # assert os.path.isfile(os.path.join(args.working_directory, 'probes.csv')) \
        #     or os.path.isfile(os.path.join(args.working_directory, 'probes_mipgen.csv')) \
        #     or os.path.isfile(os.path.join(args.working_directory, 'probes_heatseq.tsv')), 'probes.csv, probes_mipgen.csv, or probes_heatseq.tsv file missing'

        #check some basic settings
        aligners = ['naive', 'bwa', 'bowtie2', 'star'] #allowed values for the aligner
        if not config['align']['aligner'] in aligners:
            raise Exception('align: aligner must be one of {}!'.format(','.join(aligners)))

        callers = ['gatk', 'platypus'] #allowed values for the variant caller
        if not config['variants']['caller'] in callers:
            raise Exception('variants: caller must be one of {}!'.format(','.join(callers)))

        if config['parse_reads']['quality_trim_threshold'] != False:
            if not isinstance(config['parse_reads']['quality_trim_threshold'], float):
                raise Exception('quality_trim_threshold must be a decimal number!')
            if not config['parse_reads']['quality_trim_threshold'] > 0 and config['parse_reads']['quality_trim_threshold'] < 1:
                raise Exception('quality_trim_threshold must be either "false" or above 0 and below 1!')

        if not (config['general']['umi_min_consensus_percentage'] >= 0 and config['general']['umi_min_consensus_percentage'] <= 100):
            raise Exception('umi_min_consensus_percentage must be between 0 and 100')

        if not (config['parse_reads']['min_percentage_good'] >= 0 and config['parse_reads']['min_percentage_good'] <= 100):
            raise Exception('min_percentage_good must be between 0 and 100')

        if not (config['parse_reads']['umi_one'] >= 0 and config['parse_reads']['umi_two'] >= 0):
            raise Exception('umi_one and umi_two must be 0 or greater')

        if config['annotate']['annovar']['protocols'].count(',') != config['annotate']['annovar']['operations'].count(','):
            raise Exception('The number of comma-separated protocols and operations under `annotate: annovar:` must match!')

        #if we don't have UMIs (either on reads or as bam tag) we definitely have to ignore them
        #this makes it possible to "auto-detect" whether we need to ignore_umis or not
        if (config['parse_reads']['umi_one'] + config['parse_reads']['umi_two'] == 0) and config['general']['umi_tag_name'] == "":
            config['general']['ignore_umis'] = True

        #check we have proper paths
        if not config['general']['genome_name'] in config['paths'] or not isinstance(config['paths'][config['general']['genome_name']], dict):
            raise Exception('Could not find list of paths for genome_name: "{}". Please add the paths to your default configuration or your local config.yaml file.'.format(config['general']['genome_name']))

        for name, path in config['paths'][config['general']['genome_name']].items():
            if path.startswith('/PATH/TO/'):
                raise Exception('Path for {} reference is set to {}, which is probably incorrect. Please set the correct path in your default configuration or your local config.yaml file, or leave it empty.'.format(
                    name, path))

        #check input files
        sys.stderr.write('Checking input files...\n')
        if os.path.isfile(os.path.join(args.working_directory, 'probes.csv')):
            read_new_probe_design(os.path.join(args.working_directory, 'probes.csv'), reference_type = 'genome')
        if os.path.isfile(os.path.join(args.working_directory, 'probes_mipgen.csv')):
            process_probe_design(read_and_convert_mipgen_probes(os.path.join(args.working_directory, 'probes_mipgen.csv')))
        if os.path.isfile(os.path.join(args.working_directory, 'probes_heatseq.tsv')):
            process_probe_design(read_and_convert_heatseq_probes(os.path.join(args.working_directory, 'probes_heatseq.tsv')))
        if os.path.isfile(os.path.join(args.working_directory, 'targets.csv')):
            #note: this will fail on overlapping targets
            read_targets(os.path.join(args.working_directory, 'targets.csv'), check_overlaps=True, reference_type = 'genome', file_type = 'csv')
        if os.path.isfile(os.path.join(args.working_directory, 'targets.bed')):
            #note: this will fail on overlapping targets
            read_targets(os.path.join(args.working_directory, 'targets.bed'), check_overlaps=True, reference_type = 'genome', file_type = 'bed')
        if os.path.isfile(os.path.join(args.working_directory, 'snps.txt')):
            read_snps_txt(os.path.join(args.working_directory, 'snps.txt'), reference_type = 'genome')

        #this will be used to (very hackily) make sure amplimap can be imported as amplimap.xxx
        #by adding the parent dir to the top of sys.path in the Snakefile
        config['general']['amplimap_parent_dir'] = os.path.dirname(basedir)

        #check if analysis dir exists already
        analysis_dir = os.path.join(args.working_directory, 'analysis')
        configfile = os.path.join(analysis_dir, 'config_used.yaml')
        used_versions_path = os.path.join(analysis_dir, 'versions.yaml')

        #the analysis dir may exist just because we did a dry run, but once the versions exist we actually executed snakemake!
        if os.path.exists(analysis_dir) and os.path.exists(used_versions_path):
            if not args.resume:
                raise Exception('An analysis directory already exists. Please rename it or set --resume to reuse it and possibly overwrite existing files.')
            else:
                #check version
                if os.path.isfile(used_versions_path):
                    with open(used_versions_path, 'r') as used_versions_file:
                        used_versions = yaml.safe_load(used_versions_file.read())
                        if used_versions['_amplimap'] != str(__version__):
                            sys.stderr.write('This analysis was performed with {} {} but this is {} {}!\n\n'.format(__title__, used_versions['_amplimap'], __title__, __version__))
                            sys.stderr.write('Please use the correct version of {} or start a new analysis.\n'.format(__title__))
                            return 1
                        else:
                            sys.stderr.write('{} version checked.\n'.format(__title__))

                #check used config file
                if os.path.isfile(configfile):
                    with open(configfile, 'r') as used_config_file:
                        used_config = yaml.safe_load(used_config_file.read())
                        differences = compare_config_dicts(config, used_config)
                        if len(differences) > 0:
                            sys.stderr.write('config_used.yaml in analysis directory differs from current config.yaml in working directory! Please rename or delete the old analysis directory to restart analysis with the new configuration.\n')
                            sys.stderr.write('Different settings:\n')
                            for diff in differences:
                                sys.stderr.write('\t- {}\n'.format(':'.join(diff)))
                            return 1
                        else:
                            sys.stderr.write('Config files checked.\n')

                #check hashes of input files
                if not args.skip_file_check:
                    used_file_hashes_path = os.path.join(analysis_dir, 'file_hashes.yaml')
                    if os.path.isfile(used_file_hashes_path):
                        with open(used_file_hashes_path, 'r') as used_file_hashes_file:
                            used_file_hashes = yaml.safe_load(used_file_hashes_file.read())

                            from .reader import get_file_hashes
                            for fn, current_hash in get_file_hashes(args.working_directory).items():
                                if used_file_hashes[fn] != current_hash:
                                    sys.stderr.write('File {} seems to have changed since the last run!\n\n'.format(fn))
                                    sys.stderr.write('To ensure consistent results, you should rename or delete the old analysis directory and start a new analysis.\n')
                                    sys.stderr.write('To ignore this error, add the --skip-file-check parameter.\n')
                                    return 1
                            sys.stderr.write('Input files checked.\n')
                else:
                    sys.stderr.write('Warning: Skipping input file check.\n')

        #ensure analysis dir exists now
        try:
            os.makedirs(analysis_dir)
        except OSError as e:
            pass

        #write config to analysis directory, and then use that for snakemake
        with open(configfile, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)

        #set up cluster commands
        cluster_command_nosync = None
        cluster_command_sync = None

        if args.cluster:
            if args.cluster in config['clusters'] and isinstance(config['clusters'][args.cluster], dict):
                if 'command_sync' in config['clusters'][args.cluster]:
                    cluster_command_sync = config['clusters'][args.cluster]['command_sync']
                elif 'command_nosync' in config['clusters'][args.cluster]:
                    cluster_command_nosync = config['clusters'][args.cluster]['command_nosync']
                else:
                    raise Exception('Invalid cluster configuration -- need either command_sync or command_nosync for: {}'.format(args.cluster))
            else:
                raise Exception('Cluster type not found in config: {}'.format(args.cluster))

            sys.stderr.write('Running in cluster mode {} with {} parallel jobs\n'.format(args.cluster, args.njobs))
            sys.stderr.write('cluster_command_nosync={}\n'.format(cluster_command_nosync))
            sys.stderr.write('cluster_command_sync={}\n'.format(cluster_command_sync))

            #make sure cluster log directory exists (this assumed the cluster command is using this as a parameter)
            cluster_logs = os.path.join(args.working_directory, 'cluster_log')
            try:
                os.makedirs(cluster_logs)
            except OSError as e:
                pass
            sys.stderr.write('Will write cluster logs to: {}\n'.format(cluster_logs))
        else:
            sys.stderr.write('Running locally with {} cores\n'.format(args.ncores))

        success = snakemake.snakemake(
            snakefile = os.path.join(basedir, "Snakefile"),
            configfile = configfile,
            cores = args.ncores, #ignored if cluster
            nodes = args.njobs, #ignored if not cluster
            workdir = args.working_directory,
            targets = args.TARGET,
            dryrun = not args.run,
            cluster = cluster_command_nosync,
            cluster_sync = cluster_command_sync,
            jobname = "{}.{{rulename}}.{{jobid}}.sh".format(__title__),
            unlock = args.unlock,
            latency_wait = args.latency_wait,
            #debug_dag = args.debug_dag,
            )

        if success:
            if args.unlock:
                sys.stderr.write('Unlocked working directory. Run without --unlock to start.\n')
            elif not args.run:
                sys.stderr.write('{} {} dry run successful. Set --run to run!\n'.format(__title__,  __version__))
            else:
                sys.stderr.write('{} {} finished!\n'.format(__title__, __version__))
            return 0
        else:
            if args.cluster:
                sys.stderr.write('{} {} failed! Please see output above or cluster logs for details.\n'.format(__title__,  __version__))
            else:
                sys.stderr.write('{} {} failed! Please see output above for details.\n'.format(__title__,  __version__))

            return 1
    except AmplimapReaderException as e:
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.stderr.write(str(e))
        sys.stderr.write('{} {} failed!\n'.format(__title__, __version__))
        return 2
    except Exception as e:
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.stderr.write('\nERROR: {}\n\n'.format(e))
        sys.stderr.write('{} {} failed!\n'.format(__title__, __version__))
        return 1

if __name__ == '__main__':
    sys.exit(main())

