#!/usr/bin/env python
import os
import sys

import snakemake
import argparse
import yaml

import readline
import collections

from .version import __title__, __version__

class Question:
    def __init__(self, question, default = None, options = None):
        self.question = question
        self.default = default
        self.options = options

        if self.default and self.options:
            assert self.default in self.options

    def ask(self):
        if self.options is not None:
            if self.default is not None:
                question = '{} ([{}], {})'.format(self.question, self.default, ', '.join([o for o in self.options if o != default]))
            else:
                question = '{} ({})'.format(self.question, ', '.join(self.options))
        elif self.default is not None:
            question = '{} [{}]'.format(self.question, self.default)

        answer = ''
        while len(answer) == 0:
            answer = input(question)
            if len(answer) == 0:
                #use the default if we have no answer
                if self.default is not None:
                    answer = self.default
                else:
                    sys.stdout.write('Please provide a value!\n')
            else:
                #check valid answer if we have options
                if self.options is not None:
                    if not answer in self.options:
                        sys.stdout.write('Please enter one of: {}\n'.format(', '.join(self.options)))
                        answer = ''

        return answer

def main(argv = None):
    """
    Run amplimap setup wizard.
    """
    try:
        basedir = os.path.dirname(os.path.realpath(__file__))
        
        #parse the arguments, which will be available as properties of args (e.g. args.probe)
        parser = argparse.ArgumentParser(
            description = "amplimap v{} setup wizard".format(__version__),
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)

        #specify parameters
        parser.add_argument("operation", help="type of setup to perform: paths / indices", action="store_true")
        parser.add_argument("--debug", help="debug mode", action="store_true")
        if argv is None:
            args = parser.parse_args()
        else:
            args = parser.parse_args(argv)

        raise Exception('Not implemented yet.')

        if args.operation == 'paths':
            #TODO:
            #- ask aligner, caller
            #- ask about modules:
            #- check path for:
            #aligner, caller, bedtools, samtools (error if not found)
            #annovar, picard, bcftools (warning if not found)

            answers = collections.ordereddict()
            for key, q in questions.items():
                answers[key] = q.ask()

            # print('Your answers:')
            # for key, answer in answers.items():
            #     print('{}: {}'.format(
            #         key, answer
            #     ))

            tools_required = [
                'bedtools', 'samtools'
            ]

            tools_for_variants = [
                'annovar', 'bcftools'
            ]

            if answers['aligner'] == 'bwa':
                tools_required.append('bwa')
            elif answers['aligner'] == 'bowtie2':
                tools_required.append('bwa')
            elif answers['aligner'] == 'star':
                tools_required.append('bwa')
        elif args.operation == 'indices':
            #TODO:
            #- ask path to reference genome
            #- build fasta index, default aligner index
            #- ask to enter name of other index to create (bwa, bowtie2, star)
            pass
        else:
            raise Exception('Please specify a valid operation!')
    except KeyboardInterruptException:
        sys.stderr.write('Interrupted.\n')
        return 1
    except EOFError:
        sys.stderr.write('Aborted.\n')
        return 1
    except Exception as e:
        sys.stderr.write('\nERROR: {}\n\n'.format(e))
        sys.stderr.write('{} {} failed!\n'.format(__title__, __version__))
        return 1

if __name__ == '__main__':
    sys.exit(main())
