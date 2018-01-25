# amplimap
*amplimap* is a pipeline to process and analyse short-read sequencing data from **smMIP-panels** or **PCR-based amplicons**. It was designed to support both **germline variant calling** and the **quantification of somatic mutations** using ultra-high-coverage pileups with and without **UMIs**.

*amplimap* takes fastq or bam files generated from an *Illumina* sequencer (tested with MiSeq/HiSeq/NextSeq), assigns each read pair to a probe/amplicon, trims and aligns the reads and then generates a set of basic statistics. After that, two different types of analysis can be performed:

1. Germline variant calling using *Platypus* or *GATK*, including some further processing to generate a simplified table of potentially pathogenic variants.

2. Generation of per-basepair "pileup" tables with nucleotide counts of each base at each target basepair for each sample, to help estimate somatic mutation frequencies.

Built on top of *Snakemake* and *Python 3*, *amplimap* is entirely automated and can easily be run on a single machine as well as a cluster (eg. LSF, SGD).
