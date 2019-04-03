Quickstart
----------

-  Create a new working directory
-  Put the :ref:`probes-csv` file into the working directory
   (:download:`download probes.csv example <examples/probes.csv>`)
-  Put the :ref:`targets-csv` file into the working directory
   (:download:`download targets.csv example <examples/targets.csv>`)
-  Create a subdirectory called :ref:`reads-in` in the working directory
   and copy your zipped sample FASTQ files (ending in ``.fastq.gz``)
   into it
-  Open a terminal and change into the working directory
-  Do a dry run: ``amplimap``
-  Run amplimap: ``amplimap --run``

Pipeline overview
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Clip off UMIs (if any), identify probe by matching known probe arms
   to the beginning of both reads, trim off probe arms, optionally trim
   low-quality bases (→ parsed FASTQ files)
2. Align parsed reads (without arms) to reference genome (→ BAM files)
3. Calculate alignment stats
4. Germline variant calling and annotation (e.g. to call variants in resequencing data):

   1. Calculate coverage across target regions
   2. Call variants on raw reads in target regions
   3. Call variants on UMI deduplicated reads (not consensus) in target
      regions
   4. Annotate coverage and variant tables with sample information

5. Per-basepair pileups (e.g. to find low-frequency somatic mutations):

   1. Calculate consensus pileup for target regions
   2. Calculate consensus pileup for known SNPs (if provided)
