# PacBio Subread Stats Pipeline #

Computes summary statistics for unaligned PacBio reads including total throughput and subread length statistics.

Subread statistics come in two categories and are prefixed with "SUB" or "ROI".

* SUB: Calculated on all subreads.
  * Estimates total throughput
  * May be skewed by short inserts that were sequenced several times (generated more than one subread), which creates
    "mutual-information" in the sequencing data set. 
* ROI: Calculated on the longest subread per read.
  * Reads of Insert (ROI)
  * A better estimation of the input library


## Run ##

### Input file ###

The only required input file is a list of sequence data files. These may be ".bax.h5" files from RS II, or
".subreads.bam" files from Sequel.

This pipeline is up to date as of Sequel 2.2 (1M chip). It is not sensitive to the sequencing platform or chemistry
as long as the structure of the input files remains consistent.

The file name should be "SAMPLE.fofn" where "SAMPLE" is the sample name. If it follows this convention, then the
pipeline can get the sample name from the FOFN file name instead of requiring that it be specified. Using sample names
allows several samples to be run in the same run directory.

This file name is provided to the pipeline after `--config` on the command line (see below).

### Execute ###

Define a variable that gives the full path to the pipeline code, which is the directory that contains `Snakefile`
and this `README.md` file. The pipeline itself does not use the variable, but commands in this README will.

Example:
`PIPELINE_DIR=/net/eichler/vol27/projects/structural_variation/nobackups/pipelines/subread_stats/201905`

Load required modules (may work with later versions of these modules):
```
module load miniconda/4.5.12
```

Run distributed:
`mkdir -p log; snakemake -s ${PIPELINE_DIR}/Snakefile --ri -j 30 -k --jobname "{rulename}.{jobid}" --drmaa " -V -cwd -e ./log -o ./log -pe serial {cluster.cpu} -l mfree={cluster.mem} -l h_rt={cluster.rt} -l gpfsstate=0 -j y -w n -S /bin/bash" -w 60 -u ${PIPELINE_DIR}/config/cluster.json --config fofn=/path/to/SAMPLE.fofn`

...where "SAMPLE.fofn" is the input file name.


#### Specifying the sample name ####

The sample name may be explicitly provided on the command line by specifying `sample=SAMPLE` on the command line
anywhere after `--config`.
