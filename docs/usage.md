# EpiDiverse-EWAS Usage
This document describes the parameter options used by the pipeline.

* [Running the pipeline](#running-the-pipeline)
* [Inputs and outputs](#inputs-and-outputs)
    * [`--input`](#--input-arg-required)
    * [`--reference`](#--reference-arg)
    * [`--output`](#--output-arg)
* [Modifiers](#modifiers)
    * [`--variants`](#--variants)
    * [`--clusters`](#--clusters)
* [Variant calling](#variant-calling)
    * [`--coverage`](#--coverage-arg)
    * [`--ploidy`](#--ploidy-arg) 
    * [`--regions`](#--regions-arg)
* [Additional Parameters](#Additional-Parameters)
    * [`--debug`](#--debug)
    * [`--version`](#--version)
    * [`--help`](#--help)
* [Software dependencies](#software-dependencies)
    * [`-profile`](#-profile)
    * [`-with-conda`](#-with-conda)
    * [`-with-docker`](#-with-docker)
    * [`-with-singularity`](#-with-singularity)
* [Other command line parameters](#other-command-line-parameters)
    * [`-work-dir`](#-work-dir)
    * [`-params-file`](#-params-file)
    * [`-config`](#-config)
    * [`-resume`](#-resume)
    * [`-name`](#-name)

## Workflow

![EpiDiverse/snp Workflow](/docs/images/workflow.png)

## Running the pipeline
The main command for running the pipeline is as follows:

```bash
nextflow run epidiverse/snp [OPTIONS]
```

Note that the pipeline will create files in your working directory:

```bash
work/           # Directory containing the nextflow working files
snps/           # Finished results (configurable, see below)
.nextflow.log   # Log file from Nextflow
.nextflow/      # Nextflow cache and history information
```

## Inputs and Outputs

### `--input <ARG>` [REQUIRED]
Specify input path for the directory containing outputs from the WGBS pipeline. The pipeline searches for bam files in '\*/{sample_name}/{sample_name}.bam' format.

### `--reference <ARG>`
Specify the path to the input reference genome file in fasta format. REQUIRED for the variant calling aspect of the pipeline, along with a corresponding fasta index *.fai file in the same location.

### `--output <ARG>`
A string that will be used as the name for the output results directory, which will be generated in the working directory. [default: snps]


## Modifiers

If neither of the following two options are specified, the pipeline will run only the double-masking procedure and provide output bam files. 

### `--variants`
Run pipeline in variant calling mode. [default: off]

### `--clusters`
Run pipeline in clustering mode. [default: off]


## Variant Calling 

### `--coverage <ARG>`
Require at least this coverage to process a variant site. [default: 0]

### `--ploidy <ARG>`
Specify the expected ploidy count for the genome [default: 2]

### `--regions <ARG>`
Variant calling with Freebayes parallel will split bam files into chunks of alignments over reference regions of N bases as defined by this parameter [default: 100000]


## Additional Parameters

### `--debug`
Specify in order to prevent Nextflow from clearing the work dir cache following a successful pipeline completion. [default: off]

### `--version`
When called with `nextflow run epidiverse/snp --version` this will display the pipeline version and quit.

### `--help`
When called with `nextflow run epidiverse/snp --help` this will display the parameter options and quit.

## Software Dependencies

There are different ways to provide the required software dependencies for the pipeline. The recommended method is to use the Conda, Docker or Singularity profiles as provided by the pipeline. 

### `-profile`
Use this parameter to choose a preset configuration profile. See the [installation documentation](https://app.gitbook.com/@epidiverse/s/project/epidiverse-pipelines/installation) for more information about profiles.

Profiles available with the pipeline are:

* `standard`
    * The default profile, used if `-profile` is not specified.
    * Uses sensible resource allocation for , runs using the `local` executor (native system calls) and expects all software to be installed and available on the `$PATH`.
    * This profile is mainly designed to be used as a starting point for other configurations and is inherited by most of the other profiles below.
* `conda`
    * Builds a conda environment from the environment.yml file provided by the pipeline
    * Requires conda to be installed on your system.
* `docker`
    * Launches a docker image pulled from epidiverse/snp
    * Requires docker to be installed on your system. 
* `singularity`
    * Launches a singularity image pulled from epidiverse/snp
    * Requires singularity to be installed on your system.
* `epi|diverse`
    * Designed to be used on the [EpiDiverse](http://epidiverse.eu/) clusters `epi` or `diverse`
    * Launches jobs using the `SLURM` executor.
    * Uses pre-built conda environments to provide all software requirements.
* `custom`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config for process resource allocation.

If you wish to provide your own package containers it is possible to do so by setting the `standard` or `custom` profile, and then providing your custom package with the command line flags below. These are not required with the the other profiles.

### `-with-conda <ARG>`
Flag to enable conda. You can provide either a pre-built environment or a *.yml file.

### `-with-docker <ARG>`
Flag to enable docker. The image will automatically be pulled from Dockerhub.

### `-with-singularity <ARG>`
Flag to enable use of singularity. The image will automatically be pulled from the internet. If running offline, follow the option with the full path to the image file.

## Other command line parameters

### `-work-dir <ARG>`
Specify the path to a custom work directory for the pipeline to run with (eg. on a scratch directory)

### `-params-file <ARG>`
Provide a file with specified parameters to avoid typing them out on the command line. This is useful for carrying out repeated analyses. A template params file [`assets/params.txt`](../assets/params.txt) has been made available in the pipeline repository.

### `-config <ARG>`
Provide a custom config file for adapting the pipeline to run on your own computing infrastructure. A template config file [`assets/custom.config`](../assets/custom.config) has been made available in the pipeline repository. This file can be used as a boilerplate for building your own custom config.

### `-resume [<ARG>]`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. Give a specific pipeline name as an argument to resume it, otherwise Nextflow will resume the most recent. NOTE: This will not work if the specified run finished successfully and the cache was automatically cleared. (see: [`--debug`](#--debug))

### `-name <ARG>`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
