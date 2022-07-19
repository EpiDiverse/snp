[<img width="200" align="right" src="docs/images/euflagbetter.jpg">](https://ec.europa.eu/programmes/horizon2020/en)
[<img width="200" align="right" src="docs/images/epidiverse-logo.jpg">](https://epidiverse.eu)

EpiDiverse-SNP Pipeline
========================

[![Nextflow](https://img.shields.io/badge/nextflow-20.07.1-45818e.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-45818e.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/cloud/automated/epidiverse/snp.svg)](https://hub.docker.com/r/epidiverse/snp)
[![Release](https://img.shields.io/github/v/release/epidiverse/snp.svg?colorB=45818e)](https://github.com/EpiDiverse/snp/releases/latest)
[![Publication](https://img.shields.io/badge/Publication-in%20prep.-45818e.svg?colorB=45818e&style=popout)]()
[![Twitter](https://img.shields.io/twitter/follow/epidiverse?style=social)](https://twitter.com/intent/follow?screen_name=epidiverse)

**EpiDiverse/snp** is a bioinformatics analysis pipeline for calling single nucleotide polymorphism variants from bisulfite sequencing data and/or for clustering of eg. environmental plant samples according to their methylation profiles while masking the genomic variation.

The workflow pre-processes a collection of bam files from the [EpiDiverse/WGBS](https://github.org/epidiverse/wgbs) pipeline using [samtools](https://github.com/samtools/samtools), then masks genomic and/or bisulfite variation relative to the reference using custom scripts. Genomic masked alignments are then extracted into fastq format and tested for kmer diversity using [kWIP](https://kwip.readthedocs.io/en/latest/) for clustering groups. Bisulfite-masked alignments are taken forward for variant calling using a combination of [Freebayes](https://github.com/ekg/freebayes) and post-call filtering with [bcftools](https://github.com/samtools/bcftools).

> See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/)

2. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

3. Start running your own analysis!

```bash
NXF_VER=20.07.1 nextflow run epidiverse/snp -profile <docker|singularity|conda> \
--input /path/to/wgbs/bam --reference /path/to/reference.fa
```

> See the [usage documentation](docs/usage.md) for all of the available options when running the pipeline.

### Test data

A minimal example dataset for testing purposes can be found in the [EpiDiverse/datasets](https://github.com/EpiDiverse/datasets/tree/snp) repository. You can either download the files manually and run the pipeline above as intended, or you can directly run the pipeline using the test profile option which will automatically download the data for you:

```bash
NXF_VER=20.07.1 nextflow run epidiverse/snp -profile test,<docker|singularity|conda>
```

## Wiki Documentation

The EpiDiverse/snp pipeline is part of the [EpiDiverse Toolkit](https://epidiverse.gitbook.io/project/-MfxkdBDZggX_vc_sG5l/epidiverse-pipelines/best-practice-pipelines), a best practice suite of tools intended for the study of [Ecological Plant Epigenetics](https://epidiverse.gitbook.io/project/-MfxkdBDZggX_vc_sG5l/). Links to general guidelines and pipeline-specific documentation can be found below:

1. [Installation](https://epidiverse.gitbook.io/project/-MfxkdBDZggX_vc_sG5l/epidiverse-pipelines/installation#1-install-nextflow)
2. Pipeline configuration
    * [Local installation](https://epidiverse.gitbook.io/project/-MfxkdBDZggX_vc_sG5l/epidiverse-pipelines/installation#2-install-the-pipeline)
    * [Adding your own system config](https://epidiverse.gitbook.io/project/-MfxkdBDZggX_vc_sG5l/epidiverse-pipelines/installation#3-pipeline-configuration)
    * [EpiDiverse infrastructure](https://epidiverse.gitbook.io/project/-MfxkdBDZggX_vc_sG5l/epidiverse-pipelines/installation#appendices)
3. [Running the pipeline](docs/usage.md)
4. [Understanding the results](docs/output.md)
5. [Runtime and memory usage guidelines](docs/runtime.md)
6. [Troubleshooting](https://epidiverse.gitbook.io/project/-MfxkdBDZggX_vc_sG5l/epidiverse-pipelines/troubleshooting)

### Credits

These scripts were originally written for use by the [EpiDiverse European Training Network](https://epidiverse.eu/), by Adam Nunn ([@bio15anu](https://github.com/bio15anu)).

This project has received funding from the European Union’s Horizon 2020 research and innovation
programme under the Marie Skłodowska-Curie grant agreement No 764965

## Citation

If you use epidiverse/snp for your analysis, please cite it using the following doi: <placeholder>