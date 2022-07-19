# EpiDiverse-SNP Runtime and memory usage guidelines
This document describes the default CPUs, RAM, and Time allocation specified for each pipeline process in the default configuration of the pipeline. Configuration was optimised on a HPC cluster with 64 CPUs and 256 Gb RAM, using the collection of plant population datasets provided by EpiDiverse. All values can be adjusted to suit individual needs.

|process|CPUs|RAM / Gb|Time / h|Retries|[errorStrategy](https://www.nextflow.io/docs/latest/process.html#errorstrategy)|
|-------|----|--------|--------|-------|-----------------|
|preprocessing|2|1|3|3|finish|
|masking|8|4|3|3|finish|
|extracting|2|1|1|3|finish|
|khmer|8|4|1|3|finish|
|kwip|8|4|1|3|finish|
|clustering|2|0.5|1|3|finish|
|sorting|2|1|1|3|finish|
|freebayes|8|8|4|3|finish|
|bcftools|2|1|1|3|finish|
|plot_vcfstats|2|1|1|3|ignore|