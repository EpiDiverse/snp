# EpiDiverse-SNP Output
This document describes the output produced by the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Pre-processing](#pre-processing) -Sample filtering and parsing the samplesheet
* [Masking](#masking) -Combining all samples into single files
* [khmer and kWIP](#khmer-and-kwip) -Intersecting methylated positions based on DMPs and/or DMRs
* [Variant calling](#variant-calling) -Calculating average methylation values for given DMRs
* [Post-processing](#post-processing) -Filtering and merging input variant call file(s)

### Output Directory Structure
![Output Directory Structure](/docs/images/directory.png)

## Pre-processing
The pipeline requires individual bedGraphs (in each specified methylation context) and a samplesheet with corresponding sample names, environmental trait values, and covariate values in order to run. The samplesheet is processed into the format required to run GEM, and individual sample bedGraphs are filtered according to coverage. If DMP or DMR comparisons are given then these will be filtered for a user-specified significance threshold on positions / regions prior to downstream analysis.

**Output directory: `ewas/input/`**

* `cov.txt`
* `env.txt`
* `gxe.txt`
  * **NB:** Only saved if GxE model is enabled during the pipeline run.


## Bedtools unionbedg
Following sample pre-processing, the entire collection of each input type is merged into single files per each methylation context. From *.bedGraph files each position denotes the methylation value for each input sample, and all files denote the presence/absence of a given position/region for individual samples/comparisons by the use of "NA".

**Output directory: `ewas/input/bed`**

* `{CpG,CHG,CHH}.bedGraph.bed`

```
chrom	start	end	sample1	sample2	sample3	sample4	sample5	sample6	sample7	sample8	sample9
scaffold_53	390	391	0.50	1.00	NA	NA	0.33	1.00	1.00	0.50	0.00
scaffold_53	392	393	NA	1.00	NA	NA	NA	0.00	NA	NA	NA
scaffold_53	581	582	0.66	0.75	1.00	1.00	0.66	1.00	0.66	0.75	1.00
scaffold_53	583	584	0.87	1.00	1.00	1.00	0.75	0.83	0.77	0.71	1.00
scaffold_53	671	672	0.50	1.00	1.00	1.00	1.00	1.00	1.00	1.00	1.00
scaffold_53	673	674	0.87	0.93	0.66	0.50	0.85	0.83	0.90	1.00	1.00
...
```

* `{CpG,CHG,CHH}.DMPs.bed`
  * **NB:** Only saved if DMPs are given during the pipeline run.
* `{CpG,CHG,CHH}.DMRs.bed`
  * **NB:** Only saved if DMRs are given during the pipeline run.

```
chrom   start   end     g1_vs_g2    g1_vs_g3    g2_vs_g3
scaffold_53     166683  166807  NA      NA      0.037
scaffold_53     227390  227644  NA      NA      0.006
scaffold_53     309090  309149  NA      0.000   NA
scaffold_53     309149  309180  0.017   0.000   NA
scaffold_53     309180  309262  0.017   NA      NA
scaffold_53     309535  309715  NA      0.000   0.000
...
```


## GxE model
The GxE model tests for the interaction between methQTLs and the environmental trait. As the matrix of methylated positions vs SNPs is orders of magnitude larger than methylated positions alone, this analysis is divided among individual scaffolds and combined at the end for FDR calculation.

**Output directory: `ewas/{positions,regions}/GxE`**

* `*.txt`
  * The full results from GEM GxE model output
* `*.filtered_*_FDR.txt`
  * The results from GEM GxE model filtered by FDR threshold

```
cpg                snp     beta        stats     pvalue        FDR
MA_1063600_3380    SNP962  0.17808859  42.28204  1.482360e-111 1.482360e-106
MA_130823_2338     SNP700 -0.21534752 -18.43573  1.761554e-47  8.807769e-43
MA_124616_3396     SNP578 -0.15171656 -16.70323  9.169281e-42  3.056427e-37
MA_659042_4763     SNP690  0.10567235  13.47239  5.237893e-31  1.309473e-26
MA_101037_18934    SNP589  0.07781375  13.07099  1.112935e-29  2.225870e-25
MA_45879_4444      SNP703  0.13979006  12.55871  5.390763e-28  8.984606e-24
...
```

* `*/*.png`
  * Plots for the top K most significant interactions and the associations with the environmental trait for major allele homozygote (AA), heterozygote (AB) and minor allele homozygote (BB) for all SNPs across all samples.

<img align="center" alt="Plot for a single interaction of SNP and methylated position" src="images/kplot.png">


## Pipeline Info
Nextflow has several built-in reporting tools that give information about the pipeline run.

**Output directory: `template/`**

* `dag.svg`
  * DAG graph giving a diagrammatic view of the pipeline run.
  * NB: If [Graphviz](http://www.graphviz.org/) was not installed when running the pipeline, this file will be in [DOT format](http://www.graphviz.org/content/dot-language) instead of SVG.
* `report.html`
  * Nextflow report describing parameters, computational resource usage and task bash commands used.
* `timeline.html`
  * A waterfall timeline plot showing the running times of the workflow tasks.
* `trace.txt`
  * A text file with machine-readable statistics about every task executed in the pipeline.
