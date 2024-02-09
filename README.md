# Snakemake workflow for single-cell DNA sequencing analysis (_scAbsolute_, _scUnique_).

## Overview
This workflow performs absolute copy number calling and detection of recent copy number aberrations (rCNAs) in single-cell DNA sequencing data as described in [scAbsolute](https://www.biorxiv.org/content/10.1101/2022.11.14.516440v2) and [scUnique](missing). 
Please make sure to cite our publications if you use this pipeline:
- scAbsolute:
- scUnique: 

Here is how to use the approach:
1. Initially, run the per-cell part of the pipeline (_scAbsolute_)
2. Inspect results, run quality control and outlier detection across all cells in each sample. If necessary, rerun ploidy analysis with an updated ploidy window. We highly recommend to manually inspect results at this stage as per-sample variations are considerable in our experience.
3. Run the second stage of the pipeline (_scUnique_), combining the results from the per-cell analysis and creating a joint analysis of the dataset including an analysis of rCNAs.

## Dependencies
You must install [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (we support version 8 and above) and [singularity/apptainer](https://apptainer.org/docs/admin/main/installation.html).
We recommend using a cluster environment for the first part of the analysis. _scAbsolute_ is easy to parallelize across cells, and the speedup is linear in the number of CPUs.

## Getting started

```bash
git clone https://github.com/markowetzlab/scDNAseq-workflow.git
```

A toy dataset of one hundred bam files each (PEO1/PEO4) can be downloaded from [here](https://drive.google.com/drive/folders/1402zegR4H7tWFluc2el9lyUr9H8rMXX6?usp=sharing). Move the data to `./data/aligned/PEO1` and `./data/aligned/PEO4`. The workflow can produce rCNAs and copy number profiles for both data sets as described in the _scUnique_ manuscript.
For reasons of simplicity, we only provide a subset of the PEO1/PEO4 data.

### Project Structure
The project has the following structure:

    .                               # main folder of the project
    ├── data                        
    │   └── aligned
    │       └── sample_name         # folder containing all bam files for one sample / to be created by the user
    ├── config                    
    │   ├── config.yaml             # configuration / edit as appropriate for every run
    │   ├── samples.tsv             # file containing sample metadata / to be added by the user
    ├── results                     # the output directory / created by workflow
    ├── workflow                    
        ├── rules                   
        ├── scripts                 # scripts that are our main analysis / edit for additional customization
        └── Snakefile               


## Usage

To configure this workflow, modify config/ according to your needs.
* Add per-sample folder with bam files to data/aligned folder (as with the PEO1/PEO4 folders).
* Create per sample configuration files to the config folder (one file per sample, see PEO1/PEO4.tsv examples).
* Edit variables in config/config.yaml as appropriate.
* after qc, exclude samples in config/exclude.tsv

The pipeline will jointly call all samples that are listed in samples.tsv, with the exception of the samples in exclude.tsv.

Given that the workflow has been properly configured, it can be executed as follows.

```bash
~/scDNAseq-workflow$ export SINGULARITY_DOCKER_USERNAME=AWS SINGULARITY_DOCKER_PASSWORD=$(aws ecr-public get-login-password --region us-east-1)
~/scDNAseq-workflow$ snakemake --cores 1 --use-singularity --use-conda results/scale/500/out.rds
```
Please take the time to analyse the data (we have an example qc script available at scripts/run-qc.R)

run second stage
```bash
~/scDNAseq-workflow$ snakemake --cores 1 --use-singularity --use-conda results/scale/500/out.rds
```

Snakemake will automatically detect the main Snakefile in the workflow subfolder and execute the workflow module.

Note, you should not move the output of the first step before running the second step of the pipeline, as this will result in errors in the output file detection.
