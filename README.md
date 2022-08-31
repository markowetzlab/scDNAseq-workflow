# Snakemake workflow for single-cell DNA sequencing analysis (scAbsolute, scUnique).

This workflow performs absolute copy number calling in single cell DNA sequencing data as described in TODO

### Checkout repository

```bash
git clone https://github.com/markowetzlab/scDNAseq-workflow.git
```

## Example
A toy dataset of three .bam files can be downloaded from [here](https://drive.google.com/drive/folders/1402zegR4H7tWFluc2el9lyUr9H8rMXX6?usp=sharing). This workflow runs scAbsolute on the toy dataset and merges the outputs into one QDNASeq object named "out.rds".

## Project Structure
We recommend to use the following structure of the project.

    .                               # main folder of the project
    ├── data                        # contains all of the project’s original/raw data
    ├── config                    
    │   ├── config.yaml             # runtime specific configuration
    │   ├── samples.tsv             # metadata and experiement information
    ├── results                     # the output directory
    ├── workflow                    
    │   ├── rules                   # building blocks of the workflows
    │   ├── scripts                 # scripts that are our main analysis
    │   └── Snakefile               # the steps needed that run all the scripts in the project
    ├── LICENSE
    └── README.md



## Usage

To configure this workflow, modify config/ according to your needs. 
* Add samples to config/samples.tsv
* Change binSize in config/config.yaml

The pipeline will jointly call all samples that are defined.

Given that the workflow has been properly configured, it can be executed as follows.

```bash
~/scDNAseq-workflow$ export SINGULARITY_DOCKER_USERNAME=AWS SINGULARITY_DOCKER_PASSWORD=$(aws ecr get-login-password)
~/scDNAseq-workflow$ snakemake --cores 1 --use-singularity --use-conda results/scale/500/predict/out.rds
```

Snakemake will automatically detect the main Snakefile in the workflow subfolder and execute the workflow module.
