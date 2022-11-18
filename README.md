## Snakemake workflow for single-cell DNA sequencing analysis (scAbsolute, scUnique).

This workflow performs absolute copy number calling in single cell DNA sequencing data as described in the [scAbsolute manuscript
](https://www.biorxiv.org/content/10.1101/2022.11.14.516440v2)

=======
This workflow performs absolute copy number calling in single cell DNA sequencing data as described in TODO
Please make sure to cite our publications if you use this program.
Link to citation

We recommend the following approach.
1. Initially, run the per-cell part of the pipeline (scAbsolute, source code available at
   https://github.com/markowetzlab/scAbsolute).
2. Check results, run qc and exclude outlier cells. Check ploidy, and if necessary rerun script with modified ploidy limits for selected cells. This part of the pipeline requires manual inspection of results.
3. Run the second stage of the pipeline (scUnique, source code available at https://github.com/markowetzlab/scUnique), combining the results from the per-cell analysis and creating a joint analysis of the dataset.

### Dependencies
You will have to install snakemake and singularity/docker (depending on the level of system access you have and where you want to run your analysis)

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
* after qc, exclude samples in config/exclude.tsv

The pipeline will jointly call all samples that are listed in samples.tsv, with the exception of the samples in exclude.tsv.

Given that the workflow has been properly configured, it can be executed as follows.

```bash
~/scDNAseq-workflow$ export SINGULARITY_DOCKER_USERNAME=AWS SINGULARITY_DOCKER_PASSWORD=$(aws ecr get-login-password)
~/scDNAseq-workflow$ snakemake --cores 1 --use-singularity --use-conda results/scale/500/predict/out.rds
```
Please take the time to analyse the data (we have an example qc script available at scripts/run-qc.R)

run second stage
```bash
~/scDNAseq-workflow$ snakemake --cores 1 --use-singularity --use-conda results/scale/500/predict/out.rds
```

Snakemake will automatically detect the main Snakefile in the workflow subfolder and execute the workflow module.

Note, you should not move the output of the first step before running the second step of the pipeline, as this will result in errors in the output file detection.
