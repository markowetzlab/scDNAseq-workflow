# Snakemake workflow for single-cell DNA sequencing analysis (_scAbsolute_, _scUnique_)

## Overview
This workflow performs absolute copy number calling and detection of recent copy number aberrations (rCNAs) in single-cell DNA sequencing data as described in [*scAbsolute*](https://www.biorxiv.org/content/10.1101/2022.11.14.516440v2) and [*scUnique*](missing).
Please cite the relevant publications if you use this pipeline:
- *scAbsolute*:
```
doi: https://doi.org/10.1101/2022.11.14.516440
```
- *scUnique*:
```
as soon as it's out
```
- if you use *scUnique*, you should also cite [MEDICC2](https://doi.org/10.1186/s13059-022-02794-9):
```
Kaufmann, T.L., Petkovic, M., Watkins, T.B.K. et al. MEDICC2: whole-genome doubling aware copy-number phylogenies for cancer evolution. Genome Biol 23, 241 (2022). https://doi.org/10.1186/s13059-022-02794-9
```

Here is how to use the approach:
1. Initially, run the per-cell part of the pipeline
2. Inspect results and run quality control and outlier detection across all cells in each sample. If necessary, rerun the ploidy analysis with an updated ploidy window. I highly recommend manually inspecting results at this stage, as per-sample variations are considerable in our experience. Per-cell ploidy windows can be specified in the sample files (see config/sample_PEO1.tsv for reference).
3. Run the second stage of the pipeline (_scUnique_), combining the results from the per-cell analysis and creating a joint analysis of the dataset, including an analysis of rCNAs.

## Dependencies
You must install [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (we support version 8 and above) and [singularity/apptainer](https://apptainer.org/docs/admin/main/installation.html).
You must clone the *scAbsolute* and *scUnique* GitHub repositories to your local machine. We demonstrate usage with cloning to the home directory, but it is also possible to use another directory.
```bash
git clone https://github.com/markowetzlab/scAbsolute.git
```
```bash
git clone https://github.com/markowetzlab/scUnique.git
```

We recommend using a cluster environment for the first part of the analysis, in particular if the data set comprises in the order of 100s of cells. 
_scAbsolute_ is easy to parallelize across cells, and the speedup is linear in the number of CPUs. 

## Getting started

```bash
git clone https://github.com/markowetzlab/scDNAseq-workflow.git && cd scDNAseq-workflow
```

A realistic, but small dataset of one hundred bam files each (PEO1/PEO4) can be downloaded from [here](https://drive.google.com/drive/folders/1402zegR4H7tWFluc2el9lyUr9H8rMXX6?usp=sharing). Move the data to `./data/aligned/PEO1` and `./data/aligned/PEO4`. The workflow can produce rCNAs and copy number profiles for both data sets as described in the _scUnique_ manuscript.

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
    └── vignettes                   # documentation and tutorials


## Usage

To configure this workflow, modify config/ according to your needs.
* Add per-sample folder with bam files to data/aligned folder (as with the PEO1/PEO4 folders).
* Create per sample configuration files to the config folder (one file per sample, see PEO1/PEO4.tsv examples).
* Edit variables in config/config.yaml as appropriate.

### Per-cell analysis

Given that the workflow has been properly configured, the **scAbsolute** part can be executed as follows:
```bash
export SINGULARITY_DOCKER_USERNAME=AWS SINGULARITY_DOCKER_PASSWORD=$(aws ecr-public get-login-password --region us-east-1)
```
```bash
snakemake --cores 32 --snakefile workflow/Snakefile_absolute --software-deployment-method conda apptainer --use-conda --use-singularity --singularity-args "-B /home/${USER}/.cache -B /home/${USER}/scAbsolute:/opt/scAbsolute"
```

### QC analysis

Please take the time to analyze the data (the qc script to be used for this step is available at scripts/qc-script.R)
```bash
conda env create -f envs/copynumber.yml
```
```bash
conda activate copynumber
```
Run the qc script available in workflow/scripts/qc-script.R, ensuring the cutoffs are reasonable given the plots generated. Some information about the process can be found in the 
vignette (vignettes/vignette-rCNA). Depending on the outcome (e.g. the ploidy estimates), it might be worth re-running the initial pipeline with more stringent ploidy constraints or an updated bin size.
All files passing the qc step will be added to a sample-specific TSV file in results/pass_qc. Only cells listed in any one of the files in the pass_qc folder will be in the subsequent joint analysis step.

### Joint copy number analysis and recent copy number aberration (rCNA) detection

Note that you should not move the output of the first step before running the second step of the pipeline, as this will result in errors in the output file detection. The second part of the pipeline (**scUnique**) can be run as follows:
```bash
export SINGULARITY_DOCKER_USERNAME=AWS SINGULARITY_DOCKER_PASSWORD=$(aws ecr-public get-login-password --region us-east-1)
```
```bash
snakemake --cores 16 --snakefile workflow/Snakefile_unique --software-deployment-method conda apptainer --use-conda --use-singularity --singularity-args "-B /home/${USER}/.cache -B /home/${USER}/scUnique:/opt/scUnique -B /home/${USER}/scAbsolute:/opt/scAbsolute
```

Results are then available in results/sample_name. See vignette-rCNAs for examples of how to analyze the results in practice.


## FAQ

**Q: Can I use a different directory structure for the project?**

**A:** In this case, the BASEDIR variable in the qc-script has to be adapted, and the path to the singularity bind dirs (-B) has to be adapted on the host side.


**Q: What bin size should I choose?**

**A:** I would recommend choosing a bin size such that the effective reads per bin and copy number (the rpc variable in the code/rho variable in the manuscript) are larger than 50 (and ideally around 100). The bin sizes in this package are limited by the available bin sizes in QDNAseq (1, 5, 10, 15, 30, 50, 100, 500, and 1000 kbp). We have added support for 200kb in our implementation, and it's possible to add other custom bin sizes as per the QDNAseq instructions. 


**Q: My data is very shallow. Can I use *absolute* and the rest of the pipeline?**

**A:** It's not possible to run the workflow on very sparse datasets, as *scAbsolute* requires a sufficient number of bins to fit the Gaussians. The maximum bin size supported by the entire workflow is 1MB. The minimum bin size mainly depends on and is limited by the sequencing depth. From experience, any cells with substantially less than 300,000 reads are unsuitable for this approach, and should be excluded from the analysis.


**Q: My data is single-cell RNAseq not DNAseq. Can I use this workflow?**

**A:** No, this method has been specifically developed for single-cell DNA sequencing data.


## Other question?

Please open an issue or contact the author directly via email.
