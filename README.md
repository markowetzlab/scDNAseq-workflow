# Snakemake workflow for single-cell DNA sequencing analysis (_scAbsolute_, _scUnique_)

## Overview
This workflow performs absolute copy number calling and detection of recent copy number aberrations (rCNAs) in single-cell DNA sequencing data as described in [*scAbsolute*](https://www.biorxiv.org/content/10.1101/2022.11.14.516440v2) and [*scUnique*](missing).
Please cite the relevant publications if you use this pipeline:
- *scAbsolute*:
```
Schneider, Michael P., et al. 2024.
scAbsolute: Measuring Single-Cell Ploidy and Replication Status.
Genome Biology 25 (1): 62.
https://doi.org/10.1186/s13059-024-03204-y.
```
- *scUnique*:
```
we will release the code to use scUnique upon publication. please contact the author if you need earlier access
```
- if you use *scUnique*, you should also cite [MEDICC2](https://doi.org/10.1186/s13059-022-02794-9):
```
Kaufmann, T.L., Petkovic, M., Watkins, T.B.K. et al.  
MEDICC2: whole-genome doubling aware copy-number phylogenies for cancer evolution.  
Genome Biol 23, 241 (2022).  
https://doi.org/10.1186/s13059-022-02794-9
```

Here is how to use the approach:
1. Initially, run the per-cell part of the pipeline (*scAbsolute* + single-cell copy number calling)
2. Inspect results and run quality control and outlier detection across all cells in each sample. If necessary, rerun the ploidy analysis (step 1) with an updated ploidy window. Per-cell ploidy windows can be specified in the sample files (see config/sample_PEO1.tsv for reference). I highly recommend manually inspecting results at this stage, as per-sample variations are considerable in our experience and cutoffs may need to be adjusted on a per-sample basis. 
3. Run the second stage of the pipeline (_scUnique_), combining the results from the per-cell analysis and creating a joint analysis of the dataset, including an analysis of recent copy number aberrations (rCNAs).

## Dependencies
You must install [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (we support version 8 and above) and [singularity/apptainer](https://apptainer.org/docs/admin/main/installation.html). To access the containers hosted on AWS ECR public registry, you can use the [AWS CLI application](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).
You must clone the *scAbsolute* and *scUnique* GitHub repositories to your local machine. We demonstrate usage with cloning to the home directory, but it is also possible to use another directory, but you will have to modify the BASEDIR variable.
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

A realistic, but small dataset of one hundred (sorted) bam files each (PEO1/PEO4) can be downloaded from [here](https://drive.google.com/drive/folders/1402zegR4H7tWFluc2el9lyUr9H8rMXX6?usp=sharing). Move the data to `./data/aligned/PEO1` and `./data/aligned/PEO4`. The workflow can produce rCNAs and copy number profiles for both data sets as described in the _scUnique_ manuscript.

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
        ├── scripts                 # main analysis scripts / edit for additional customization
        └── Snakefile               
    └── vignettes                   # documentation and tutorials


## Usage

To configure this workflow, modify config/ according to your needs.
* Add per-sample folder with **coordinate-sorted and duplicate-marked** bam files (one bam file per cell) to data/aligned folder (as with the PEO1/PEO4 folders).
* Create per sample configuration files to the config folder (one file per sample, see PEO1/PEO4.tsv examples).
* Edit variables in config/config.yaml as appropriate.

### Per-cell analysis

Given that the workflow has been properly configured, the **scAbsolute** part can be executed as follows:
```bash
snakemake --cores 16 --snakefile workflow/Snakefile_absolute \
    --software-deployment-method conda apptainer --use-conda --use-singularity \
    --singularity-args "-B /home/${USER}/.cache -B /home/${USER}/scAbsolute:/opt/scAbsolute"
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
vignette (vignettes/vignette-rCNA). Depending on the outcome (the ploidy and rho/rpc estimates, respectively), it might be worth re-running the initial pipeline with more stringent ploidy constraints or an updated bin size.
All files passing the qc step will be added to a sample-specific TSV file in results/pass_qc. Only cells listed in any one of the files in the pass_qc folder will be included in the subsequent joint analysis step.
For the tutorial run with the PEO1_subset sample, you can use the example file `config/pass_all_PEO1_subset.tsv` and move it to `results/pass_qc/pass_all.tsv`, to execute the next step of the analysis on all samples in the PEO1_subset dataset. These are all relatively high-quality cells, so it is not really necessary to remove any of them.
While you will ideally run the `qc-script.R` interactively, it is possible to run it from the command line:
```bash
Rscript workflow/scripts/qc-script.R "~/scDNAseq-workflow/results/500/PEO1_subset_500.rds" "PEO1_subset"
```


### Joint copy number analysis and recent copy number aberration (rCNA) detection

Note that you should not move the output of the first step before running the second step of the pipeline, as this will result in errors in the output file detection. The second part of the pipeline (**scUnique**) can be run as follows:
```bash
snakemake --cores 16 --snakefile workflow/Snakefile_unique \
    --software-deployment-method conda apptainer --use-conda --use-singularity \
    --singularity-args "-B /home/${USER}/.cache -B /home/${USER}/scUnique:/opt/scUnique -B /home/${USER}/scAbsolute:/opt/scAbsolute"
```

Results are then available in results/sample_name. 

### More examples and explanations

Please see the vignettes folder for more examples and in-depth explanations, particularly the vignette-rCNAs file for examples of how to analyze the results and some examples of plotting the data.



## FAQ

**Q: I have issues with the container authentication. What can I do?**

**A:** Try running the command below and ensure you are not pulling from the wrong container registry. First, we recently moved to the public docker hub registry from AWS, so make sure you are pulling from the correct registry.


**Q: Why are all your bam files named the way they are?**

**A:**  We have a custom naming convention for our data. We recommend following this convention with your own data for compatibility reasons.
We include the following information in every sample name:
- Unique sample identifier
- Library identifier
- Cell number (padded to 6 digits)
- Cell tag, i.e., all cell specific information, such as plate and well numbers, custom sequences, and or experimental information
The different fields are separated by an underscore, and the file name then looks like these examples:
```
  [unique sample identifier]_[library identifier]_[cell id]_[cell tag].bam
  UID-10X-Fibroblast-cell_SLX-00000_000001_AAACCTGAGCAGCCTC-1.bam
  UID-JBL-CIOV1_SLX-17170_000001_AAAGCAAGTAGAACAT-1
  UID-FML-PEO1-2N_SLX-23003_000288_SINCEL-194-Plate-236-P21-CGACAAGGAT-AAGGCCACCT.bam
```


**Q: What is the input for the programme?**

**A:** Our workflow requires **coordinate-sorted and duplicate-marked** bam files (one file per cell) with a minimum read depth of about 200,000 reads per cell. The actual required read depth depends on the sample ploidy and bin size. Knowing these parameters you can then easily compute the required number of reads.


**Q: How do I compute the number of reads I should aim for in my experiment?**

**A:** Three elements determine the effective number of reads you want to obtain per cell. Bin size (we support 1, 5, 10, 15, 30, 50, 100, 200, 500, and 1000 kb bins), sample ploidy (defined as the average copy number state across the genome), and the target number of reads per copy and per bin (we refer to this in the publication as ρ and in the code as rpc).
Let's say you want to check if your sample harbours rCNAs at a bin size of 500kb, and you know your sample has a ploidy of 2.5. Assume the quality of your sequencing is similar to the DLP sequencing assay, and you would be okay with an FDR of about 10%. In this case, you can target a value of ρ = 50. This gives us the required number of reads as n_reads = ρ * ploidy * n_bins, where n_bins is about 5000 (Note that the genome at 500kb resolution has closer to 6000 bins, but we lose about 1000 bins due to quality issues, and we recommend ignoring the sex chromosomes for rCNA analysis). 
This gives us a target per-cell number of 50 * 2.5 * 5000 ~ 600,000 reads. Note that it makes sense to target a slightly higher number (at least an additional 10-20%, due to duplicated reads, the fact that we blacklist a certain proportion of bins, and the fact that we obtain a distribution of reads over the cells and we should aim to keep most of the cells at a ρ value of 50 and over). 
Similarly, if we have a sample with a ploidy of 1.5, and we want to do an analysis at 100kb resolution at ρ=50, we then end up with close to 2 million reads per cell (n_reads = 50 * 1.5 * 25000 = 1,875,000).
Finally, another scenario. Assume we have sequenced a precious tumour sample with an unknown ploidy at an average sequencing depth of 3 million reads per cell. Running the first part of the workflow we determine the average cell ploidy to be 2.7. Again, given the sequencing quality, we need at least ρ=50 to have some trust in the per-cell copy number calls. We then can ask what bin size we can analyze the data at. n_bins = n_reads / (ρ * ploidy) = 3e6 / (50 * 2.7) ~ 22 000. We could, therefore, probably work with 100kb or more cautiously with 200kb bin size. 


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
