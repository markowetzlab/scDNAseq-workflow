#!/bin/bash
#!
#! SLURM job script for clust1
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#SBATCH --job-name=SLX-25393
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --time=168:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shadi.shafighi@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH --output=logs/log.%A.log
#SBATCH --error=errs/%x_slurm.%A.%a.err


###############################################################
### Execute script                                         ####
###############################################################

#!echo "cores: $1"
#!echo "output: $2"
dataset=SLX-25393

#echo -e "sample_name" > config/$dataset.tsv
#find data/aligned/$dataset -type f -name '*.bam' -exec basename {} \; | sed 's/\.bam$//' >> config/$dataset.tsv
#snakemake -s workflow/Snakefile_absolute --cores 400 --nolock --keep-going --rerun-incomplete --use-singularity  --singularity-args "--bind /mnt" --debug-dag --use-conda results/500/${dataset}_500.rds
snakemake --cores 32 --resources mem_mb=150000  --nolock --snakefile workflow/Snakefile_unique --software-deployment-method conda apptainer --use-conda --use-singularity --singularity-args "-B /home/darvis01/.cache -B /mnt/scratchc/fmlab/darvis01/scDNAseq-workflow/scAbsolute:/opt/scAbsolute -B /mnt/scratchc/fmlab/darvis01/scDNAseq-workflow/scUnique-edited:/opt/scUnique"

