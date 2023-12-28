rule combine:
    input:
        rds=expand("results/scale/"+ str(config["binSize"]) + "/predict/{sample}.rds", sample=SAMPLE_FILES)
    output:
        merged_rds="results/scale/"+ str(config["binSize"]) + "/out.rds"
    singularity:
        IMAGE
    message:
        "Combining per cell objects into a per library object"
    threads:
        8
    shell:
        """
        set +eu
        . /opt/conda/etc/profile.d/conda.sh
        conda activate conda_runtime
        which conda
        set -eu
        export OMP_NUM_THREADS=4
        Rscript --vanilla "workflow/scripts/merge.R" "{output.merged_rds}" {input.rds}
        """
