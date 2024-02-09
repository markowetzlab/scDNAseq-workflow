rule combine:
    input:
        rds=expand("results/" + str(config["binSize"]) + "/" + str(config["sampleName"]) + "/" + "{sample}.rds", sample=SAMPLE_FILES)
    output:
        merged_rds="results/" + str(config["binSize"]) + "/" + str(config["sampleName"]) + "_" + str(config["binSize"]) + ".rds"
    container:
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
        type conda
        set -eu
        export OMP_NUM_THREADS=4
        Rscript --vanilla "workflow/scripts/merge.R" "{output.merged_rds}" {input.rds}
        """
