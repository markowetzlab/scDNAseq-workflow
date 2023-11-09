rule scale_scAbsolute:
    input:
        bam="data/align/{sample}.bam",
        flagstat="data/align/{sample}.flagstat",
        bai="data/align/{sample}.bam.bai",
        position="data/align/{sample}.position.tsv"
    params:
        prefix=lambda wildcards, input: os.path.dirname(input.bam),
        filefix=lambda wildcards, input: os.path.basename(input.bam)
    output:
        rds="results/scale/"+ str(config["binSize"]) + "/individual/{sample}.rds"
    singularity:
        IMAGE
    message:
        "Calling absolute copy number profile for {input}"
    threads:
        1
    shell:
        """
        set +eu
        . /opt/conda/etc/profile.d/conda.sh
        conda activate conda_runtime
        type conda
        set -eu
        export RETICULATE_PYTHON=/opt/conda/envs/conda_runtime/bin/python
        export MKL_THREADING_LAYER=sequential
        export OMP_NUM_THREADS=2
        type python
        python -c "import tensorflow; import numpy; import pandas;"
        touch "{output.rds}"
        Rscript -e "library(reticulate); reticulate::py_discover_config();"
        Rscript --vanilla "workflow/scripts/run_scAbsolute.R" "{config[species]}" "{config[genome]}" "{params.prefix}" "{params.filefix}" "{output.rds}" "{config[binSize]}" || true
        """
