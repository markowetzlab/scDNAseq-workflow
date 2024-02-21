rule scale_scAbsolute:
    input:
        bam="data/aligned/" + str(config["sampleName"]) + "/{sample}.bam",
        flagstat="data/aligned/" + str(config["sampleName"]) + "/{sample}.flagstat",
        bai="data/aligned/" + str(config["sampleName"]) + "/{sample}.bam.bai",
        position="data/aligned/" + str(config["sampleName"]) + "/{sample}.position.tsv" if config["estimateReadDensity"] else [],
    params:
        prefix=lambda wildcards, input: os.path.dirname(input.bam),
        filefix=lambda wildcards, input: os.path.basename(input.bam),
        minPloidy=lambda wildcards, input: LOOKUP_PLOIDY[wildcards.sample]["minPloidy"] if wildcards.sample in LOOKUP_PLOIDY else "NNULL",
        maxPloidy=lambda wildcards, input: LOOKUP_PLOIDY[wildcards.sample]["maxPloidy"] if wildcards.sample in LOOKUP_PLOIDY else "NNULL",
    output:
        rds="results/" + str(config["binSize"]) + "/" + str(config["sampleName"]) + "/" + "{sample}.rds"
    container:
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
        Rscript -e "library(reticulate); reticulate::py_discover_config();"
        if [ "{params.minPloidy}" -eq "NULL" ]; then
            Rscript --vanilla "workflow/scripts/run_scAbsolute.R" "{config[species]}" \
                                                                  "{config[genome]}" \
                                                                  "{params.prefix}" \
                                                                  "{params.filefix}" \
                                                                  "{output.rds}" \
                                                                  "{config[binSize]}" \
                                                                  "{config[estimateReadDensity]}" || true
        else
            Rscript --vanilla "workflow/scripts/run_scAbsolute.R" "{config[species]}" \
                                                                  "{config[genome]}" \
                                                                  "{params.prefix}" \
                                                                  "{params.filefix}" \
                                                                  "{output.rds}" \
                                                                  "{config[binSize]}" \
                                                                  "{config[estimateReadDensity]}" \
                                                                  "{params.minPloidy}" \
                                                                  "{params.maxPloidy}" || true
        fi
        """
