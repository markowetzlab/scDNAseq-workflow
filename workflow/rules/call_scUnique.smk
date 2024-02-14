rule call_scUnique:
    input:
        rds="results/" + str(config["binSize"]) + "/" + str(config["sampleName"]) + "_" + str(config["binSize"]) + ".rds",
    params:
        longpath=lambda wildcards, input: os.path.abspath(input.rds),
        resultpath=lambda wildcards, input: os.path.abspath("./results")
    output:
        rdata="results/" + str(config["sampleName"]) + "_" + str(config["binSize"]) + "/" + str(config["sampleName"]) + "_" + str(config["binSize"]) + ".RData" ,
    container:
        IMAGE_scUnique
    message:
        "Calling recent copy number aberrations and joint copy number profiles for {input}"
    threads:
        config["NCORES"]
    shell:
        """
        set +eu
        . /opt/conda/etc/profile.d/conda.sh
        conda activate conda_runtime
        set -eu
        export RETICULATE_PYTHON=/opt/conda/envs/conda_runtime/bin/python
        python -c "import tensorflow; import numpy; import pandas;"
        export NCORES="{config[NCORES]}"
        export MKL_THREADING_LAYER=sequential
        export OMP_NUM_THREADS=1
        export PROCESSX_NOTIFY_OLD_SIGCHLD=true
        Rscript --vanilla "workflow/scripts/run_scUnique.R" {params.longpath} {params.resultpath} 0.05
        """
