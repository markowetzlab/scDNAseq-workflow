"""
Author: M. P Schneider
Affiliation: Cancer Research UK, CI
Program: Workflow for single cell or sWGS DNAseq absolute copy number calling - unknown ploidy
Date: April 2020
"""
import numpy as np
import pandas as pd
import itertools
import sys
import glob, os
from os.path import expandvars
from pprint import pprint
from snakemake.utils import validate
workdir: expandvars("$HOME")

ROOTDIR = expandvars("$HOME/Data" + "/" + config["PROJECT"] + "/")
SAMPLE_DESCRIPTIONS = expandvars("$HOME/Data/" + config["PROJECT"] + "/metadata/")
# NOTE the bind argument in snakemake is evaluated after the script is run -> i.e. paths outside of rules refer to host, inside refer to container
INPUT_DIR = "20-align"
print(config["binSize"])
OUTPUT_DIR = "30-scale/" + str(config["binSize"]) + "/" + "predict"
IMAGE = "docker://770395835384.dkr.ecr.eu-central-1.amazonaws.com/scabsolute:v2.5.0"
assert("UID" in config)
assert(isinstance(config["binSize"], int))
config["WALLTIME"] = int(config["WALLTIME"])

# allow format of variables as list or as single string instance in snakeconfig file
def input_format(variables):
    for var in variables:
        if var in config and (isinstance(config[var], str) or isinstance(config[var], int)):
            config[var] = [config[var]]

input_format(["UID", "LIBRARY", "ALIGNMENT"])

## Query and select all samples
sample_dfs = []
for file in glob.glob(expandvars(SAMPLE_DESCRIPTIONS)+"copynumber/**/*.tsv", recursive=True):
    # print(file)
    samples = pd.read_table(file, sep="\t",
        names=["sample", "library", "filename", "path", "readcounts", "readcounts.dedup"])
    validate(samples, "samples.schema.yaml")
    sample_dfs.append(samples)

SAMPLES = pd.concat(sample_dfs)
EXCLUSION = pd.read_table(expandvars(SAMPLE_DESCRIPTIONS)+"copynumber/exclusion.list", sep="\t", names=["sample",
    "library", "filename", "path", "readcounts", "readcounts.dedup"])["filename"]
FAILED = pd.read_table(expandvars(SAMPLE_DESCRIPTIONS)+"copynumber/failed.list", sep="\t", names=["sample", "library", "filename", "path", "readcounts", "readcounts.dedup"])["filename"]
EXCLUSION = pd.concat([EXCLUSION, FAILED])

SAMPLES = SAMPLES.loc[SAMPLES['sample'].isin(config["UID"]) & \
    (~((SAMPLES['filename'].isin(EXCLUSION)) | \
      (SAMPLES['readcounts.dedup'] < config["lowread"])))]
if "LIBRARY" in config:
    SAMPLES = SAMPLES.loc[SAMPLES['library'].isin(config["LIBRARY"])]
SAMPLE_FILES = glob_wildcards(ROOTDIR+INPUT_DIR+"/{samples}.bam")[0]
OUTPUT = SAMPLES["sample"].unique() + "_" + str(config["binSize"])
INPUT = [os.path.splitext(i)[0] for i in SAMPLES['path'].tolist()]

print(list(set(INPUT) - set(SAMPLE_FILES)))
#print([elem in SAMPLE_FILES for elem in INPUT])
assert(all(elem in SAMPLE_FILES for elem in INPUT))

ruleorder: qc > position_search > scale_scAbsolute > combine

rule all:
    input:
        all=expand(ROOTDIR+OUTPUT_DIR+"/{out}.rds",zip,out=OUTPUT)
# pprint(rules.all.input)
# sys.exit(2)

rule qc:
    input:
        bam=ROOTDIR+INPUT_DIR+"/{sample}.bam"#, sample=INPUT)
    output:
        flagstat=ROOTDIR+INPUT_DIR+"/{sample}.flagstat",
        bai=ROOTDIR+INPUT_DIR+"/{sample}.bam.bai"
    params:
        project=config["PROJECT"],
    singularity:
        "docker://770395835384.dkr.ecr.eu-central-1.amazonaws.com/biocontainers/samtools:v1.9-4-deb_cv1"
    message:
        "Quality control on bam file [{input}]"
    threads:
        1
    resources:
        mem_mb=(lambda wildcards, attempt: attempt * 3420),
        walltime=config["WALLTIME"]
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        sorted=$(samtools stats {input.bam} | grep "is sorted:" | cut -f 3)
        if [ $sorted -eq "1" ]; then
            samtools index {input.bam}
        else
            exit 1;
        fi
        """

rule scale_scAbsolute:
    input:
        bam=ROOTDIR+INPUT_DIR+"/{sample}.bam",
        flagstat=ROOTDIR+INPUT_DIR+"/{sample}.flagstat",
        bai=ROOTDIR+INPUT_DIR+"/{sample}.bam.bai",
        position=ROOTDIR+INPUT_DIR+"/{sample}.position.tsv"
        # bam=rules.qc.input.bam,
        # flagstat=rules.qc.output.flagstat
    params:
        project=config["PROJECT"],
        prefix=lambda wildcards, input: os.path.dirname(input.bam),
        filefix=lambda wildcards, input: os.path.basename(input.bam)
    output:
        rds=temp(ROOTDIR+OUTPUT_DIR+"/{sample}.rds")
    singularity:
        IMAGE
    message:
        "Calling absolute copy number profile for {input}"
    threads:
        1
    resources:
        mem_mb=(lambda wildcards, attempt: attempt * 3420),
        walltime=config["WALLTIME"]
    shell:
        """
        set +eu
        . /opt/conda/etc/profile.d/conda.sh
        conda activate conda_runtime
        which conda
        set -eu
        # note, not compatible with conda_runtime
        # export PYTHONPATH=/opt/conda/bin/python
        # export SINGULARITYENV_PYTHONPATH=/opt/conda/bin/python
        export RETICULATE_PYTHON=/opt/conda/envs/conda_runtime/bin/python
        export MKL_THREADING_LAYER=sequential
        export OMP_NUM_THREADS=2
        which python
        python -c "import tensorflow; import numpy; import pandas;"
        Rscript -e "library(reticulate); reticulate::py_discover_config();"
        Rscript --vanilla "scAbsolute/scripts/10-scale-scdata.R" "{params.prefix}" "{params.filefix}" "{output.rds}" "{config[binSize]}"
        """

def tmp_write(files):
# BUGFIX https://bitbucket.org/snakemake/snakemake/issues/878/errno-7-argument-list-too-long-path-to-bin
    import tempfile
    with tempfile.NamedTemporaryFile(delete=False, dir=os.path.expandvars("$HOME/.cache/"), mode="w") as fn:
        for item in files:
            fn.write("%s\n" % item)

    return fn.name

def collect(wildcards):
    # collect bam files based on mappings, i.e. apply input/output mapping
    output_sample = wildcards.out.split("_")[0]
    samples = SAMPLES[SAMPLES["sample"] == output_sample]
    samples = samples["path"].tolist()
    samples = [ROOTDIR + OUTPUT_DIR + "/" + os.path.splitext(i)[0]+'.rds' for i in samples]
    # remove duplicates
    samples = list(set(samples))
    return samples

rule combine:
    input:
        rds=collect
    params:
        project=config["PROJECT"],
        # NOTE, can be removed once snakemake bug is fixed
        rds_file=lambda wildcards, input: tmp_write(collect(wildcards)),
    output:
        merged_rds=ROOTDIR+OUTPUT_DIR+"/{out}.rds"#,out=OUTPUT)
    singularity:
        IMAGE
    message:
        "Combining per cell objects into a per library object"
    threads:
        8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 64000,
        walltime=config["WALLTIME"]
    shell:
        """
        set +eu
        . /opt/conda/etc/profile.d/conda.sh
        conda activate conda_runtime
        which conda
        set -eu
        export OMP_NUM_THREADS=4
        cat {params.rds_file}
        Rscript --vanilla scAbsolute/scripts/merge.R "{output.merged_rds}" $(cat {params.rds_file})
        rm -f {params.rds_file}
        """
