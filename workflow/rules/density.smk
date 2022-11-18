rule density:
    input:
        bam="data/align/{sample}.bam"
    output:
        position="data/align/{sample}.position.tsv",
    singularity:
        IMAGE
    conda:
        "envs/position_search.yaml"
    message:
        "Create position table of reads for [{input}]"
    threads:
        1
    shell:
        """
        type conda
        /bin/bash /opt/scAbsolute/data/readPosition/extract-start-sites.sh {input.bam} /opt/scAbsolute/data/readPosition/assembly.tsv;
        """
||||||| parent of dfdbefc (working version to reproduce scAbsolute)
=======
rule density:
    input:
        bam="data/align/{sample}.bam"
    output:
        position="data/align/{sample}.position.tsv",
    singularity:
        IMAGE
    conda:
        "envs/position_search.yaml"
    message:
        "Create position table of reads for [{input}]"
    threads:
        1
    shell:
        """
        which conda
        /bin/bash /opt/scAbsolute/data/readPosition/extract-start-sites.sh {input.bam} /opt/scAbsolute/data/readPosition/assembly.tsv;
        """
