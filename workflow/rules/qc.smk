rule qc:
    input:
        bam="data/align/{sample}.bam",
    output:
        flagstat="data/align/{sample}.flagstat",
        bai="data/align/{sample}.bam.bai"
    singularity:
        config["qc_img"]
    message:
        "Quality control on bam file [{input}]"
    threads:
        1
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
