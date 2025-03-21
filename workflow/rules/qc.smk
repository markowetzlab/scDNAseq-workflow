rule qc:
    input:
        bam="data/aligned/" + str(config["sampleName"]) + "/{sample}.bam"
    output:
        flagstat="data/aligned/" + str(config["sampleName"]) + "/{sample}.flagstat",
        bai="data/aligned/" + str(config["sampleName"]) + "/{sample}.bam.bai"
    container:
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
            echo "We required sorted bam files as input"
            exit 1;
        fi
        """
