
if config["estimateReadDensity"]:
    rule density:
        input:
            bam="data/aligned/{sample}.bam"
        output:
            position="data/aligned/{sample}.position.tsv",
        container:
            IMAGE
        conda:
            "envs/position_search.yaml"
        message:
            "Create position table of reads for [{input}]"
        threads:
            1
        shell:
            """
            if [[ "{config[genome]}" == "hg19" ]] || [[ "{config[genome]}" == "GRCh37" ]]; then
                /bin/bash /opt/scAbsolute/data/readPosition/extract-start-sites.sh {input.bam} /opt/scAbsolute/data/readPosition/assembly.tsv
            elif [[ "{config[genome]}" == "mm10" ]] || [[ "{config[genome]}" == "GRCm38" ]]; then
                /bin/bash /opt/scAbsolute/data/readPosition/extract-start-sites.sh {input.bam} /opt/scAbsolute/data/readPosition/assembly_mouse.tsv
            else
                echo "Pick a valid genome."
                exit 1;
            fi
            """

