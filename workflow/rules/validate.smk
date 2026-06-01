rule validate_bams:
    """Pre-flight BAM integrity check. Runs samtools quickcheck on every BAM
    listed in the sample sheet and fails fast with a full list of bad files
    before any qc / scale_scAbsolute job is scheduled. Common failure mode
    this catches: BAMs that were truncated during copy from the source path.

    On failure, writes the bad cells into the SAME failed_cells.csv that the
    combine / merge step uses downstream, with failure_reason = "truncated_bam".
    This gives a single canonical place to look for any failed cell across the
    whole pipeline (truncated input, missing output, process crash, or QC
    failure during scAbsolute).
    """
    input:
        bams=expand("data/aligned/" + str(config["sampleName"]) + "/{sample}.bam",
                    sample=SAMPLE_FILES)
    output:
        report="data/aligned/" + str(config["sampleName"]) + "/.bam_validation.passed"
    params:
        failed_csv="results/" + str(config["binSize"]) + "/" + str(config["sampleName"]) + "_" + str(config["binSize"]) + "_failed_cells.csv"
    container:
        config["qc_img"]
    message:
        "Pre-flight BAM integrity check (samtools quickcheck)"
    threads:
        1
    shell:
        r"""
        set -uo pipefail
        bad_tmp=$(mktemp)
        n_total=0
        n_bad=0
        for bam in {input.bams}; do
            n_total=$((n_total+1))
            if ! samtools quickcheck "$bam" 2>/dev/null; then
                # Strip .bam extension to match cell-name convention used downstream
                base=$(basename "$bam" .bam)
                echo "$base" >> "$bad_tmp"
                n_bad=$((n_bad+1))
            fi
        done

        if [ "$n_bad" -gt 0 ]; then
            # Write failed cells into the canonical failed_cells.csv (same file
            # merge.R uses) so the user has ONE place to look for failures.
            mkdir -p "$(dirname {params.failed_csv})"
            {{
                echo '"name","failure_reason"'
                while IFS= read -r cell; do
                    printf '"%s","truncated_bam"\n' "$cell"
                done < "$bad_tmp"
            }} > {params.failed_csv}

            echo "" >&2
            echo "============================================================" >&2
            echo "ERROR: $n_bad / $n_total BAM file(s) failed integrity check:" >&2
            echo "============================================================" >&2
            cat "$bad_tmp" >&2
            echo "============================================================" >&2
            echo "Full failure report written to:" >&2
            echo "  {params.failed_csv}" >&2
            echo "(same CSV that downstream failures use; failure_reason = truncated_bam)" >&2
            echo "Common cause: incomplete copy from source filesystem." >&2
            echo "Action: re-copy listed files from source and re-run." >&2
            echo "============================================================" >&2
            rm -f "$bad_tmp"
            exit 1
        fi
        rm -f "$bad_tmp"
        echo "All $n_total BAMs passed samtools quickcheck on $(date)" > {output.report}
        """
