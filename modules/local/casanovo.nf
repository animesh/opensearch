process CASANOVO {
    errorStrategy 'ignore'
    tag "${sample_id}"
    label 'process_high'

    publishDir "${params.outdir}/${params.casanovo_bin.toString().tokenize('/')[-1]}", mode: 'copy'

    input:
    tuple val(sample_id), val(mzml_path)

    output:
    tuple val(sample_id), path("*.${params.casanovo_bin.toString().tokenize('/')[-1]}"), emit: casanovo_dir
    tuple val(sample_id), path("${sample_id}_casanovo_mqc.tsv"), emit: mqc
    path "versions.yml", emit: versions

    script:
    def tool_name = params.casanovo_bin.toString().tokenize('/')[-1]
    """
    set -euo pipefail

    OUTDIR="${sample_id}.${params.casanovo_bin.toString().tokenize('/')[-1]}"
    mkdir -p "\$OUTDIR"
    MZML="${mzml_path}"
    OUTPUT_ROOT="${sample_id}"
    LOG="\$OUTDIR/\${OUTPUT_ROOT}.log"
    MZTAB="\$OUTDIR/\${OUTPUT_ROOT}.mztab"

    rc=0
    status="SUCCESS"
    message="${params.casanovo_bin.toString().tokenize('/')[-1]} completed successfully using its installed default model"

    [[ -s "\$MZML" ]] || { rc=2; status="FAILED"; message="Casanovo input mzML does not exist: \$MZML"; }

    if [[ \$rc -eq 0 ]]; then
        set +e
        "${params.casanovo_bin}" sequence "\$MZML" --output_dir "\$OUTDIR" --output_root "\$OUTPUT_ROOT"
        rc=\$?
        set -e
        if [[ \$rc -ne 0 ]]; then
            status="FAILED"
            message="Casanovo exited with code \$rc. See \$LOG."
        fi
    fi

    sequenced=0; score50=0; score90=0
    if [[ -s "\$LOG" ]]; then
        sequenced=\$(grep -oP 'Sequenced \\K[0-9]+' "\$LOG" | tail -1 || echo 0)
        score50=\$(grep -oP '([0-9]+) spectra \\(.*?\\) scored ≥ 0\\.50' "\$LOG" | grep -oP '^[0-9]+' | tail -1 || echo 0)
        score90=\$(grep -oP '([0-9]+) spectra \\(.*?\\) scored ≥ 0\\.90' "\$LOG" | grep -oP '^[0-9]+' | tail -1 || echo 0)
    fi
    if [[ \$status == SUCCESS && ! -s "\$MZTAB" ]]; then
        status="FAILED"
        message="Casanovo returned success but produced no mzTab result."
    fi

    printf 'sample\\ttool\\tstatus\\texit_code\\tmessage\\n' > "\$OUTDIR/status.tsv"
    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "${sample_id}" "${params.casanovo_bin.toString().tokenize('/')[-1]}" "\$status" "\$rc" "\$message" >> "\$OUTDIR/status.tsv"

    printf '# id: casanovo_stats\\n# plot_type: generalstats\\n# pconfig:\\n#   sequenced_spectra:\\n#     title: Sequenced\\n#   score_ge_50pct:\\n#     title: "Score>=0.5"\\n#   score_ge_90pct:\\n#     title: "Score>=0.9"\\nSample\\tsequenced_spectra\\tscore_ge_50pct\\tscore_ge_90pct\\n%s\\t%s\\t%s\\t%s\\n' \\
        "${sample_id}" "\$sequenced" "\$score50" "\$score90" > "${sample_id}_casanovo_mqc.tsv"

    printf 'tool\\tsample\\trole\\tpath\\n' > "\$OUTDIR/opensearch_manifest.tsv"
    printf '%s\\t%s\\tstatus\\tstatus.tsv\\n' "${params.casanovo_bin.toString().tokenize('/')[-1]}" "${sample_id}" >> "\$OUTDIR/opensearch_manifest.tsv"
    [[ -s "\$LOG" ]] && printf '%s\\t%s\\tlog\\t%s\\n' "${params.casanovo_bin.toString().tokenize('/')[-1]}" "${sample_id}" "\${LOG#\$OUTDIR/}" >> "\$OUTDIR/opensearch_manifest.tsv"
    [[ -s "\$MZTAB" ]] && printf '%s\\t%s\\tmztab\\t%s\\n' "${params.casanovo_bin.toString().tokenize('/')[-1]}" "${sample_id}" "\${MZTAB#\$OUTDIR/}" >> "\$OUTDIR/opensearch_manifest.tsv"

    CASANOVO_VERSION=\$("${params.casanovo_bin}" --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || true)
    CASANOVO_VERSION=\${CASANOVO_VERSION:-unknown}
    cat > versions.yml <<END_VERSIONS
    "${task.process}":
        ${params.casanovo_bin.toString().tokenize('/')[-1]}: \$CASANOVO_VERSION
    END_VERSIONS

    exit 0
    """
}
