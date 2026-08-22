process CASANOVO {
    errorStrategy 'ignore'
    tag "${sample_id}"
    label 'process_high'

    input:
    tuple val(sample_id), path(mzml_file)

    output:
    tuple val(sample_id), path("${sample_id}.${params.casanovo_workdir_suffix}"), emit: casanovo_dir
    tuple val(sample_id), path("${sample_id}_casanovo_mqc.tsv"), emit: mqc
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    OUTDIR="${sample_id}.${params.casanovo_workdir_suffix}"
    mkdir -p "\$OUTDIR"

    MZML="\$(realpath "${mzml_file}")"
    rc=0
    status="SUCCESS"
    message="Casanovo completed successfully using its installed default model"

    echo "[OpenSearch] Casanovo input: \$MZML"
    # Deliberately do NOT pass --model. Casanovo's installed default model is
    # used for both Orbitrap and timsTOF data.
    set +e
    "${params.casanovo_bin}" sequence "\$MZML" --output_dir "\$OUTDIR"
    rc=\$?
    set -e

    if [[ \$rc -ne 0 ]]; then
        status="FAILED"
        message="Casanovo exited with code \$rc. See Casanovo log / Nextflow .command.log."
        echo "WARNING: ${sample_id}: \$message" >&2
    fi

    log_file="\$(find "\$OUTDIR" -maxdepth 1 -type f -name 'casanovo_*.log' -print -quit)"
    mztab_file="\$(find "\$OUTDIR" -maxdepth 1 -type f -name 'casanovo_*.mztab' -print -quit)"
    sequenced=0; score50=0; score90=0

    if [[ -n "\$log_file" ]]; then
        sequenced=\$(grep -oP 'Sequenced \\K[0-9]+' "\$log_file" | tail -1 || echo 0)
        score50=\$(grep -oP '([0-9]+) spectra \\(.*?\\) scored ≥ 0\\.50' "\$log_file" | grep -oP '^[0-9]+' | tail -1 || echo 0)
        score90=\$(grep -oP '([0-9]+) spectra \\(.*?\\) scored ≥ 0\\.90' "\$log_file" | grep -oP '^[0-9]+' | tail -1 || echo 0)
    fi

    if [[ \$status == SUCCESS && -z "\$mztab_file" ]]; then
        status="FAILED"
        message="Casanovo returned success but produced no mzTab result."
    fi

    printf 'sample\\ttool\\tstatus\\texit_code\\tmessage\\n' > "\$OUTDIR/status.tsv"
    printf '%s\\tCasanovo\\t%s\\t%s\\t%s\\n' "${sample_id}" "\$status" "\$rc" "\$message" >> "\$OUTDIR/status.tsv"

    printf '# id: casanovo_stats\\n# plot_type: generalstats\\n# pconfig:\\n#   sequenced_spectra:\\n#     title: Sequenced\\n#   score_ge_50pct:\\n#     title: "Score>=0.5"\\n#   score_ge_90pct:\\n#     title: "Score>=0.9"\\nSample\\tsequenced_spectra\\tscore_ge_50pct\\tscore_ge_90pct\\n%s\\t%s\\t%s\\t%s\\n' \\
        "${sample_id}" "\$sequenced" "\$score50" "\$score90" > "${sample_id}_casanovo_mqc.tsv"

    CASANOVO_VERSION=\$("${params.casanovo_bin}" --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || true)
    CASANOVO_VERSION=\${CASANOVO_VERSION:-unknown}
    cat > versions.yml <<END_VERSIONS
    "${task.process}":
        casanovo: \$CASANOVO_VERSION
    END_VERSIONS

    exit 0
    """
}
