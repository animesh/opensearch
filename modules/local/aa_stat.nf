process AA_STAT {
    errorStrategy 'ignore'
    tag "${sample_id}"
    label 'process_high'

    input:
    tuple val(sample_id), path(mzml_file), path(pepxml_file)

    output:
    tuple val(sample_id), path("${sample_id}.${params.aastat_workdir_suffix}"), emit: aa_stat_dir
    tuple val(sample_id), path("${sample_id}_aastat_mqc.tsv"), emit: mqc
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    OUTDIR="${sample_id}.${params.aastat_workdir_suffix}"
    mkdir -p "\$OUTDIR"

    MZML="\$(realpath "${mzml_file}")"
    PEPXML="\$(realpath "${pepxml_file}")"
    rc=0
    status="SUCCESS"
    message="AA_stat completed successfully"

    echo "[OpenSearch] AA_stat mzML: \$MZML"
    echo "[OpenSearch] AA_stat pepXML: \$PEPXML"
    set +e
    "${params.aa_stat_bin}" \\
        -n ${task.cpus} \\
        --mzml "\$MZML" \\
        --pepxml "\$PEPXML" \\
        --dir "\$OUTDIR"
    rc=\$?
    set -e

    if [[ \$rc -ne 0 ]]; then
        status="FAILED"
        message="AA_stat exited with code \$rc. See Nextflow .command.log."
        echo "WARNING: ${sample_id}: \$message" >&2
    fi

    n_shifts=\$(grep -c '^[+-]' "\$OUTDIR/aa_statistics_table.csv" 2>/dev/null || true)
    n_loc=\$(tail -n +2 "\$OUTDIR/localization_statistics.csv" 2>/dev/null | wc -l || true)

    printf 'sample\\ttool\\tstatus\\texit_code\\tmessage\\n' > "\$OUTDIR/status.tsv"
    printf '%s\\tAA_stat\\t%s\\t%s\\t%s\\n' "${sample_id}" "\$status" "\$rc" "\$message" >> "\$OUTDIR/status.tsv"

    printf '# id: aastat_stats\\n# plot_type: generalstats\\n# pconfig:\\n#   n_mass_shifts:\\n#     title: Mass Shifts\\n#   n_localizations:\\n#     title: Localized bins\\nSample\\tn_mass_shifts\\tn_localizations\\n%s\\t%s\\t%s\\n' \\
        "${sample_id}" "\$n_shifts" "\$n_loc" > "${sample_id}_aastat_mqc.tsv"

    AA_STAT_VERSION=\$("${params.aa_stat_bin}" --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || true)
    AA_STAT_VERSION=\${AA_STAT_VERSION:-unknown}
    cat > versions.yml <<END_VERSIONS
    "${task.process}":
        AA_stat: \$AA_STAT_VERSION
    END_VERSIONS

    exit 0
    """
}
