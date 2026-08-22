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
    mkdir -p ${sample_id}.${params.aastat_workdir_suffix}
    rc=0
    status="SUCCESS"
    message="AA_stat completed successfully"

    AA_STAT_VERSION=\$(${params.aa_stat_bin} --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || true)
    AA_STAT_VERSION=\${AA_STAT_VERSION:-unknown}
    cat > versions.yml <<END_VERSIONS
    "${task.process}":
        AA_stat: \$AA_STAT_VERSION
    END_VERSIONS

    echo "[OpenSearch] Running AA_stat for ${sample_id}"
    set +e
    ${params.aa_stat_bin} \\
        -n ${task.cpus} \\
        --mzml ${mzml_file} \\
        --pepxml ${pepxml_file} \\
        --dir ${sample_id}.${params.aastat_workdir_suffix}
    rc=\$?
    set -e

    if [[ \$rc -ne 0 ]]; then
        status="FAILED"
        message="AA_stat exited with code \$rc. See the Nextflow .command.log."
        echo "WARNING: ${sample_id}: \$message" >&2
    fi

    n_shifts=\$(grep -c '^[+-]' ${sample_id}.${params.aastat_workdir_suffix}/aa_statistics_table.csv 2>/dev/null || true)
    n_loc=\$(tail -n +2 ${sample_id}.${params.aastat_workdir_suffix}/localization_statistics.csv 2>/dev/null | wc -l || true)
    printf 'sample\\ttool\\tstatus\\texit_code\\tmessage\\n${sample_id}\\tAA_stat\\t%s\\t%s\\t%s\\n' "\$status" "\$rc" "\$message" > ${sample_id}.${params.aastat_workdir_suffix}/status.tsv
    printf '# id: aastat_stats\\n# plot_type: generalstats\\n# pconfig:\\n#   n_mass_shifts:\\n#     title: Mass Shifts\\n#   n_localizations:\\n#     title: Localized bins\\nSample\\tn_mass_shifts\\tn_localizations\\n${sample_id}\\t'\$n_shifts'\\t'\$n_loc'\\n' > ${sample_id}_aastat_mqc.tsv

    # Optional downstream tools are failure-isolated; status.tsv retains the
    # real exit code so the integrated report can distinguish failure from no data.
    exit 0
    """
}
