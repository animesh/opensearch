process AA_STAT {
    errorStrategy 'ignore'
    tag "${sample_id}"
    label 'process_high'

    publishDir "${params.outdir}/${params.aa_stat_bin.toString().tokenize('/')[-1]}", mode: 'copy'

    input:
    tuple val(sample_id), val(mzml_path), path(pepxml_file)

    output:
    tuple val(sample_id), path("*.${params.aa_stat_bin.toString().tokenize('/')[-1]}"), emit: aa_stat_dir
    tuple val(sample_id), path("${sample_id}_aastat_mqc.tsv"), emit: mqc
    path "versions.yml", emit: versions

    script:
    def tool_name = params.aa_stat_bin.toString().tokenize('/')[-1]
    """
    set -euo pipefail

    OUTDIR="${sample_id}.${params.aa_stat_bin.toString().tokenize('/')[-1]}"
    mkdir -p "\$OUTDIR"
    MZML="${mzml_path}"
    PEPXML="\$(realpath "${pepxml_file}")"
    rc=0
    status="SUCCESS"
    message="${params.aa_stat_bin.toString().tokenize('/')[-1]} completed successfully"

    set +e
    "${params.aa_stat_bin}" -n ${task.cpus} --mzml "\$MZML" --pepxml "\$PEPXML" --dir "\$OUTDIR"
    rc=\$?
    set -e

    if [[ \$rc -ne 0 ]]; then
        status="FAILED"
        message="AA_stat exited with code \$rc. See Nextflow .command.log."
    fi

    n_shifts=0; n_loc=0
    [[ -s "\$OUTDIR/aa_statistics_table.csv" ]] && n_shifts=\$(grep -c '^[+-]' "\$OUTDIR/aa_statistics_table.csv" || true)
    [[ -s "\$OUTDIR/localization_statistics.csv" ]] && n_loc=\$(tail -n +2 "\$OUTDIR/localization_statistics.csv" | wc -l || true)

    printf 'sample\\ttool\\tstatus\\texit_code\\tmessage\\n' > "\$OUTDIR/status.tsv"
    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "${sample_id}" "${params.aa_stat_bin.toString().tokenize('/')[-1]}" "\$status" "\$rc" "\$message" >> "\$OUTDIR/status.tsv"

    printf '# id: aastat_stats\\n# plot_type: generalstats\\n# pconfig:\\n#   n_mass_shifts:\\n#     title: Mass Shifts\\n#   n_localizations:\\n#     title: Localized bins\\nSample\\tn_mass_shifts\\tn_localizations\\n%s\\t%s\\t%s\\n' \\
        "${sample_id}" "\$n_shifts" "\$n_loc" > "${sample_id}_aastat_mqc.tsv"

    printf 'tool\\tsample\\trole\\tpath\\n' > "\$OUTDIR/opensearch_manifest.tsv"
    printf '%s\\t%s\\tstatus\\tstatus.tsv\\n' "${params.aa_stat_bin.toString().tokenize('/')[-1]}" "${sample_id}" >> "\$OUTDIR/opensearch_manifest.tsv"
    [[ -s "\$OUTDIR/report.html" ]] && printf '%s\\t%s\\treport\\treport.html\\n' "${params.aa_stat_bin.toString().tokenize('/')[-1]}" "${sample_id}" >> "\$OUTDIR/opensearch_manifest.tsv"
    [[ -s "\$OUTDIR/aa_statistics_table.csv" ]] && printf '%s\\t%s\\tstats\\taa_statistics_table.csv\\n' "${params.aa_stat_bin.toString().tokenize('/')[-1]}" "${sample_id}" >> "\$OUTDIR/opensearch_manifest.tsv"
    [[ -s "\$OUTDIR/interpretations.json" ]] && printf '%s\\t%s\\tinterpretations\\tinterpretations.json\\n' "${params.aa_stat_bin.toString().tokenize('/')[-1]}" "${sample_id}" >> "\$OUTDIR/opensearch_manifest.tsv"
    [[ -s "\$OUTDIR/localization_statistics.csv" ]] && printf '%s\\t%s\\tlocalizations\\tlocalization_statistics.csv\\n' "${params.aa_stat_bin.toString().tokenize('/')[-1]}" "${sample_id}" >> "\$OUTDIR/opensearch_manifest.tsv"

    AA_STAT_VERSION=\$("${params.aa_stat_bin}" --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || true)
    AA_STAT_VERSION=\${AA_STAT_VERSION:-unknown}
    cat > versions.yml <<END_VERSIONS
    "${task.process}":
        ${params.aa_stat_bin.toString().tokenize('/')[-1]}: \$AA_STAT_VERSION
    END_VERSIONS

    exit 0
    """
}
