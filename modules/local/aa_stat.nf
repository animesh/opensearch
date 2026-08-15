process AA_STAT {
    tag "${sample_id}"
    label 'process_medium'

    input:
    tuple val(sample_id), path(mzml_file), path(pepxml_file)

    output:
    tuple val(sample_id), path("${sample_id}.AA_stat_v2p5p6hum"), emit: aa_stat_dir
    tuple val(sample_id), path("${sample_id}_aastat_mqc.tsv"), emit: mqc
    path "versions.yml", emit: versions

    script:
    """
    ${params.aa_stat_bin} \
        -n ${task.cpus} \
        --mzml ${mzml_file} \
        --pepxml ${pepxml_file} \
        --dir ${sample_id}.AA_stat_v2p5p6hum

    # Extract number of mass shifts and recommended variable mods for general stats
    n_shifts=\$(grep -c '^[+-]' ${sample_id}.AA_stat_v2p5p6hum/aa_statistics_table.csv 2>/dev/null || echo 0)
    n_loc=\$(tail -n +2 ${sample_id}.AA_stat_v2p5p6hum/localization_statistics.csv 2>/dev/null | wc -l || echo 0)
    printf '# id: aastat_stats\\n# plot_type: generalstats\\n# pconfig:\\n#   n_mass_shifts:\\n#     title: Mass Shifts\\n#   n_localizations:\\n#     title: Localized Shifts\\nSample\\tn_mass_shifts\\tn_localizations\\n${sample_id}\\t'\$n_shifts'\\t'\$n_loc'\\n' > ${sample_id}_aastat_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AA_stat: \$(${params.aa_stat_bin} --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo unknown)
    END_VERSIONS
    """
}