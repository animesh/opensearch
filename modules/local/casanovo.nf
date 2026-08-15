process CASANOVO {
    tag "${sample_id}"
    label 'process_medium'

    input:
    tuple val(sample_id), path(mzml_file)

    output:
    tuple val(sample_id), path("${sample_id}.DNv5p1p2"), emit: casanovo_dir
    tuple val(sample_id), path("${sample_id}_casanovo_mqc.tsv"), emit: mqc
    path "versions.yml", emit: versions

    script:
    """
    ${params.casanovo_bin} sequence ${mzml_file} --output_dir ${sample_id}.DNv5p1p2

    log_file=\$(ls ${sample_id}.DNv5p1p2/casanovo_*.log | head -1)
    total=\$(grep -oP 'dataset contains \\K[0-9]+' "\$log_file" || echo 0)
    sequenced=\$(grep -oP 'Sequenced \\K[0-9]+' "\$log_file" || echo 0)
    score50=\$(grep -oP '([0-9]+) spectra \\(.*?\\) scored ≥ 0\\.50' "\$log_file" | grep -oP '^[0-9]+' || echo 0)
    score90=\$(grep -oP '([0-9]+) spectra \\(.*?\\) scored ≥ 0\\.90' "\$log_file" | grep -oP '^[0-9]+' || echo 0)
    printf '# id: casanovo_stats\\n# plot_type: generalstats\\n# pconfig:\\n#   total_spectra:\\n#     title: Total Spectra\\n#   sequenced_spectra:\\n#     title: Sequenced\\n#   score_ge_50pct:\\n#     title: "Score>=0.5"\\n#   score_ge_90pct:\\n#     title: "Score>=0.9"\\nSample\\ttotal_spectra\\tsequenced_spectra\\tscore_ge_50pct\\tscore_ge_90pct\\n${sample_id}\\t'\$total'\\t'\$sequenced'\\t'\$score50'\\t'\$score90'\\n' > ${sample_id}_casanovo_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        casanovo: \$(${params.casanovo_bin} --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo unknown)
    END_VERSIONS
    """
}