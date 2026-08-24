process OPENSEARCH_SUMMARY {
    tag 'integrated summary'
    label 'process_low'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: '*.tsv'

    input:
    path result_dirs

    output:
    path '*_mqc.json', emit: mqc_json
    path '*_mqc.html', emit: mqc_html
    path 'summary.tsv', emit: summary
    path 'provenance.tsv', emit: provenance
    path 'peptide_comparison.tsv', emit: peptide_comparison, optional: true
    path 'peptide_comparison_summary.tsv', emit: peptide_comparison_summary, optional: true

    script:
    def fragpipe_dir_name = new File(params.fragpipe_bin.toString()).name
    def casanovo_dir_name = new File(params.casanovo_bin.toString()).name
    def aastat_dir_name = new File(params.aa_stat_bin.toString()).name
    """
    opensearch_summary.py \\
        --dirs ${result_dirs.join(' ')} \\
        --fragpipe_dir_name ${fragpipe_dir_name} \\
        --casanovo_dir_name ${casanovo_dir_name} \\
        --aastat_dir_name ${aastat_dir_name}
    """
}
