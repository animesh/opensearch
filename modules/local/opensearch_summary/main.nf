process OPENSEARCH_SUMMARY {
    tag 'integrated summary'
    label 'process_low'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'summary.tsv'

    input:
    path fragpipe_dirs
    path casanovo_dirs
    path aastat_dirs

    output:
    path '*_mqc.json', emit: mqc_json
    path '*_mqc.html', emit: mqc_html
    path 'summary.tsv', emit: summary

    script:
    """
    opensearch_summary.py \\
        --fragpipe ${fragpipe_dirs.join(' ')} \\
        --casanovo ${casanovo_dirs.join(' ')} \\
        --aastat ${aastat_dirs.join(' ')}
    """
}
