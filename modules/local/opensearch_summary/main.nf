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

    script:
    """
    opensearch_summary.py \
        --dirs ${result_dirs.join(' ')} \
        --fragpipe_suffix .${params.fragpipe_workdir_suffix} \
        --casanovo_suffix .${params.casanovo_workdir_suffix} \
        --aastat_suffix .${params.aastat_workdir_suffix}
    """
}
