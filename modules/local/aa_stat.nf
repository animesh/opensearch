process AA_STAT {
    tag "${sample_id}"
    label 'process_medium'

    input:
    tuple val(sample_id), path(mzml_file), path(pepxml_file)

    output:
    tuple val(sample_id), path("${sample_id}.AA_stat_v2p5p6hum"), emit: aa_stat_dir

    script:
    """
    ${params.aa_stat_bin} \
        -n ${task.cpus} \
        --mzml ${mzml_file} \
        --pepxml ${pepxml_file} \
        --dir ${sample_id}.AA_stat_v2p5p6hum
    """
}