process CASANOVO {
    tag "${sample_id}"
    label 'process_medium'

    input:
    tuple val(sample_id), path(mzml_file)

    output:
    tuple val(sample_id), path("${sample_id}.DNv5p1p2"), emit: casanovo_dir

    script:
    """
    ${params.casanovo_bin} sequence ${mzml_file} --output_dir ${sample_id}.DNv5p1p2
    """
}