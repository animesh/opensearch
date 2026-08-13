process FRAGPIPE {
    tag "${sample_id}"
    label 'process_high'
    stageInMode 'copy'

    input:
    tuple val(sample_id), path(raw_input)
    path manifest_template
    path workflow_file

    output:
    tuple val(sample_id), path("*.FPv22hum"), emit: fp_dir
    tuple val(sample_id), path("*_calibrated.mzML"), emit: mzml
    tuple val(sample_id), path("*.FPv22hum/*/*.pepXML"), emit: pepxml

    script:
    raw_name = raw_input.name
    raw_stem = raw_name.replaceFirst(/\.[^.]+$/, '')
    """
    RAWDIR="\$(realpath ${raw_input})"
    RAWFILE="${raw_stem}"
    
    awk -v RAWDIR="\$RAWDIR" -v RAWFILE="\$RAWFILE" \
        '{ gsub("RAWDIR",RAWDIR); gsub("RAWFILE",RAWFILE); print }' \
        ${manifest_template} > fp.generated.manifest.txt

    ${params.fragpipe_bin} \
        --headless \
        --threads ${params.fragpipe_threads} \
        --ram ${params.fragpipe_ram_gb} \
        --workflow ${workflow_file} \
        --manifest fp.generated.manifest.txt \
        --workdir ${sample_id}.FPv22hum
    """
}