process FRAGPIPE {
    tag "${sample_id}"
    label 'process_high'
    stageInMode 'copy'

    input:
    tuple val(sample_id), path(raw_input)
    path manifest_template
    path workflow_file

    output:
    tuple val(sample_id), path("*.FPv24hum"), emit: fp_dir
    tuple val(sample_id), path("*_calibrated.mzML"), emit: mzml
    tuple val(sample_id), path("*.FPv24hum/*/*.pepXML"), emit: pepxml
    tuple val(sample_id), path("${sample_id}_fragpipe_mqc.tsv"), emit: mqc
    path "versions.yml", emit: versions

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
        --workdir ${sample_id}.FPv24hum

    psm_count=\$(tail -n +2 ${sample_id}.FPv24hum/*/psm.tsv 2>/dev/null | wc -l || echo 0)
    pep_count=\$(tail -n +2 ${sample_id}.FPv24hum/*/peptide.tsv 2>/dev/null | wc -l || echo 0)
    prot_count=\$(tail -n +2 ${sample_id}.FPv24hum/*/protein.tsv 2>/dev/null | wc -l || echo 0)
    printf '# id: fragpipe_stats\\n# plot_type: generalstats\\n# pconfig:\\n#   psm_count:\\n#     title: PSMs\\n#   peptide_count:\\n#     title: Peptides\\n#   protein_count:\\n#     title: Proteins\\nSample\\tpsm_count\\tpeptide_count\\tprotein_count\\n${sample_id}\\t'\$psm_count'\\t'\$pep_count'\\t'\$prot_count'\\n' > ${sample_id}_fragpipe_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fragpipe: \$(${params.fragpipe_bin} --version 2>&1 | grep -oP 'FragPipe[\\s-]+\\K[\\d.]+' || echo unknown)
    END_VERSIONS
    """
}