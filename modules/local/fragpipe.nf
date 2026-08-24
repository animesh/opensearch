process FRAGPIPE {
    tag "${sample_id}"
    label 'process_high'
    errorStrategy 'ignore'

    publishDir "${params.outdir}/${params.fragpipe_bin.toString().tokenize('/')[-1]}", mode: 'copy'

    input:
    tuple val(sample_id), val(raw_source), val(input_model)
    path workflow_file
    path spectrum_counter
    path workflow_renderer

    output:
    tuple val(sample_id), path("*.${params.fragpipe_bin.toString().tokenize('/')[-1]}"), emit: fp_dir
    tuple val(sample_id), path("*_calibrated.path"), emit: mzml_path
    tuple val(sample_id), path("*.${params.fragpipe_bin.toString().tokenize('/')[-1]}/*/*.pepXML", optional: true), emit: pepxml
    tuple val(sample_id), path("${sample_id}_fragpipe_mqc.tsv"), emit: mqc
    path "versions.yml", emit: versions

    script:
    def launch_dir = workflow.launchDir.toString()
    def tool_name = params.fragpipe_bin.toString().tokenize('/')[-1]
    """
    set -euo pipefail

    RAW_SOURCE="${raw_source}"
    RAW_NAME=\$(basename "\$RAW_SOURCE")
    RAW_STEM=\"\${RAW_NAME%.*}\"
    INPUT_COPY="${launch_dir}/\$RAW_NAME"
    FPDIR="${sample_id}.${params.fragpipe_bin.toString().tokenize('/')[-1]}"
    GENERATED_MANIFEST="fp.generated.manifest.txt"
    GENERATED_WORKFLOW="fp.generated.workflow"
    FP_RC=0

    if [[ "\$RAW_SOURCE" != "\$INPUT_COPY" ]]; then
        rm -rf "\$INPUT_COPY"
        cp -a "\$RAW_SOURCE" "\$INPUT_COPY"
    fi
    [[ -e "\$INPUT_COPY" ]] || { echo "ERROR: copied input does not exist: \$INPUT_COPY" >&2; FP_RC=2; }

    printf '%s\\t%s\\t\\tDDA\\n' "\$INPUT_COPY" "\$RAW_STEM" > "\$GENERATED_MANIFEST"

    if [[ \$FP_RC -eq 0 ]]; then
        python3 "${workflow_renderer}" \\
            --source "${workflow_file}" \\
            --destination "\$GENERATED_WORKFLOW" \\
            --database "${params.fragpipe_database}" \\
            --tools-folder "${params.fragpipe_tools_folder}" \\
            --diann "${params.fragpipe_diann}" \\
            --python "${params.fragpipe_python}" \\
            --im-ms "${input_model == 'timstof' ? 'true' : 'false'}" \\
            --crystalc "${params.fragpipe_crystalc == null ? 'auto' : params.fragpipe_crystalc.toString()}" \\
            --write-calibrated || FP_RC=\$?
    fi

    if [[ \$FP_RC -eq 0 ]]; then
        [[ -s "${params.fragpipe_database}" ]] || { echo "ERROR: FASTA database missing: ${params.fragpipe_database}" >&2; FP_RC=3; }
        [[ -d "${params.fragpipe_tools_folder}" ]] || { echo "ERROR: FragPipe tools directory missing: ${params.fragpipe_tools_folder}" >&2; FP_RC=4; }
        [[ -x "${params.fragpipe_bin}" ]] || { echo "ERROR: FragPipe executable missing: ${params.fragpipe_bin}" >&2; FP_RC=5; }
    fi

    if [[ \$FP_RC -eq 0 ]]; then
        echo "[OpenSearch] Running FragPipe for ${sample_id} (${input_model})"
        "${params.fragpipe_bin}" \\
            --headless \\
            --threads ${params.cpus} \\
            --ram ${params.ram_gb} \\
            --workflow "\$GENERATED_WORKFLOW" \\
            --manifest "\$GENERATED_MANIFEST" \\
            --workdir "\$FPDIR" \\
            --config-tools-folder "${params.fragpipe_tools_folder}" \\
            --config-diann "${params.fragpipe_diann}" \\
            --config-python "${params.fragpipe_python}" || FP_RC=\$?
    fi

    OUTPUT_NAME=\$(cut -f2 "\$GENERATED_MANIFEST")
    FP_RESULTS="\$FPDIR/\$OUTPUT_NAME"
    CALIBRATED="${launch_dir}/\${RAW_STEM}_calibrated.mzML"
    printf '%s\\n' "\$CALIBRATED" > "${sample_id}_calibrated.path"

    if [[ -s "\$CALIBRATED" ]]; then
        python3 "${spectrum_counter}" "\$CALIBRATED" "\$FPDIR/spectrum_count.tsv" || true
    fi

    PSM="\$FP_RESULTS/psm.tsv"
    PEPTIDE="\$FP_RESULTS/peptide.tsv"
    PROTEIN="\$FP_RESULTS/protein.tsv"
    PTM="\$FPDIR/ptm-shepherd-output/global.modsummary.tsv"
    SPECTRA="\$FPDIR/spectrum_count.tsv"

    psm_count=0; pep_count=0; prot_count=0
    [[ -s "\$PSM" ]] && psm_count=\$(awk 'NR>1 {n++} END {print n+0}' "\$PSM")
    [[ -s "\$PEPTIDE" ]] && pep_count=\$(awk 'NR>1 {n++} END {print n+0}' "\$PEPTIDE")
    [[ -s "\$PROTEIN" ]] && prot_count=\$(awk 'NR>1 {n++} END {print n+0}' "\$PROTEIN")

    status="SUCCESS"
    message="FragPipe completed and core tables were produced."
    if [[ \$FP_RC -ne 0 ]]; then
        status="FAILED"
        message="FragPipe exited with code \$FP_RC. Partial artifacts were retained."
    elif [[ \$psm_count -eq 0 && \$pep_count -eq 0 && \$prot_count -eq 0 ]]; then
        status="PARTIAL"
        message="FragPipe returned success but produced no core PSM/peptide/protein tables."
    fi
    [[ -s "\$CALIBRATED" ]] || message="\$message Calibrated mzML was not produced."

    printf 'sample\\ttool\\tstatus\\texit_code\\tmessage\\n' > "\$FPDIR/status.tsv"
    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "${sample_id}" "${params.fragpipe_bin.toString().tokenize('/')[-1]}" "\$status" "\$FP_RC" "\$message" >> "\$FPDIR/status.tsv"

    printf '# id: fragpipe_stats\\n# plot_type: generalstats\\n# pconfig:\\n#   psm_count:\\n#     title: PSMs\\n#   peptide_count:\\n#     title: Peptides\\n#   protein_count:\\n#     title: Proteins\\nSample\\tpsm_count\\tpeptide_count\\tprotein_count\\n%s\\t%s\\t%s\\t%s\\n' \\
        "${sample_id}" "\$psm_count" "\$pep_count" "\$prot_count" > "${sample_id}_fragpipe_mqc.tsv"

    printf 'tool\\tsample\\trole\\tpath\\n' > "\$FPDIR/opensearch_manifest.tsv"
    printf '%s\\t%s\\tstatus\\tstatus.tsv\\n' "${params.fragpipe_bin.toString().tokenize('/')[-1]}" "${sample_id}" >> "\$FPDIR/opensearch_manifest.tsv"
    [[ -s "\$PSM" ]] && printf '%s\\t%s\\tpsm\\t%s\\n' "${params.fragpipe_bin.toString().tokenize('/')[-1]}" "${sample_id}" "\${PSM#\$FPDIR/}" >> "\$FPDIR/opensearch_manifest.tsv"
    [[ -s "\$PEPTIDE" ]] && printf '%s\\t%s\\tpeptide\\t%s\\n' "${params.fragpipe_bin.toString().tokenize('/')[-1]}" "${sample_id}" "\${PEPTIDE#\$FPDIR/}" >> "\$FPDIR/opensearch_manifest.tsv"
    [[ -s "\$PROTEIN" ]] && printf '%s\\t%s\\tprotein\\t%s\\n' "${params.fragpipe_bin.toString().tokenize('/')[-1]}" "${sample_id}" "\${PROTEIN#\$FPDIR/}" >> "\$FPDIR/opensearch_manifest.tsv"
    [[ -s "\$PTM" ]] && printf '%s\\t%s\\tptm\\t%s\\n' "${params.fragpipe_bin.toString().tokenize('/')[-1]}" "${sample_id}" "\${PTM#\$FPDIR/}" >> "\$FPDIR/opensearch_manifest.tsv"
    [[ -s "\$SPECTRA" ]] && printf '%s\\t%s\\tspectrum_count\\t%s\\n' "${params.fragpipe_bin.toString().tokenize('/')[-1]}" "${sample_id}" "\${SPECTRA#\$FPDIR/}" >> "\$FPDIR/opensearch_manifest.tsv"

    FP_VERSION=\$("${params.fragpipe_bin}" --version 2>&1 | grep -oP 'FragPipe[\\s-]+\\K[\\d.]+' | head -1 || true)
    FP_VERSION=\${FP_VERSION:-unknown}
    cat > versions.yml <<END_VERSIONS
    "${task.process}":
        ${params.fragpipe_bin.toString().tokenize('/')[-1]}: \$FP_VERSION
    END_VERSIONS

    exit 0
    """
}
