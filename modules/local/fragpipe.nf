process FRAGPIPE {
    tag "${sample_id}"
    label 'process_high'

    // FragPipe is the primary analysis, but one bad sample must not stop the run.
    // The real FragPipe exit code is recorded in status.tsv.
    errorStrategy 'ignore'

    // IMPORTANT: copy the complete input into the task directory.
    // This is deliberately simple and works for both Thermo .raw and Bruker .d.
    stageInMode 'copy'

    input:
    tuple val(sample_id), path(raw_input), val(input_model)
    path workflow_file
    path spectrum_counter
    path workflow_renderer

    output:
    tuple val(sample_id), path("*.${params.fragpipe_workdir_suffix}"), emit: fp_dir
    tuple val(sample_id), val(input_model), path("${sample_id}_calibrated.mzML", optional: true), emit: mzml
    tuple val(sample_id), path("*.${params.fragpipe_workdir_suffix}/*/*.pepXML", optional: true), emit: pepxml
    tuple val(sample_id), path("${sample_id}_fragpipe_mqc.tsv"), emit: mqc
    path "versions.yml", emit: versions

    script:
    raw_name = raw_input.name
    raw_stem = raw_name.replaceFirst(/\.[^.]+$/, '')
    is_im_ms = input_model == 'timstof'
    crystal_mode = params.fragpipe_crystalc == null ? 'auto' : params.fragpipe_crystalc.toString()
    """
    set -euo pipefail

    RAW_INPUT="${raw_name}"
    RAWFILE="${raw_stem}"
    RAW_PATH="\$(realpath "\$RAW_INPUT")"
    FPDIR="${sample_id}.${params.fragpipe_workdir_suffix}"
    GENERATED_MANIFEST="fp.generated.manifest.txt"
    GENERATED_WORKFLOW="fp.generated.workflow"
    FP_RC=0

    # FragPipe headless mode expects the four-column manifest.
    printf '%s\\t%s\\t\\tDDA\\n' "\$RAW_PATH" "\$RAWFILE" > "\$GENERATED_MANIFEST"

    # Render only the environment/input-specific FragPipe settings. The supplied
    # workflow remains the source of truth for all search parameters.
    python3 "${workflow_renderer}" \\
        --source "${workflow_file}" \\
        --destination "\$GENERATED_WORKFLOW" \\
        --database "${params.fragpipe_database}" \\
        --tools-folder "${params.fragpipe_tools_folder}" \\
        --diann "${params.fragpipe_diann}" \\
        --python "${params.fragpipe_python}" \\
        --im-ms "${is_im_ms ? 'true' : 'false'}" \\
        --crystalc "${crystal_mode}" \\
        --write-calibrated || FP_RC=\$?

    if [[ \$FP_RC -eq 0 ]]; then
        if [[ ! -s "${params.fragpipe_database}" ]]; then
            echo "ERROR: FASTA database missing or empty: ${params.fragpipe_database}" >&2
            FP_RC=2
        elif [[ ! -d "${params.fragpipe_tools_folder}" ]]; then
            echo "ERROR: FragPipe tools directory missing: ${params.fragpipe_tools_folder}" >&2
            FP_RC=3
        elif [[ ! -x "${params.fragpipe_bin}" ]]; then
            echo "ERROR: FragPipe executable missing/not executable: ${params.fragpipe_bin}" >&2
            FP_RC=4
        fi
    fi

    if [[ \$FP_RC -eq 0 ]]; then
        echo "[OpenSearch] FragPipe: ${sample_id} (${input_model})"
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

    # FragPipe writes the calibrated mzML beside the input in the task directory.
    # Find it by exact sample name first, then fall back to any calibrated mzML.
    CALIBRATED="${sample_id}_calibrated.mzML"
    FOUND_MZML=""
    if [[ -s "\$CALIBRATED" ]]; then
        FOUND_MZML="\$CALIBRATED"
    else
        FOUND_MZML="\$(find "\$PWD" -type f -iname "${raw_stem}_calibrated.mzML" -size +0c -print -quit 2>/dev/null || true)"
    fi

    if [[ -n "\$FOUND_MZML" && "\$FOUND_MZML" != "\$CALIBRATED" ]]; then
        cp -f "\$FOUND_MZML" "\$CALIBRATED"
    fi

    if [[ -s "\$CALIBRATED" ]]; then
        echo "[OpenSearch] Calibrated mzML: \$(realpath "\$CALIBRATED")"
        if ! python3 "${spectrum_counter}" "\$(realpath "\$CALIBRATED")" "\$FPDIR/spectrum_count.tsv"; then
            echo "WARNING: could not count spectra in calibrated mzML; FragPipe result remains usable." >&2
            rm -f "\$FPDIR/spectrum_count.tsv"
        fi
    else
        echo "WARNING: no calibrated mzML produced for ${sample_id}" >&2
    fi

    psm_count=\$(find "\$FPDIR" -type f -name psm.tsv -size +0c -print0 2>/dev/null | xargs -0 -r awk 'NR>1 {n++} END {print n+0}')
    pep_count=\$(find "\$FPDIR" -type f -name peptide.tsv -size +0c -print0 2>/dev/null | xargs -0 -r awk 'NR>1 {n++} END {print n+0}')
    prot_count=\$(find "\$FPDIR" -type f -name protein.tsv -size +0c -print0 2>/dev/null | xargs -0 -r awk 'NR>1 {n++} END {print n+0}')
    psm_count=\${psm_count:-0}; pep_count=\${pep_count:-0}; prot_count=\${prot_count:-0}

    status="SUCCESS"
    message="FragPipe completed and core tables were produced."
    if [[ \$FP_RC -ne 0 ]]; then
        status="FAILED"
        message="FragPipe exited with code \$FP_RC. Partial FragPipe artifacts were retained when available."
    elif [[ \$psm_count -eq 0 && \$pep_count -eq 0 && \$prot_count -eq 0 ]]; then
        status="PARTIAL"
        message="FragPipe returned success but produced no core PSM/peptide/protein tables."
    fi

    printf 'sample\\ttool\\tstatus\\texit_code\\tmessage\\n' > "\$FPDIR/status.tsv"
    printf '%s\\tFragPipe\\t%s\\t%s\\t%s\\n' "${sample_id}" "\$status" "\$FP_RC" "\$message" >> "\$FPDIR/status.tsv"

    printf '# id: fragpipe_stats\\n# plot_type: generalstats\\n# pconfig:\\n#   psm_count:\\n#     title: PSMs\\n#   peptide_count:\\n#     title: Peptides\\n#   protein_count:\\n#     title: Proteins\\nSample\\tpsm_count\\tpeptide_count\\tprotein_count\\n%s\\t%s\\t%s\\t%s\\n' \\
        "${sample_id}" "\$psm_count" "\$pep_count" "\$prot_count" > "${sample_id}_fragpipe_mqc.tsv"

    cat > versions.yml <<END_VERSIONS
    "${task.process}":
        fragpipe: \$(${params.fragpipe_bin} --version 2>&1 | grep -oP 'FragPipe[\\s-]+\\K[\\d.]+' | head -1 || echo unknown)
    END_VERSIONS

    # Always emit the task products so the integrated report can distinguish
    # complete, partial and failed samples without stopping the whole batch.
    exit 0
    """
}
