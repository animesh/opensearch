process FRAGPIPE {
    tag "${sample_id}"
    label 'process_high'
    // A single problematic sample (for example a blank with insufficient
    // high-confidence PSMs for PeptideProphet/ProteinProphet) must not cancel
    // every other sample. The status.tsv in the result directory is the
    // authoritative per-sample status.
    errorStrategy { params.fragpipe_allow_partial ? 'ignore' : 'terminate' }
    stageInMode params.fragpipe_stage_mode

    input:
    tuple val(sample_id), path(raw_input), val(input_model)
    path workflow_file
    path spectrum_counter
    path workflow_renderer

    output:
    tuple val(sample_id), path("*.${params.fragpipe_workdir_suffix}"), emit: fp_dir
    tuple val(sample_id), val(input_model), path("${sample_id}_analysis.mzML", optional: true), emit: mzml
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
    PYTHON_CONFIG="${params.fragpipe_python}"
    if [[ -f "\$PYTHON_CONFIG" ]]; then PYTHON_CONFIG="\$(dirname "\$PYTHON_CONFIG")"; fi
    RAW_PATH="\$PWD/\$RAW_INPUT"
    FPDIR="${sample_id}.${params.fragpipe_workdir_suffix}"
    GENERATED_MANIFEST="fp.generated.manifest.txt"
    GENERATED_WORKFLOW="fp.generated.workflow"
    FP_RC=0

    mkdir -p "\$FPDIR"

    # Always leave a version file, even when FragPipe itself fails. This keeps
    # Nextflow's optional/failure-isolated reporting deterministic.
    FRAGPIPE_VERSION=\$(${params.fragpipe_bin} --version 2>&1 | grep -oP 'FragPipe[\\s-]+\\K[\\d.]+' | head -1 || true)
    FRAGPIPE_VERSION=\${FRAGPIPE_VERSION:-unknown}
    cat > versions.yml <<END_VERSIONS
    "${task.process}":
        fragpipe: \$FRAGPIPE_VERSION
    END_VERSIONS

    if [[ ! -e "\$RAW_INPUT" ]]; then
        echo "ERROR: staged raw input does not exist: \$RAW_INPUT" >&2
        FP_RC=2
    elif [[ "\$RAW_INPUT" == *.d && ! -d "\$RAW_INPUT" ]]; then
        echo "ERROR: staged timsTOF input is not a directory: \$RAW_INPUT" >&2
        FP_RC=3
    fi

    if [[ \$FP_RC -eq 0 ]]; then
        # One FragPipe task handles one file, so this is a one-run experiment.
        # The four-column manifest format is required by FragPipe headless mode.
        printf '%s\\t%s\\t\\tDDA\\n' "\$RAW_PATH" "\$RAWFILE" > "\$GENERATED_MANIFEST"
        echo "[OpenSearch] Generated FragPipe manifest:"
        cat -A "\$GENERATED_MANIFEST"

        python3 ${workflow_renderer} \\
            --source "${workflow_file}" \\
            --destination "\$GENERATED_WORKFLOW" \\
            --database "${params.fragpipe_database}" \\
            --tools-folder "${params.fragpipe_tools_folder}" \\
            --diann "${params.fragpipe_diann}" \\
            --python "${params.fragpipe_python}" \\
            --im-ms "${is_im_ms ? 'true' : 'false'}" \\
            --crystalc "${crystal_mode}" \\
            --write-calibrated || FP_RC=\$?

        echo "[OpenSearch] FragPipe workflow rendered for ${is_im_ms ? 'Bruker timsTOF IM-MS' : 'regular MS'} input"
        echo "[OpenSearch] Database: ${params.fragpipe_database}"
        echo "[OpenSearch] Tools: ${params.fragpipe_tools_folder}"
        echo "[OpenSearch] DIA-NN: ${params.fragpipe_diann}"
        echo "[OpenSearch] Python: ${params.fragpipe_python}"
        echo "[OpenSearch] Crystal-C mode: ${crystal_mode}"

        if [[ ! -s "${params.fragpipe_database}" ]]; then
            echo "ERROR: FragPipe FASTA database is missing or empty: ${params.fragpipe_database}" >&2
            FP_RC=4
        elif [[ ! -d "${params.fragpipe_tools_folder}" ]]; then
            echo "ERROR: FragPipe tools directory is missing: ${params.fragpipe_tools_folder}" >&2
            FP_RC=5
        elif [[ ! -x "${params.fragpipe_bin}" ]]; then
            echo "ERROR: FragPipe executable is missing or not executable: ${params.fragpipe_bin}" >&2
            FP_RC=6
        fi
    fi

    if [[ \$FP_RC -eq 0 ]]; then
        echo "[OpenSearch] Running FragPipe for ${sample_id} with ${params.cpus} CPUs / ${params.ram_gb} GB RAM"
        ${params.fragpipe_bin} \\
            --headless \\
            --threads ${params.cpus} \\
            --ram ${params.ram_gb} \\
            --workflow "\$GENERATED_WORKFLOW" \\
            --manifest "\$GENERATED_MANIFEST" \\
            --workdir "\$FPDIR" \\
            --config-tools-folder "${params.fragpipe_tools_folder}" \\
            --config-diann "${params.fragpipe_diann}" \\
            --config-python "\$PYTHON_CONFIG" || FP_RC=\$?
    fi

    # FragPipe's Open workflow may fail after MSFragger has already produced a
    # usable pepXML/mzML. Preserve those artifacts so Casanovo/AA_stat can run
    # and the integrated report can distinguish partial from complete results.
    CALIBRATED="${raw_stem}_calibrated.mzML"
    UNCALIBRATED="${raw_stem}_uncalibrated.mzML"
    copy_if_found() {
        local target="\$1"; shift
        local candidate
        for candidate in "\$@"; do
            if [[ -s "\$candidate" ]]; then
                cp -f "\$candidate" "\$target"
                return 0
            fi
        done
        return 1
    }

    if copy_if_found "\$CALIBRATED" \\
        "\$CALIBRATED" \\
        "\$RAW_INPUT/../\$CALIBRATED" \\
        "\$FPDIR/\$CALIBRATED" \\
        "\$(find . -type f -iname "\${RAWFILE}_calibrated.mzML" -size +0c -print -quit 2>/dev/null || true)"; then
        cp -f "\$CALIBRATED" "${sample_id}_analysis.mzML"
        echo "[OpenSearch] Calibrated mzML staged as \$CALIBRATED"
        python3 ${spectrum_counter} "${sample_id}_analysis.mzML" "\$FPDIR/spectrum_count.tsv" || \
            echo "WARNING: spectrum counting failed for \$CALIBRATED" >&2
    elif copy_if_found "\$UNCALIBRATED" \\
        "\$UNCALIBRATED" \\
        "\$RAW_INPUT/../\$UNCALIBRATED" \\
        "\$FPDIR/\$UNCALIBRATED" \\
        "\$(find . -type f -iname "\${RAWFILE}_uncalibrated.mzML" -size +0c -print -quit 2>/dev/null || true)"; then
        cp -f "\$UNCALIBRATED" "${sample_id}_analysis.mzML"
        echo "WARNING: calibrated mzML was unavailable; using uncalibrated mzML for downstream QC/tools" >&2
        python3 ${spectrum_counter} "${sample_id}_analysis.mzML" "\$FPDIR/spectrum_count.tsv" || \
            echo "WARNING: spectrum counting failed for \$UNCALIBRATED" >&2
    else
        echo "WARNING: no calibrated or uncalibrated mzML was found for ${sample_id}" >&2
    fi

    psm_count=\$(find "\$FPDIR" -type f -name psm.tsv -size +0c -print0 2>/dev/null | xargs -0 -r awk 'NR>1 {n++} END {print n+0}')
    pep_count=\$(find "\$FPDIR" -type f -name peptide.tsv -size +0c -print0 2>/dev/null | xargs -0 -r awk 'NR>1 {n++} END {print n+0}')
    prot_count=\$(find "\$FPDIR" -type f -name protein.tsv -size +0c -print0 2>/dev/null | xargs -0 -r awk 'NR>1 {n++} END {print n+0}')
    psm_count=\${psm_count:-0}; pep_count=\${pep_count:-0}; prot_count=\${prot_count:-0}

    status="SUCCESS"; message="FragPipe completed and core tables were produced."
    if [[ \$FP_RC -ne 0 ]]; then
        if [[ \$psm_count -gt 0 || \$pep_count -gt 0 || \$prot_count -gt 0 ]]; then
            status="PARTIAL"
            message="FragPipe exited with code \$FP_RC after producing partial core tables; intermediate artifacts were retained."
        else
            status="FAILED"
            message="FragPipe exited with code \$FP_RC and produced no usable core tables; intermediate artifacts were retained when available."
        fi
    elif [[ \$psm_count -eq 0 && \$pep_count -eq 0 && \$prot_count -eq 0 ]]; then
        status="PARTIAL"
        message="FragPipe returned success but produced no core PSM/peptide/protein tables."
    fi

    printf 'sample\\ttool\\tstatus\\texit_code\\tmessage\\n${sample_id}\\tFragPipe\\t%s\\t%s\\t%s\\n' "\$status" "\$FP_RC" "\$message" > "\$FPDIR/status.tsv"
    printf '# id: fragpipe_stats\\n# plot_type: generalstats\\n# pconfig:\\n#   psm_count:\\n#     title: PSMs\\n#   peptide_count:\\n#     title: Peptides\\n#   protein_count:\\n#     title: Proteins\\nSample\\tpsm_count\\tpeptide_count\\tprotein_count\\n${sample_id}\\t%s\\t%s\\t%s\\n' "\$psm_count" "\$pep_count" "\$prot_count" > "${sample_id}_fragpipe_mqc.tsv"

    # With failure isolation enabled, return success to Nextflow so the
    # partial result directory and status record are emitted. The real tool
    # exit code remains in status.tsv. Set --fragpipe_allow_partial false to
    # propagate the FragPipe exit code and fail the pipeline.
    if [[ "${params.fragpipe_allow_partial}" == "true" ]]; then
        exit 0
    fi
    exit \$FP_RC
    """
}
