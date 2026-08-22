process CASANOVO {
    errorStrategy 'ignore'
    tag "${sample_id}"
    label 'process_high'

    input:
    tuple val(sample_id), path(mzml_file)

    output:
    tuple val(sample_id), path("${sample_id}.${params.casanovo_workdir_suffix}"), emit: casanovo_dir
    tuple val(sample_id), path("${sample_id}_casanovo_mqc.tsv"), emit: mqc
    path "versions.yml", emit: versions

    script:
    def casanovo_dir = "${sample_id}.${params.casanovo_workdir_suffix}"
    def mzml_name = mzml_file.getName()

    """
    set -euo pipefail

    mkdir -p "${casanovo_dir}"

    CASANOVO_VERSION=\$(${params.casanovo_bin} --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || true)
    CASANOVO_VERSION=\${CASANOVO_VERSION:-unknown}

    echo "[OpenSearch] Running Casanovo for ${sample_id} using the default model"
    echo "[OpenSearch] Input: ${mzml_name}"
    echo "[OpenSearch] Output: ${casanovo_dir}"

    rc=0
    set +e

    ${params.casanovo_bin} sequence \
        "${mzml_name}" \
        --output_dir "${casanovo_dir}"

    rc=\$?
    set -e

    status="SUCCESS"
    message="Casanovo completed successfully using the default model"

    if [[ \$rc -ne 0 ]]; then
        status="FAILED"
        message="Casanovo exited with code \$rc; see the Casanovo log."
        echo "WARNING: ${sample_id}: \$message" >&2
    fi

    mztab_file=\$(find "${casanovo_dir}" \
        -maxdepth 1 \
        -type f \
        -name '*.mztab' \
        -size +0c \
        -print -quit)

    if [[ "\$status" == "SUCCESS" && -z "\$mztab_file" ]]; then
        status="FAILED"
        message="Casanovo returned success but produced no usable mzTab result."
        echo "WARNING: ${sample_id}: \$message" >&2
    fi

    printf 'sample\\ttool\\tstatus\\texit_code\\tmessage\\n' > status.tsv
    printf '%s\\tCasanovo\\t%s\\t%s\\t%s\\n' \
        "${sample_id}" "\$status" "\$rc" "\$message" >> status.tsv

    sequenced=0
    score50=0
    score90=0

    printf '# id: casanovo_stats
# plot_type: generalstats
# pconfig:
#   sequenced_spectra:
#     title: Sequenced
#   score_ge_50pct:
#     title: "Score>=0.5"
#   score_ge_90pct:
#     title: "Score>=0.9"
Sample\\tsequenced_spectra\\tscore_ge_50pct\\tscore_ge_90pct
${sample_id}\\t\$sequenced\\t\$score50\\t\$score90
' > "${sample_id}_casanovo_mqc.tsv"

    cat > versions.yml <<EOF
"${task.process}":
    casanovo: \$CASANOVO_VERSION
EOF

    exit 0
    """
}