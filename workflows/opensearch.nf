include { FRAGPIPE } from '../modules/local/fragpipe'
include { CASANOVO } from '../modules/local/casanovo'
include { AA_STAT  } from '../modules/local/aa_stat'
include { MULTIQC  } from '../modules/nf-core/multiqc/main'
include { OPENSEARCH_SUMMARY } from '../modules/local/opensearch_summary'

workflow OPENSEARCH {
    take:
    raw_ch
    workflow_ch

    main:
    ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_files  = Channel.empty()
    spectrum_counter = file("${projectDir}/bin/count_mzml_spectra.py")
    workflow_renderer = file("${projectDir}/bin/render_fragpipe_workflow.py")

    FRAGPIPE(raw_ch, workflow_ch, spectrum_counter, workflow_renderer)
    ch_versions = FRAGPIPE.out.versions

    // A failed FragPipe post-processing step must not discard a usable mzML or
    // pepXML produced earlier in the same task. FRAGPIPE is intentionally
    // failure-isolated; status.tsv records the exact result state.
    def run_casanovo = params.run_casanovo.toString().toBoolean()
    def run_aa_stat = params.run_aa_stat.toString().toBoolean()
    def casanovo_available = params.casanovo_bin && new File(params.casanovo_bin as String).canExecute()
    def aa_stat_available = params.aa_stat_bin && new File(params.aa_stat_bin as String).canExecute()

    if (run_casanovo && !casanovo_available) {
        log.warn "CASANOVO requested but binary is missing or not executable: '${params.casanovo_bin}'. Skipping CASANOVO."
    }
    if (run_aa_stat && !aa_stat_available) {
        log.warn "AA_STAT requested but binary is missing or not executable: '${params.aa_stat_bin}'. Skipping AA_STAT."
    }

    if (run_casanovo && casanovo_available) {
        CASANOVO(FRAGPIPE.out.mzml.map { sample, model, mzml -> tuple(sample, mzml) })
        ch_versions = ch_versions.mix(CASANOVO.out.versions)
    }

    if (run_aa_stat && aa_stat_available) {
        aa_stat_input = FRAGPIPE.out.mzml.map { sample, model, mzml -> tuple(sample, mzml) }.join(FRAGPIPE.out.pepxml, by: 0)
        AA_STAT(aa_stat_input)
        ch_versions = ch_versions.mix(AA_STAT.out.versions)
    }

    // Integrated reporting receives every result directory, including failed
    // or partial FragPipe directories, so it can report status and provenance.
    summary_dirs = FRAGPIPE.out.fp_dir.map { it[1] }
    if (run_casanovo && casanovo_available) {
        summary_dirs = summary_dirs.mix(CASANOVO.out.casanovo_dir.map { it[1] })
    }
    if (run_aa_stat && aa_stat_available) {
        summary_dirs = summary_dirs.mix(AA_STAT.out.aa_stat_dir.map { it[1] })
    }
    OPENSEARCH_SUMMARY(summary_dirs.collect())
    ch_multiqc_files = ch_multiqc_files.mix(OPENSEARCH_SUMMARY.out.mqc_json)
    ch_multiqc_files = ch_multiqc_files.mix(OPENSEARCH_SUMMARY.out.mqc_html)

    ch_versions
        .collectFile(name: 'software_versions_mqc.yaml', storeDir: "${params.outdir}/pipeline_info")

    def multiqc_config_path = file("${projectDir}/assets/multiqc_config.yml")
    MULTIQC(
        ch_multiqc_files.collect()
            .map { files -> [ [id: 'multiqc'], files, multiqc_config_path, [], [], [] ] }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { meta, report -> report }
}
