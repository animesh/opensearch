include { FRAGPIPE } from '../modules/local/fragpipe'
include { CASANOVO } from '../modules/local/casanovo'
include { AA_STAT  } from '../modules/local/aa_stat'
include { MULTIQC  } from '../modules/nf-core/multiqc/main'
include { OPENSEARCH_SUMMARY } from '../modules/local/opensearch_summary'

workflow OPENSEARCH {
    take:
    raw_ch
    manifest_ch
    workflow_ch

    main:
    ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_files  = Channel.empty()

    FRAGPIPE(raw_ch, manifest_ch, workflow_ch)
    ch_multiqc_files = ch_multiqc_files.mix(FRAGPIPE.out.mqc.map { it[1] })
    ch_versions = FRAGPIPE.out.versions

    // Coerce CLI-supplied flags to proper Groovy booleans.
    // Without this, --run_aa_stat false on the command line arrives as the
    // *string* "false", which is truthy in Groovy and bypasses every guard.
    def run_casanovo    = params.run_casanovo.toString().toBoolean()
    def run_aa_stat     = params.run_aa_stat.toString().toBoolean()
    def casanovo_available = params.casanovo_bin && new File(params.casanovo_bin as String).canExecute()
    def aa_stat_available  = params.aa_stat_bin && new File(params.aa_stat_bin as String).canExecute()

    if (run_casanovo && !casanovo_available) {
        log.warn "CASANOVO requested but binary is missing or not executable: '${params.casanovo_bin}'. Skipping CASANOVO step."
    }
    if (run_aa_stat && !aa_stat_available) {
        log.warn "AA_STAT requested but binary is missing or not executable: '${params.aa_stat_bin}'. Skipping AA_STAT step."
    }

    if (run_casanovo && casanovo_available) {
        CASANOVO(FRAGPIPE.out.mzml)
        ch_multiqc_files = ch_multiqc_files.mix(CASANOVO.out.mqc.map { it[1] })
        ch_versions = ch_versions.mix(CASANOVO.out.versions)
    }

    if (run_aa_stat && aa_stat_available) {
        aa_stat_input = FRAGPIPE.out.mzml.join(FRAGPIPE.out.pepxml, by: 0)
        AA_STAT(aa_stat_input)
        ch_multiqc_files = ch_multiqc_files.mix(AA_STAT.out.mqc.map { it[1] })
        ch_versions = ch_versions.mix(AA_STAT.out.versions)
    }

    // Build one integrated summary from whichever result directories exist.
    // This makes --run_casanovo false / --run_aa_stat false work naturally.
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

    // Collect all versions into a single file for MultiQC
    ch_versions
        .collectFile(name: 'software_versions_mqc.yaml', storeDir: "${params.outdir}/pipeline_info")

    // Build the MULTIQC input tuple without .combine() to avoid LinkedList flattening.
    // The module signature is a single tuple:
    //   tuple val(meta), path(multiqc_files), path(multiqc_config), path(logo), path(replace_names), path(sample_names)
    // We embed the config path directly via file() inside the map closure.
    def multiqc_config_path = file("${projectDir}/assets/multiqc_config.yml")
    MULTIQC(
        ch_multiqc_files.collect()
            .map { files -> [ [id: 'multiqc'], files, multiqc_config_path, [], [], [] ] }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { meta, report -> report }
}
