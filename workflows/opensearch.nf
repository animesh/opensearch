include { FRAGPIPE } from '../modules/local/fragpipe'
include { CASANOVO } from '../modules/local/casanovo'
include { AA_STAT  } from '../modules/local/aa_stat'
include { MULTIQC  } from '../modules/nf-core/multiqc/main'

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

    def casanovo_available = params.casanovo_bin && new File(params.casanovo_bin as String).canExecute()
    def aa_stat_available  = params.aa_stat_bin && new File(params.aa_stat_bin as String).canExecute()

    if (params.run_casanovo && !casanovo_available) {
        log.warn "CASANOVO requested but binary is missing or not executable: '${params.casanovo_bin}'. Skipping CASANOVO step."
    }
    if (params.run_aa_stat && !aa_stat_available) {
        log.warn "AA_STAT requested but binary is missing or not executable: '${params.aa_stat_bin}'. Skipping AA_STAT step."
    }

    if (params.run_casanovo && casanovo_available) {
        CASANOVO(FRAGPIPE.out.mzml)
        ch_multiqc_files = ch_multiqc_files.mix(CASANOVO.out.mqc.map { it[1] })
        ch_versions = ch_versions.mix(CASANOVO.out.versions)
    }

    if (params.run_aa_stat && aa_stat_available) {
        aa_stat_input = FRAGPIPE.out.mzml.join(FRAGPIPE.out.pepxml, by: 0)
        AA_STAT(aa_stat_input)
        ch_multiqc_files = ch_multiqc_files.mix(AA_STAT.out.mqc.map { it[1] })
        ch_versions = ch_versions.mix(AA_STAT.out.versions)
    }

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
