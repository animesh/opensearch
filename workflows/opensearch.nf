include { FRAGPIPE } from '../modules/local/fragpipe'
include { CASANOVO } from '../modules/local/casanovo'
include { AA_STAT  } from '../modules/local/aa_stat'

workflow OPENSEARCH {
    take:
    raw_ch
    manifest_ch
    workflow_ch

    main:
    FRAGPIPE(raw_ch, manifest_ch, workflow_ch)

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
    }

    if (params.run_aa_stat && aa_stat_available) {
        aa_stat_input = FRAGPIPE.out.mzml.join(FRAGPIPE.out.pepxml, by: 0)
        AA_STAT(aa_stat_input)
    }

    emit:
    // Placeholder until a project-specific MultiQC aggregation step is added.
    multiqc_report = Channel.empty()
}