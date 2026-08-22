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
    spectrum_counter = file("${projectDir}/bin/count_mzml_spectra.py", checkIfExists: true)
    workflow_renderer = file("${projectDir}/bin/render_fragpipe_workflow.py", checkIfExists: true)
    multiqc_config_path = file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)

    // One simple primary task: copy the input, run FragPipe, and expose the
    // calibrated mzML and pepXML produced in that same task.
    FRAGPIPE(raw_ch, workflow_ch, spectrum_counter, workflow_renderer)

    // Only samples for which a calibrated mzML actually exists enter Casanovo.
    // Casanovo itself is deliberately called without --model, so its installed
    // default model is used for both Orbitrap and timsTOF.
    if (params.run_casanovo.toString().toBoolean() && params.casanovo_bin) {
        cas_mzml = FRAGPIPE.out.mzml.filter { it.size() >= 3 && it[2] != null }
        CASANOVO(cas_mzml.map { sample, model, mzml -> tuple(sample, mzml) })
    }

    // AA_stat needs both the calibrated mzML and pepXML from FragPipe.
    if (params.run_aa_stat.toString().toBoolean() && params.aa_stat_bin) {
        aa_mzml = FRAGPIPE.out.mzml.filter { it.size() >= 3 && it[2] != null }
        aa_input = aa_mzml.map { sample, model, mzml -> tuple(sample, mzml) }
            .join(FRAGPIPE.out.pepxml, by: 0)
            .map { sample, mzml, pepxml -> tuple(sample, mzml, pepxml) }
        AA_STAT(aa_input)
    }

    // Summary receives every FragPipe result directory plus successful/failed
    // optional-tool directories. Missing optional directories therefore appear
    // explicitly as NOT_RUN in the integrated report.
    summary_dirs = FRAGPIPE.out.fp_dir.map { it[1] }
    if (params.run_casanovo.toString().toBoolean() && params.casanovo_bin) {
        summary_dirs = summary_dirs.mix(CASANOVO.out.casanovo_dir.map { it[1] })
    }
    if (params.run_aa_stat.toString().toBoolean() && params.aa_stat_bin) {
        summary_dirs = summary_dirs.mix(AA_STAT.out.aa_stat_dir.map { it[1] })
    }

    OPENSEARCH_SUMMARY(summary_dirs.collect())

    multiqc_inputs = OPENSEARCH_SUMMARY.out.mqc_json
        .mix(OPENSEARCH_SUMMARY.out.mqc_html)
        .collect()
        .map { files -> [[id: 'multiqc'], files, multiqc_config_path, [], [], []] }

    MULTIQC(multiqc_inputs)


    emit:
    multiqc_report = MULTIQC.out.report.map { meta, report -> report }
}
