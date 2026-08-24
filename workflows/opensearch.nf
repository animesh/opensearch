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

    FRAGPIPE(raw_ch, workflow_ch, spectrum_counter, workflow_renderer)

    // FRAGPIPE creates a tiny marker containing the absolute path of the
    // calibrated mzML beside the copied input in the launch directory.
    // Downstream tools receive that absolute path as a value, not as a staged
    // Nextflow file.
    calibrated = FRAGPIPE.out.mzml_path
        .map { sample, marker -> tuple(sample, marker.text.trim()) }
        .filter { sample, mzml -> mzml && new File(mzml).isFile() && new File(mzml).length() > 0 }

    if (params.run_casanovo.toString().toBoolean() && params.casanovo_bin) {
        CASANOVO(calibrated)
    }

    if (params.run_aa_stat.toString().toBoolean() && params.aa_stat_bin) {
        aa_input = calibrated
            .join(FRAGPIPE.out.pepxml, by: 0)
            .map { sample, mzml, pepxml -> tuple(sample, mzml, pepxml) }
        AA_STAT(aa_input)
    }

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
