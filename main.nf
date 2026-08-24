#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/opensearch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/opensearch
    Website: https://nf-co.re/opensearch
    Slack  : https://nfcore.slack.com/channels/opensearch
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { OPENSEARCH  } from './workflows/opensearch'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_opensearch_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_opensearch_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_OPENSEARCH {

    take:
    samplesheet // channel: tuple(sample_id, raw_input_path, input_model)

    main:

    required = [
        fragpipe_workflow: params.fragpipe_workflow,
        fragpipe_database: params.fragpipe_database,
        fragpipe_tools_folder: params.fragpipe_tools_folder,
        fragpipe_diann: params.fragpipe_diann,
        fragpipe_bin: params.fragpipe_bin,
    ]
    if (params.run_casanovo.toString().toBoolean()) required.casanovo_bin = params.casanovo_bin
    if (params.run_aa_stat.toString().toBoolean()) required.aa_stat_bin = params.aa_stat_bin
    missing = required.findAll { key, value -> !value }.keySet()
    if (missing) {
        error "Missing required tool configuration: ${missing.join(', ')}"
    }

    workflow_ch = Channel.fromPath(params.fragpipe_workflow, checkIfExists: true).first()

    //
    // WORKFLOW: Run pipeline
    //
    OPENSEARCH (
        samplesheet,
        workflow_ch,
    )
    emit:
    multiqc_report = OPENSEARCH.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_OPENSEARCH (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_OPENSEARCH.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
