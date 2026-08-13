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
    samplesheet // channel: tuple(sample_id, raw_input_path)

    main:

    if (!params.fragpipe_manifest && !params.scripts_dir) {
        error "Please provide FragPipe inputs via either --scripts_dir <dir-with-fp.manifest.txt-and-fp.dl.workflow.txt> or --fragpipe_manifest <file> plus --fragpipe_workflow <file>."
    }
    if (params.fragpipe_manifest && !params.fragpipe_workflow) {
        error "When using --fragpipe_manifest you must also provide --fragpipe_workflow."
    }
    if (params.fragpipe_workflow && !params.fragpipe_manifest) {
        error "When using --fragpipe_workflow you must also provide --fragpipe_manifest."
    }

    manifest_ch = Channel.fromPath(params.fragpipe_manifest ?: "${params.scripts_dir}/fp.manifest.txt", checkIfExists: true)
    workflow_ch = Channel.fromPath(params.fragpipe_workflow ?: "${params.scripts_dir}/fp.dl.workflow.txt", checkIfExists: true)

    //
    // WORKFLOW: Run pipeline
    //
    OPENSEARCH (
        samplesheet,
        manifest_ch,
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
