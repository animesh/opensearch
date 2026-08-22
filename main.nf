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

    // FragPipe inputs are intentionally flexible for backwards compatibility.
    // Preferred: --fragpipe_workflow <file>.
    // Compatibility: --scripts_dir <dir> containing fp.dl.workflow.txt.
    // The manifest is generated per sample inside FRAGPIPE, so a static
    // fp.manifest.txt is not required by this pipeline.
    def workflow_path = params.fragpipe_workflow?.toString()
    if (params.scripts_dir) {
        def scripts_workflow = file("${params.scripts_dir}/fp.dl.workflow.txt")
        if (!scripts_workflow.exists()) {
            error "FragPipe workflow not found in --scripts_dir: ${scripts_workflow}"
        }
        workflow_path = scripts_workflow.toString()
    }

    if (!workflow_path) {
        workflow_path = "${projectDir}/fp.dl.workflow.txt"
    }

    def workflow_file = file(workflow_path)
    if (!workflow_file.exists()) {
        error "FragPipe workflow does not exist: ${workflow_file}"
    }

    // Keep the old CLI flags accepted without requiring a manifest.
    // FRAGPIPE creates the correct task-local four-column manifest itself.
    if (params.fragpipe_manifest) {
        def manifest_file = file(params.fragpipe_manifest.toString())
        if (!manifest_file.exists()) {
            error "FragPipe manifest does not exist: ${manifest_file}"
        }
    }

    workflow_ch = Channel.value(workflow_file)

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
