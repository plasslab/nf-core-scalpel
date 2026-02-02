#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/scalpel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/scalpel
    Website: https://nf-co.re/scalpel
    Slack  : https://nfcore.slack.com/channels/scalpel
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SCALPEL  }                from './workflows/scalpel'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_scalpel_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_scalpel_pipeline'


workflow NFCORE_SCALPEL {
    take:
    samplesheet // channel: samplesheet read in from --input
    genome_fasta // channel: genome FASTA file read in from --genome
    transcriptome_fasta // channel: transcriptome FASTA file read in from --transcriptome
    genome_annotation // channel: genome annotation GTF file read in from --gtf

    main:
    SCALPEL (
        samplesheet,
        genome_fasta,
        transcriptome_fasta,
        genome_annotation)
}

workflow {
    main:
    
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.samplesheet,
        params.help,
        params.help_full,
        params.show_hidden)


    NFCORE_SCALPEL (
        PIPELINE_INITIALISATION.out.samplesheet,
        channel.fromPath(params.genome, checkIfExists: true).collect(),
        channel.fromPath(params.transcriptome, checkIfExists: true).collect(),
        channel.fromPath(params.gtf, checkIfExists: true).collect())
    
    // PIPELINE_COMPLETION (
    //     params.email,
    //     params.email_on_fail,
    //     params.plaintext_email,
    //     params.outdir,
    //     params.monochrome_logs,
    //     params.hook_url,
    // )
}