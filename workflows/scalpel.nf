
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scalpel_pipeline'
include { ANNOTATION_PROCESSING  } from '../subworkflows/local/annotation_processing/main.nf'
include { READS_PROCESSING       } from '../subworkflows/local/reads_processing/main.nf'


workflow SCALPEL {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    genome_fasta // channel: genome FASTA file read in from --genome
    transcriptome_fasta // channel: transcriptome FASTA file read in from --transcriptome
    genome_annotation // channel: genome annotation GTF file read in from --gtf
    main:

    ch_versions = channel.empty()
    
    ANNOTATION_PROCESSING(
        ch_samplesheet,
        genome_fasta,
        transcriptome_fasta,
        genome_annotation)

    READS_PROCESSING(
        ch_samplesheet,
        ANNOTATION_PROCESSING.out.all_transcripts)

    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'scalpel_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    versions = ch_versions
}