
include { SALMON_INDEX }       from '../../../modules/nf-core/salmon/index/main.nf'
include { SALMON_QUANT }       from '../../../modules/nf-core/salmon/quant/main.nf'

include { TRANSCRIPT_TPM_COUNT_AVERAGE_AND_WEIGTHING } from '../../../modules/local/transcript_tpm_count_average_and_weigthing/main'
include { ISOFORM_SELECTION                          } from '../../../modules/local/isoform_selection/main'

// main workflow
workflow ANNOTATION_PROCESSING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    genome_fasta // path: genome FASTA file read in from --genome
    transcriptome_fasta // path: transcriptome FASTA file read in from --transcriptome
    genome_annotation // path: genome annotation GTF file read in from --gtf

    main:

    ch_versions = channel.empty()

    // Indexing of transcriptome using Salmon
    SALMON_INDEX(
        genome_fasta,
        transcriptome_fasta)
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions_salmon)

    // Load fastq files (paths resolved & checked in PIPELINE_INITIALISATION)
    ch_samplesheet.map{
        sample_id, reads -> tuple(sample_id, [file(reads[0]), file(reads[1])])
    }.set{ ch_samplesheet2 }

    // Bulk Quantification of samples using Salmon
    SALMON_QUANT(
        ch_samplesheet2,
        SALMON_INDEX.out.index,
        genome_annotation,
        transcriptome_fasta,
        0,
        'ISR')
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions_salmon)

    TRANSCRIPT_TPM_COUNT_AVERAGE_AND_WEIGTHING(
        SALMON_QUANT.out.results.collate(2).transpose().toList().map{ it = it[1] },
        genome_annotation)
    ch_versions = ch_versions.mix(TRANSCRIPT_TPM_COUNT_AVERAGE_AND_WEIGTHING.out.versions)


    // Isoform Selection based on 3'end region profile
    ISOFORM_SELECTION(
        TRANSCRIPT_TPM_COUNT_AVERAGE_AND_WEIGTHING.out.gtf_entries.flatten(),
        file("${projectDir}/subworkflows/local/annotation_processing/bin/isoform_selection.R"),
        params.distance_3end,
        params.distance_profile)
    ch_versions = ch_versions.mix(ISOFORM_SELECTION.out.versions)
    
    // Collect all stats into a single log file
    ch_annotation_log = ISOFORM_SELECTION.out.stats
        .collectFile(name: 'collapsed_isoforms.log', storeDir: "${params.outdir}/pipeline_info", newLine: true)
    
    emit:
    all_transcripts = ISOFORM_SELECTION.out.full_gtf
    all_stats = ch_annotation_log
    versions = ch_versions
}
