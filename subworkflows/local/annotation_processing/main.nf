
include { SALMON_INDEX }       from '../../../modules/nf-core/salmon/index/main.nf'
include { SALMON_QUANT }       from '../../../modules/nf-core/salmon/quant/main.nf'
include { ISOFORM_WEIGHTING }  from '../../../modules/local/weighting/main.nf'
include { SPLIT_ANNOTATION; ISOFORM_SELECTION }  from '../../../modules/local/isoformselection/main.nf'

workflow ANNOTATION_PROCESSING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    genome_fasta // path: genome FASTA file read in from --genome
    transcriptome_fasta // path: transcriptome FASTA file read in from --transcriptome
    genome_annotation // path: genome annotation GTF file read in from --gtf

    main:

    ch_versions = Channel.empty()

    // Indexing of transcriptome using Salmon
    SALMON_INDEX(
        genome_fasta,
        transcriptome_fasta)
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

    // Bulk Quantification of samples using Salmon
    SALMON_QUANT(
        ch_samplesheet.map{ sample_id, reads -> tuple( sample_id, [reads[0], reads[1]])},
        SALMON_INDEX.out.index,
        genome_annotation,
        transcriptome_fasta,
        0,
        'ISR')
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

    // Splitting genome annotation GTF by chromosome
    SPLIT_ANNOTATION(genome_annotation)
    ch_versions = ch_versions.mix(SPLIT_ANNOTATION.out.versions)

    // Weighting of Isoform based on bulk quantification results
    ISOFORM_WEIGHTING(
        SALMON_QUANT.out.results.collate(2).transpose().toList().map{ it = it[1] })
    ch_versions = ch_versions.mix(ISOFORM_WEIGHTING.out.versions)

    // Isoform Selection based on 3'end region profile
    ISOFORM_SELECTION(
        SPLIT_ANNOTATION.out.chromosome_gtfs.flatten(),
        ISOFORM_WEIGHTING.out.weights)
    ch_versions = ch_versions.mix(ISOFORM_SELECTION.out.versions)

    // emit:
    // all_isoforms = ISOFORM_SELECTION.out.full_gtf
    // unique_isoforms = ISOFORM_SELECTION.out.unique_gtf
    // versions = ch_versions
}
