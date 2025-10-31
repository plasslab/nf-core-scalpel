// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SALMON_INDEX }       from '../../../modules/nf-core/salmon/index/main.nf'
include { SALMON_QUANT }       from '../../../modules/nf-core/salmon/quant/main.nf'

workflow ANNOTATION_PROCESSING {

    take:
    // TODO nf-core: edit input (take) channels
    genome_fasta // channel: genome FASTA file read in from --genome
    transcriptome_fasta // channel: transcriptome FASTA file read in from --transcriptome

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    // Indexing of transcriptome using Salmon
    SALMON_INDEX(
        genome_fasta,
        transcriptome_fasta)

    
    emit:
    // TODO nf-core: edit emitted channels
    versions = ch_versions                     // channel: [ versions.yml ]
}
