

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROBABILITY_DISTRIBUTION                   } from '../../../modules/local/probability_distribution/main'
include { FRAGMENT_PROBABILITY                       } from '../../../modules/local/fragment_probability/main'
include { CELL_DATA_GENERATION                       } from '../../../modules/local/cell_data_generation/main'
include { EM_ALGORITHM                               } from '../../../modules/local/em_algorithm/main'
include { DGE_GENERATION                             } from '../../../modules/local/dge_generation/main'

workflow ISOFORM_QUANTIFICATION {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    mapped_reads // channel: mapped reads from READS_PROCESSING

    main:
    ch_versions = channel.empty()

    mapped_reads.map{ it = tuple(it[0], it[2]) }.set{ mapped_reads}

    // Calculate probability distribution for each sample on transcriptomic space
    PROBABILITY_DISTRIBUTION(
        mapped_reads,
        file("${workflow.projectDir}/subworkflows/local/isoform_quantification/bin/probability_distribution.R"),
        params.gene_fraction,
        params.binsize
    )
    ch_versions = ch_versions.mix(PROBABILITY_DISTRIBUTION.out.versions)
    PROBABILITY_DISTRIBUTION.out.probabilities.join(mapped_reads).map{
        it = [it[0], it[1], it[3]]
    }.transpose().set{ reads_probs }

    // Calculate fragment probablity for each read and isoform combination
    FRAGMENT_PROBABILITY(
        reads_probs,
        file("${workflow.projectDir}/subworkflows/local/isoform_quantification/bin/fragment_probability.R")
    )
    ch_versions = ch_versions.mix(FRAGMENT_PROBABILITY.out.versions)

    // Generate cell-level data for EM algorithm
    CELL_DATA_GENERATION(FRAGMENT_PROBABILITY.out.isoform_probabilities.groupTuple())
    ch_versions = ch_versions.mix(CELL_DATA_GENERATION.out.versions)

    // Create channel of sample_id, cell_chunk
    CELL_DATA_GENERATION.out.cell_isoform_probabilities.transpose().set{ cell_chunks }

    // Perform EM algorithm for each sample and chunk of cells
    EM_ALGORITHM(
        cell_chunks,
        file("${workflow.projectDir}/subworkflows/local/isoform_quantification/bin/em_algorithm.R")
    )
    ch_versions = ch_versions.mix(EM_ALGORITHM.out.versions)

    // Extract DGE count matrix for each sample
    // files: [fastq_1, fastq_2, cranger_path, bam, bai, dge_dir, barcodes]
    ch_samplesheet.map{ meta, files -> tuple(meta.id, files[5]) }.join(
        EM_ALGORITHM.out.em_results.groupTuple(by: 0)
    ).map{
        sample_id, sample_dge_path, em_results -> tuple(sample_id, file(sample_dge_path), em_results)
    }.set{ dge_inputs }

    // Generate isoform DGE matrix (sparse, 10x format) for each sample
    DGE_GENERATION(
        dge_inputs,
        file("${workflow.projectDir}/subworkflows/local/isoform_quantification/bin/dge_generation.R")
    )
    ch_versions = ch_versions.mix(DGE_GENERATION.out.versions)

    emit:
    idge = DGE_GENERATION.out.dge_matrices    // channel: [ val(meta), path(iDGE dir) ]
    versions = ch_versions                     // channel: [ versions.yml ]

}
