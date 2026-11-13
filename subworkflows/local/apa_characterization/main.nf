
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DIFFERENTIAL_ISOFORM_USAGE                 } from '../../../modules/local/differential_isoform_usage/main'
include { FILTERED_BAMS                              } from '../../../modules/local/filtered_bams/main'

workflow APA_CHARACTERIZATION {

    take:
    ch_idge // channel: [ val(meta), path(iDGE dir) ] from ISOFORM_QUANTIFICATION
    ch_filtered_reads // channel: [ val(meta), val(seqname), path(ip_filtered.fst) ] from READS_PROCESSING
    ch_chr_bams // channel: [ val(seqname), val(meta), path(chr.bam) ] from READS_PROCESSING
    ch_annots // channel: [ val(seqname), path(annots.gz) ] selected isoforms (3' ends)
    ch_probabilities // channel: [ val(meta), path(prob), path(pdf) ] read distance distributions (library QC)

    main:
    ch_versions = channel.empty()

    // Optional clusters annotation (orig.ident, Barcodes, Clusters)
    def clusters_file = params.clusters ? file(params.clusters, checkIfExists: true) : file("${projectDir}/assets/NO_CLUSTERS")

    // Differential isoform usage between clusters (or samples if no clusters provided)
    DIFFERENTIAL_ISOFORM_USAGE(
        ch_idge.map{ it = it[1] }.collect(),
        clusters_file,
        ch_annots.map{ it = it[1] }.collect(),
        ch_probabilities.map{ it = it[1] }.collect(),
        params.diu_test,
        params.diu_min_delta,
        params.diu_permutations,
        file("${projectDir}/subworkflows/local/apa_characterization/bin/differential_isoform_usage.R"),
        file("${projectDir}/subworkflows/local/apa_characterization/bin/differential_isoform_usage_functions.R")
    )
    ch_versions = ch_versions.mix(DIFFERENTIAL_ISOFORM_USAGE.out.versions)
    // Warn on sample comparisons without replicates / library differences
    DIFFERENTIAL_ISOFORM_USAGE.out.warnings.subscribe{ it.readLines().each{ line -> log.warn(line) } }

    // BAM file of the reads kept by SCALPEL for each sample (e.g. coverage visualization)
    ch_chr_bams.map{
        seqname, meta, bam -> tuple(meta, bam)
    }.groupTuple(by: 0).join(
        ch_filtered_reads.map{ meta, seqname, reads -> tuple(meta, reads) }.groupTuple(by: 0)
    ).set{ ch_bams_and_reads }
    FILTERED_BAMS(ch_bams_and_reads)
    ch_versions = ch_versions.mix(FILTERED_BAMS.out.versions)

    emit:
    diu_table = DIFFERENTIAL_ISOFORM_USAGE.out.diu_table   // path: DIU_table.csv
    diu_genes = DIFFERENTIAL_ISOFORM_USAGE.out.diu_genes   // path: DIU_genes.csv
    seurat_obj = DIFFERENTIAL_ISOFORM_USAGE.out.seurat_obj // path: iDGE_seurat.RDS
    filtered_bams = FILTERED_BAMS.out.filtered_bams        // channel: [ val(meta), bam, bai ]
    versions = ch_versions                                 // channel: [ versions.yml ]
}
