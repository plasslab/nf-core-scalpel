
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

include { DGE_BARCODES                               } from '../../../modules/local/dge_barcodes/main'
include { SEQNAMES_CHECK                             } from '../../../modules/local/seqnames_check/main'
include { IP_SPLITTING                               } from '../../../modules/local/ip_splitting/main'
include { BED_CONVERSION                             } from '../../../modules/local/bed_conversion/main'
include { READS_MAPPING_AND_FILTERING                } from '../../../modules/local/reads_mapping_and_filtering/main'
include { IP_FILTERING                               } from '../../../modules/local/ip_filtering/main'

workflow READS_PROCESSING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_all_transcripts // path: all transcripts preprocessed from ANNOTATION_PROCESSING

    main:

    ch_versions = channel.empty()
    def barcode_tag = params.sequencing_platform == 'dropseq' ? 'XC' : 'CB'

    // Parse optional barcodes whitelist (paths relative to the whitelist csv)
    if (params.barcodes_whitelist) {
        def whitelist_dir = file(params.barcodes_whitelist).parent
        ch_barcodes_whitelist = channel
            .fromPath(params.barcodes_whitelist)
            .splitCsv(header: true)
            .map{ row ->
                def is_absolute = row.barcodes ==~ /^(\/|[a-zA-Z][a-zA-Z0-9+.-]*:\/\/).*/
                def barcodes_file = is_absolute ? file(row.barcodes) : whitelist_dir.resolve(row.barcodes)
                if( !barcodes_file.exists() ) {
                    error("Barcodes whitelist error -> sample '${row.sample}': '${row.barcodes}' not found (resolved to '${barcodes_file}'). Relative paths are resolved from the whitelist csv directory.")
                }
                tuple(row.sample, barcodes_file)
            }
    } else {
        ch_barcodes_whitelist = channel.empty()
    }

    // Extract BAM files, index & defaults Barcodes (resolved in PIPELINE_INITIALISATION)
    // files: [fastq_1, fastq_2, cranger_path, bam, bai, dge, barcodes]
    if( params.sequencing_platform == 'dropseq' ) {
        // Drop-seq: cells = barcodes of the DGE count matrix
        DGE_BARCODES(
            ch_samplesheet.map{ meta, files -> tuple(meta, file(files[5])) }
        )
        ch_versions = ch_versions.mix(DGE_BARCODES.out.versions)
        ch_samplesheet.map{
            meta, files -> tuple(meta.id, meta, files[3], files[4])
        }.join(
            DGE_BARCODES.out.barcodes.map{ meta, barcodes -> tuple(meta.id, barcodes) }
        ).map{
            sample_id, meta, bam, bai, barcodes -> tuple(meta, bam, bai, barcodes)
        }.set{ ch_bam_files }
    } else {
        // 10x Chromium: cells = filtered barcodes (Cell Ranger / STARsolo)
        ch_samplesheet.map{
            meta, files -> tuple(meta, files[3], files[4], file(files[6] ?: "${projectDir}/assets/NO_BARCODES"))
        }.set{ ch_bam_files }
    }

    // Join barcodes whitelist with BAM files if provided (override default barcodes)
    ch_bam_files.map{
        meta, bam, bai, barcodes -> tuple(meta.id, meta, bam, bai, barcodes)
    }.join(ch_barcodes_whitelist, by: 0, remainder: true).filter{
        it[1] != null
    }.map{
        sample_id, meta, bam, bai, barcodes, custom_barcodes -> tuple(meta, bam, bai, custom_barcodes ?: barcodes)
    }.branch{
        indexed: it[2]
        not_indexed: true
    }.set{ ch_bam_with_barcodes }

    // Index BAM files without index (e.g. STARsolo)
    SAMTOOLS_INDEX(
        ch_bam_with_barcodes.not_indexed.map{ meta, bam, bai, barcodes -> tuple(meta, file(bam)) }
    )
    ch_bam_with_barcodes.not_indexed.map{
        meta, bam, bai, barcodes -> tuple(meta.id, meta, bam, barcodes)
    }.join(
        SAMTOOLS_INDEX.out.index.map{ meta, bai -> tuple(meta.id, bai) }
    ).map{
        sample_id, meta, bam, barcodes, bai -> tuple(meta, bam, bai.toString(), barcodes)
    }.mix(ch_bam_with_barcodes.indexed).set{ ch_indexed_bams }

    // Check chromosome names (BAM vs annotation) & barcode tags before processing (#20, #26, #36)
    SEQNAMES_CHECK(
        ch_indexed_bams,
        ch_all_transcripts.map{ it = it[0] }.collect(),
        barcode_tag
    )
    ch_versions = ch_versions.mix(SEQNAMES_CHECK.out.versions)

    // Convert BAM to BED files
    BED_CONVERSION(
        SEQNAMES_CHECK.out.checked_bams.combine(ch_all_transcripts.map{ it=it[0] }),
        params.sequencing_platform
    )
    ch_versions = ch_versions.mix(BED_CONVERSION.out.versions)
    

    // Map reads in BED files to transcripts from annotation
    ch_all_transcripts.combine(
        BED_CONVERSION.out.bed_files, by: 0
    ).map{
        it = tuple(it[2], it[0], it[1], it[3])
    }.set{ reads_and_annots }


    READS_MAPPING_AND_FILTERING(
        reads_and_annots,
        file("${projectDir}/subworkflows/local/reads_processing/bin/reads_mapping_and_filtering.R"),
        params.distance_profile
    )
    ch_versions = ch_versions.mix(READS_MAPPING_AND_FILTERING.out.versions)


    // Split internal priming reference by chromosome (once, instead of per task)
    IP_SPLITTING(
        file(params.ip_reference, checkIfExists: true),
        ch_all_transcripts.map{ it = it[0] }.collect()
    )
    ch_versions = ch_versions.mix(IP_SPLITTING.out.versions)
    IP_SPLITTING.out.ip_chrs.flatten().map{
        it = tuple(it.name - '.ip.txt', it)
    }.set{ ch_ip_chrs }

    // Filter reads associated to internal priming positions
    READS_MAPPING_AND_FILTERING.out.mapped_reads.map{
        meta, seqname, mapped_reads -> tuple(seqname, meta, mapped_reads)
    }.combine(ch_ip_chrs, by: 0).map{
        seqname, meta, mapped_reads, ip_chr -> tuple(meta, seqname, mapped_reads, ip_chr)
    }.set{ reads_and_ips }

    IP_FILTERING(
        reads_and_ips,
        file("${projectDir}/subworkflows/local/reads_processing/bin/ip_filtering.R"),
        params.distance_ip
    )
    ch_versions = ch_versions.mix(IP_FILTERING.out.versions)


    emit:
    reads = IP_FILTERING.out.filtered_mapped_reads
    chr_bams = BED_CONVERSION.out.bed_files     // channel: [ val(seqname), val(meta), path(chr.bam) ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
