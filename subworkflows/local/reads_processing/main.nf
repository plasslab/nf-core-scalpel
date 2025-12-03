
include { BED_CONVERSION       } from '../../../modules/local/bedconversion/main.nf'
include { READ_MAPPING        }  from '../../../modules/local/readmapping/main.nf'

workflow READS_PROCESSING {

    take:
    ch_samples // channel: [ val(meta), [ bam ] ]
    ch_isoform_gtf // channel: path to isoform GTF file

    main:

    ch_versions = Channel.empty()

    // Extract BAM files & defaults Barcodes tag
    if( params.sequencing_platform == 'chromium' ) {
        ch_samples.map{ meta, files -> tuple(meta, [files[2]+'/outs/possorted_genome_bam.bam', 
                files[2]+'/outs/filtered_feature_bc_matrix/barcodes.tsv.gz']) }.set{ ch_bam_files }
    } else if( params.sequencing_platform == 'dropseq' ) {
        ch_samples.map{ meta, files -> tuple(meta, files[2]) }.set{ ch_bam_files }
    } else {
        log.error "Unsupported sequencing platform: ${params.sequencing_platform}. Supported platforms are: 'chromium' and 'dropseq'."
        System.exit(1)
    }
    
    // Convert BAM to BED files
    BED_CONVERSION(
        ch_bam_files)
    ch_versions = ch_versions.mix(BED_CONVERSION.out.versions)

    BED_CONVERSION.out.bed_files.combine(ch_isoform_gtf).map{ meta, bed_file, chr_id, isoform_gtf -> 
        tuple([meta.id, chr_id], bed_file, isoform_gtf) }.set { ch_bed_and_gtf_files }

    // READ_MAPPING(
    //     ch_bed_and_gtf_files)


    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}
