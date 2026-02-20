
// include { BED_CONVERSION       } from '../../../modules/local/bedconversion/main.nf'
// include { READ_MAPPING        }  from '../../../modules/local/readmapping/main.nf'

workflow READS_PROCESSING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_all_transcripts // path: all transcripts preprocessed from ANNOTATION_PROCESSING

    main:

    ch_versions = channel.empty()

    // Parse optional barcodes whitelist
    if (params.barcodes_whitelist) {
        ch_barcodes_whitelist = Channel
            .fromPath(params.barcodes_whitelist)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample, file(row.barcodes)) }
    } else {
        ch_barcodes_whitelist = Channel.empty()
    }

    // Extract BAM files & defaults Barcodes tag
    if( params.sequencing_platform == 'chromium' ) {

        ch_samplesheet.map{ 
            meta, 
            files -> tuple(meta, [files[2]+'/outs/possorted_genome_bam.bam', 
            files[2]+'/outs/filtered_feature_bc_matrix/barcodes.tsv.gz']) 
        }.set{ ch_bam_files }

    } else if( params.sequencing_platform == 'dropseq' ) {

        ch_samplesheet.map{ meta, files -> tuple(meta, files[2]) }.set{ ch_bam_files }

    } else {

        log.error "Unsupported sequencing platform: ${params.sequencing_platform}. Supported platforms are: 'chromium' and 'dropseq'."
        System.exit(1)

    }

    // Join barcodes whitelist with BAM files if provided (override default barcodes)
    if (params.barcodes_whitelist) {
        ch_bam_with_barcodes = ch_bam_files
            .map { meta, files -> tuple(meta.id, meta, files) }
            .join(ch_barcodes_whitelist, by: 0, remainder: true)
            .map { sample_id, meta, files, custom_barcodes ->
                def barcodes = custom_barcodes ?: (files instanceof List ? file(files[1]) : file("NO_BARCODES"))
                def bam = files instanceof List ? file(files[0]) : file(files)
                tuple(meta, bam, barcodes)
            }
    } else {
        ch_bam_with_barcodes = ch_bam_files.map { meta, files ->
            def barcodes = files instanceof List ? file(files[1]) : file("NO_BARCODES")
            def bam = files instanceof List ? file(files[0]) : file(files)
            tuple(meta, bam, barcodes)
        }
    }
    
    // Convert BAM to BED files
    BED_CONVERSION(
        ch_bam_with_barcodes,
        params.sequencing_platform)
    // ch_versions = ch_versions.mix(BED_CONVERSION.out.versions)

    // BED_CONVERSION.out.bed_files.combine(ch_all_transcripts).map{ meta, bed_file, chr_id, all_transcripts -> 
    //     tuple([meta.id, chr_id], bed_file, all_transcripts) }.set { ch_bed_and_gtf_files }

    // READ_MAPPING(
    //     ch_bed_and_gtf_files)

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}



process BED_CONVERSION {
    label 'process_single'
    tag "${meta.id}"

    conda "bioconda::samtools=1.23 bioconda::bedops=2.4.42 conda-forge::gawk=5.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bedops_samtools_gawk:bcf972bff00da13a' :
        'community.wave.seqera.io/library/bedops_samtools_gawk:00501a873fceef75' }"

    input:
    tuple val(meta), path(bam_file), path(barcodes_file)
    val sequencing_platform

    output:
    path("*.bed"), emit: bed_files
    path "versions.yml"

    script:
    def barcode_tag = sequencing_platform == 'dropseq' ? 'XC' : 'CB'
    def umi_tag = sequencing_platform == 'dropseq' ? 'XM' : 'UB'
    def has_barcodes = barcodes_file.name != 'NO_BARCODES'
    def barcode_input = has_barcodes ? (barcodes_file.toString().endsWith('.gz') ? "<(gunzip -c ${barcodes_file})" : "${barcodes_file}") : ""
    def filter_cmd = has_barcodes ? "-D ${barcode_tag}:${barcode_input}" : ""
    """
        # Convert BAM to BED
        samtools index ${bam_file}
        samtools view -b ${bam_file} -X ${bam_file}.bai ${filter_cmd} --keep-tag "${barcode_tag},${umi_tag}" > tmp.bam
        bam2bed --all-reads --split < tmp.bam | 
        gawk -v OFS="\\t" '{print \$1,\$2,\$3,\$6,\$4,\$14"::"\$15}' > ${meta.id}_${bam_file.baseName}.bed

        #split by chromosome
        awk -v OFS="\\t" '{print > ("${meta.id}_" \$1 ".bed")}' ${meta.id}_${bam_file.baseName}.bed
        rm tmp.bam
        rm ${meta.id}_${bam_file.baseName}.bed

        # Write versions.yml for bedtools version in bash
        echo "!${task.process}:" > versions.yml
        samtools --version >> versions.yml
        gawk --version >> versions.yml
    """
}