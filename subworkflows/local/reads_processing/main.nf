

workflow READS_PROCESSING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_all_transcripts // path: all transcripts preprocessed from ANNOTATION_PROCESSING

    main:

    ch_versions = channel.empty()

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
    // ch_bam_files.view()
    
    // Convert BAM to BED files
    BED_CONVERSION(ch_bam_files)
    ch_versions = ch_versions.mix(BED_CONVERSION.out.versions)

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
    tuple val(meta), file(sample)

    output:
    path("*.bed"), emit: bed_files
    path "versions.yml",   emit: versions

    script:
    def barcode_file = sample[1].toString().endsWith('.gz') ? "<(zcat ${sample[1]})" : "${sample[1]}"
    """
        # Convert BAM to BED
        samtools index ${sample[0]}
        samtools view -b ${sample[0]} -X ${sample[0]}.bai -D CB:${barcode_file} --keep-tag "CB,UB" > tmp.bam
        bam2bed --all-reads --split < tmp.bam | 
        gawk -v OFS="\\t" '{print \$1,\$2,\$3,\$6,\$4,\$14"::"\$15}' > ${meta.id}_${sample[0].baseName}.bed

        #split by chromosome
        awk -v OFS="\t" '{print > ${meta.id}_\$1".bed"}' ${meta.id}_${sample[0].baseName}.bed
        rm tmp.bam
        rm ${meta.id}_${sample[0].baseName}.bed

        # Write versions.yml for bedtools version in bash
        echo "!${task.process}:" > versions.yml
        samtools --version >> versions.yml
        gawk --version >> versions.yml
    """
}