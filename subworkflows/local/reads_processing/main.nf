
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
        ch_barcodes_whitelist = channel
            .fromPath(params.barcodes_whitelist)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample, file(row.barcodes)) }
    } else {
        ch_barcodes_whitelist = channel.empty()
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

        log.error 
        """
            Unsupported sequencing platform: ${params.sequencing_platform}. 
            Supported platforms are: 'chromium' and 'dropseq'.
        """
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
        ch_bam_with_barcodes.combine(ch_all_transcripts.map{ it=it[0] }),
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
        file("${projectDir}/subworkflows/local/reads_processing/bin/reads_merging_v3.R"),
        params.distance_profile
    )
    ch_versions = ch_versions.mix(READS_MAPPING_AND_FILTERING.out.versions)


    // Filter reads associated to internal priming positions
    IP_FILTERING(
        READS_MAPPING_AND_FILTERING.out.mapped_reads,
        params.ip_reference,
        file("${projectDir}/subworkflows/local/reads_processing/bin/ip_filtering.R"),
        params.distance_ip
    )
    ch_versions = ch_versions.mix(IP_FILTERING.out.versions)


    emit:
    reads = IP_FILTERING.out.filtered_mapped_reads
    versions = ch_versions                     // channel: [ versions.yml ]
}



process BED_CONVERSION {
    label 'process_single'
    tag "${meta.id}, ${annotated_seqname}"

    conda "bioconda::samtools=1.23 bioconda::bedops=2.4.42 conda-forge::gawk=5.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/35c889f53d969f55720d9292b6418744f68b405781c91e674ba98d5e8802a4c1/data' :
        'community.wave.seqera.io/library/bedops_samtools_gawk:2ceaacb8a83d256c' }"

    input:
        tuple val(meta), val(bam_file), path(barcodes_file), val(annotated_seqname)
        val sequencing_platform

    output:
        tuple val("${annotated_seqname}"), val("${meta.id}"), path("${annotated_seqname}.bam"), emit: bed_files, optional: true
        path('versions.yml'), emit: versions

    script:
    def barcode_tag = sequencing_platform == 'dropseq' ? 'XC' : 'CB'
    def umi_tag = sequencing_platform == 'dropseq' ? 'XM' : 'UB'
    def has_barcodes = barcodes_file.name != 'NO_BARCODES'
    def barcode_input = has_barcodes ? (barcodes_file.toString().endsWith('.gz') ? "<(gunzip -c ${barcodes_file})" : "${barcodes_file}") : ""
    def filter_cmd = has_barcodes ? "-D ${barcode_tag}:${barcode_input}" : ""
    """
    # Convert BAM to BED
    samtools view \
        -b ${bam_file} \
        -X ${bam_file}.bai ${filter_cmd} \
        --keep-tag "${barcode_tag},${umi_tag}" ${annotated_seqname} > ${annotated_seqname}.bam

    # Write versions.yml for bedtools version in bash
    echo "!${task.process}:" > versions.yml
    samtools --version >> versions.yml
    gawk --version >> versions.yml
    """
    stub:
    """
    touch XXX.bam
    echo "R: stub" > versions.yml
    """
}


process READS_MAPPING_AND_FILTERING {
    label 'process_single'
    tag "${meta}, ${seqname}"

    conda "conda-forge::r-argparse=2.3.1 conda-forge::r-dplyr=1.2.0 conda-forge::r-data.table=1.17.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bedops_r-argparse_r-data.table_r-dplyr_pruned:802521a287a86b2e' :
        'community.wave.seqera.io/library/bedops_r-argparse_r-data.table_r-dplyr_pruned:688bbcb08969885b' }"

    input:
        tuple val(meta), val(seqname), path(annots), val(reads)
        path(script)
        val(distance_thr)
    
    output:
        tuple val("${meta}"), val("${seqname}"), path("${meta}_${reads.baseName}.fst"), emit: mapped_reads
        path('versions.yml'), emit: versions

    script:
    """
    Rscript ${script} \
        --reads 'bam2bed --all-reads --split < ${reads}' \
        --annotations ${annots} \
        --output ${meta}_${reads.baseName}.fst \
        --distance_threshold ${distance_thr}
    
    # Write versions.yml for bedtools version in bash
    echo "!${task.process}:" > versions.yml
    """
    stub:
    """
    touch XXX.fst
    echo "R: stub" > versions.yml
    """
}


process IP_FILTERING {
    label 'process_single'
    tag "${meta}, ${seqname}"

    conda "conda-forge::r-argparse=2.3.1 conda-forge::r-dplyr=1.2.0 conda-forge::r-data.table=1.17.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-genomicranges_r-argparse_r-data.table_r-dplyr_pruned:abda756df030d27e' :
        'community.wave.seqera.io/library/bioconductor-genomicranges_r-argparse_r-data.table_r-dplyr_r-stringr:750d9136cb69c910' }"

    input:
        tuple val(meta), val(seqname), val(mapped_reads)
        path(ip_reference)
        path(script)
        val(distance_threshold)

    output:
        tuple val("${meta}"), val("${seqname}"), path("${mapped_reads.baseName}.ip_filtered.fst"), emit: filtered_mapped_reads
        path('versions.yml'), emit: versions

    script:
    """
    Rscript ${script} \
        --reads ${mapped_reads} \
        --ip_reference ${ip_reference} \
        --output ${mapped_reads.baseName}.ip_filtered.fst \
        --threshold_dist ${distance_threshold} \
        --threads ${task.cpus}
    
    # Write versions.yml for bedtools version in bash
    echo "!${task.process}:" > versions.yml
    """
    stub:
    """
    touch XXX.ip_filtered.fst
    echo "R: stub" > versions.yml
    """
}


