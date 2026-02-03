
include { SALMON_INDEX }       from '../../../modules/nf-core/salmon/index/main.nf'
include { SALMON_QUANT }       from '../../../modules/nf-core/salmon/quant/main.nf'
// include { ISOFORM_WEIGHTING }  from '../../../modules/local/weighting/main.nf'
// include { SPLIT_ANNOTATION; ISOFORM_SELECTION }  from '../../../modules/local/isoformselection/main.nf'


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

    //TPM count average
    TRANSCRIPT_TPM_COUNT_AVERAGE_AND_WEIGTHING(
        SALMON_QUANT.out.results.collate(2).transpose().toList().map{ it = it[1] },
        genome_annotation)
    ch_versions = ch_versions.mix(TRANSCRIPT_TPM_COUNT_AVERAGE_AND_WEIGTHING.out.versions)


    // Isoform Selection based on 3'end region profile
    ISOFORM_SELECTION(
        TRANSCRIPT_TPM_COUNT_AVERAGE_AND_WEIGTHING.out.gtf_entries.flatten(),
        file("${projectDir}/subworkflows/local/annotation_processing/bin/transcript_collapser.R"),
        params.distance_3end,
        params.distance_profile)
    ch_versions = ch_versions.mix(ISOFORM_SELECTION.out.versions)

    // - More info about selected isoforms
    ISOFORM_SELECTION.out.stats
        .collectFile(name: 'merged_stats.txt', newLine: false).set{ statistic_annotation_processing}
    statistic_annotation_processing.view( 
        stats -> "\n-- Statistic related to annotation processing --\nseqnames\tnb.genes\tnb.isoforms\tnb.collapseds\n" + stats.text)
    
    emit:
    all_transcripts = ISOFORM_SELECTION.out.full_gtf
    all_stats = statistic_annotation_processing
    versions = ch_versions
}



// local tasks



process TRANSCRIPT_TPM_COUNT_AVERAGE_AND_WEIGTHING{
    label 'process_single'
    tag "${genome_annotation}"

    conda "conda-forge::r-tidyr=1.3.2 conda-forge::r-dplyr=1.1.4 bioconda::bioconductor-rtracklayer=1.66.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-rtracklayer_r-argparse_r-dplyr_r-tidyr:e9a06ddef3cd9da1' :
        'community.wave.seqera.io/library/bioconductor-rtracklayer_r-argparse_r-dplyr_r-tidyr:aa98afd503eec838' }"

    input:
        file quant_files
        file genome_annotation

    output:
        path('*.bed'), emit: gtf_entries
        path('versions.yml'), emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    library(dplyr)
    options(width=300)

    #reading...
    quant_data = do.call(rbind, lapply(list.files(path = ".", pattern = "*.sf", full.names = T, recursive=T), function(x) read.table(x, sep="\\t", header=T) %>% dplyr::mutate(sample_id=x)))
    gtf = data.frame(rtracklayer::import("${genome_annotation}"))[,c("type", "seqnames", "start", "end", "width", "strand", "gene_id", "gene_name", "transcript_id", "transcript_name")] %>% subset(type=="exon")
    gtf = gtf[,-1] #remove type column

    #process salmon bulk quant data...
    quant_data\$transcript_id = stringr::str_split_fixed(quant_data\$Name, pattern="\\\\|", n=2)[,1]
    quant_data = quant_data %>% 
        dplyr::distinct(sample_id,transcript_id, TPM) %>%
        dplyr::mutate(TPM = TPM + 1) %>%                    #avoid null division
        dplyr::group_by(transcript_id) %>%
        dplyr::mutate(meanTPM = mean(TPM)) %>%
        dplyr::distinct(transcript_id, meanTPM) %>%
        dplyr::ungroup()

    #merge with gtf
    gtf = dplyr::left_join(gtf, quant_data, by=c("transcript_id"))

    #calculate transcript bulk average TPM per gene
    gtf = gtf %>%
        dplyr::group_by(gene_id) %>%
        dplyr::mutate(
            TPM_perc = ifelse(
                sum(meanTPM, na.rm = TRUE) > 0,
                round(meanTPM / sum(meanTPM, na.rm = TRUE), 2),
                0
            ),
            exon_number = ifelse(strand == "+", seq(n(), 1), 1:n())
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(seqnames, start, desc(end)) %>%
        dplyr::select(-meanTPM) %>%
        dplyr::group_by(transcript_name) %>%
        dplyr::mutate(
            exon_number = ifelse(strand == "+", n():1, 1:n()),
            exon_length = end - start
        ) %>%
        dplyr::arrange(transcript_name, exon_number) %>%
        #calculate transcriptome relative coordinates
        dplyr::mutate(
            endR = cumsum(exon_length) - 1,
            startR = lag(endR + 1, default = 0)
        ) %>%
        dplyr::select(-exon_length) %>%
        dplyr::ungroup() %>% data.frame()


    #split and write bulk average TPM per gene by chromosome
    dev = lapply(unique(gtf\$seqnames), function(chr) gtf %>% subset(seqnames==chr) %>% write.table(file=paste0(chr,"_tpm.bed"), sep="\\t", row.names=F, quote=F))

    # Write versions.yml from R
    writeLines(c(
        '\"!{task.process}\":',
        paste0('    R: \"', R.version.string, '\"'),
        paste0('    tidyr: \"', packageVersion("tidyr"), '\"'),
        paste0('    dplyr: \"', packageVersion("dplyr"), '\"'),
        paste0('    rtracklayer: \"', packageVersion("rtracklayer"), '\"')
    ), "versions.yml")
    """
    stub:
    """
    touch XXX_tpm.tsv
    echo "R: stub" > versions.yml
    """
}

process ISOFORM_SELECTION{
    label 'process_single'
    tag "${gtf_entries}"

    conda "conda-forge::r-tidyr=1.3.2 conda-forge::r-dplyr=1.1.4 conda-forge::r-argparse=2.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-argparse_r-dplyr_r-tidyr:85542ba539718d05' :
        'community.wave.seqera.io/library/r-argparse_r-dplyr_r-tidyr:844cdec398eb87bd' }"

    input:
        file gtf_entries
        file script
        val distance_3end
        val distance_profile

    output:
        path "${gtf_entries.baseName}_collapsed.bed", emit: full_gtf
        path "stats.txt", emit: stats
        path('versions.yml'), emit: versions

    script:
    """
    Rscript ${script} \
        --input ${gtf_entries} \
        --output ${gtf_entries.baseName}_collapsed.bed \
        --distance_3end ${distance_3end} \
        --distance_profile ${distance_profile}
    

    # Write versions.yml for R version in bash
    echo "!${task.process}:" > versions.yml
    echo "    R: \"\$(Rscript -e 'cat(R.version.string)')\"" >> versions.yml
    """
    stub:
    """
    touch ${gtf_entries.baseName}_collapsed.bed
    touch stats.txt
    echo "R: stub" > versions.yml
    """
}