
include { SALMON_INDEX }       from '../../../modules/nf-core/salmon/index/main.nf'
include { SALMON_QUANT }       from '../../../modules/nf-core/salmon/quant/main.nf'

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
    
    // Collect all stats into a single log file
    ch_annotation_log = ISOFORM_SELECTION.out.stats
        .collectFile(name: 'collapsed_isoforms.log', storeDir: "${params.outdir}/pipeline_info", newLine: true)
    
    emit:
    all_transcripts = ISOFORM_SELECTION.out.full_gtf
    all_stats = ch_annotation_log
    versions = ch_versions
}



// local tasks //

process TRANSCRIPT_TPM_COUNT_AVERAGE_AND_WEIGTHING{
    label 'process_single'
    tag "${genome_annotation}"

    conda "bioconda::bioconductor-rtracklayer=1.66.0 conda-forge::r-dplyr=1.2.0 conda-forge::r-tidyr=1.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-rtracklayer_r-data.table_r-dplyr_r-tidyr:95f9946a88b84d7a' :
        'community.wave.seqera.io/library/bioconductor-rtracklayer_r-data.table_r-dplyr_r-tidyr:60d194d659530407' }"

    input:
        file quant_files
        path genome_annotation

    output:
        path('*.tsv'), emit: gtf_entries
        path('versions.yml'), emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    library(dplyr)
    options(width=300)

    # - reading...
    sf_files = list.files(
        path = ".", 
        pattern = "quant.sf", 
        full.names = T, 
        recursive=T)
    
    if(length(sf_files) == 0){
        stop("No quant.sf files found in input directory [Issue with Salmon QUANT ? / input FASTA ?]")
    }

    quant_data = do.call(
        rbind, 
        lapply(
            sf_files,
            function(x) read.table(x, sep="\\t", header=T) %>% dplyr::mutate(sample_id=x)))
    
    gtf = data.frame(
        rtracklayer::import("${genome_annotation}"))[,
            c("type", "seqnames", "start", "end", "width", 
                "strand", "gene_id", "gene_name", "transcript_id", "transcript_name")] %>%
            subset(type=="exon")
    gtf = gtf[,-1] #remove type column

    # - process salmon bulk quant data...
    quant_data\$transcript_id = stringr::str_split_fixed(quant_data\$Name, pattern="\\\\|", n=2)[,1]
    quant_data = quant_data %>% 
        dplyr::distinct(sample_id,transcript_id, TPM) %>%
        dplyr::mutate(TPM = TPM + 1) %>%                    #avoid null division
        dplyr::group_by(transcript_id) %>%
        dplyr::mutate(meanTPM = mean(TPM)) %>%
        dplyr::distinct(transcript_id, meanTPM) %>%
        dplyr::ungroup()

    # - merge with gtf...
    gtf = dplyr::left_join(gtf, quant_data, by=c("transcript_id")) %>%
        # if transcriptID is not in quant_data, assign a very small constant to avoid 0 weigth score (meanTPM = 1e-12)
        dplyr::mutate(meanTPM = ifelse(is.na(meanTPM), 1e-12, meanTPM))

    #calculate transcript bulk average TPM per gene
    gtf = gtf %>%
        dplyr::group_by(gene_id) %>%
        dplyr::mutate(
            TPM_perc = round(meanTPM / sum(meanTPM), 12),
            width = width - 1
        ) %>%
        dplyr::ungroup() %>%
        arrange(seqnames,start,desc(end)) %>% 
        dplyr::group_by(transcript_name) %>%
        dplyr::mutate(exon_number = ifelse(strand == "+", seq(n(), 1), 1:n())) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(seqnames, start, desc(end)) %>%
        dplyr::select(-meanTPM) %>%
        dplyr::arrange(transcript_name,exon_number) %>% 
        dplyr::group_by(transcript_name) %>%
        group_modify(
            ~ {
                starts = .x\$start
                ends = .x\$end
                startR = rep(NA,length(starts))
                endR = rep(NA,length(ends))
                for(i in .x\$exon_number){
                    if(i==1){
                        startR[i] = 0
                        endR[i] = (ends[i] - (starts[i])) - 1
                    } else {
                        startR[i] = endR[i-1] + 1
                        endR[i] = startR[i] + (ends[i] - starts[i])
                    }
                    .x\$startR = startR
                    .x\$endR = endR
                }
                .x
            }
        ) %>%
        distinct(
            seqnames,
            start,
            end,
            width,
            strand,
            startR,
            endR,
            gene_id,
            gene_name,
            transcript_id,
            transcript_name,
            exon_number,
            TPM_perc
        ) %>%
        dplyr::ungroup() %>% data.frame()


    # - split and write bulk average TPM per gene by chromosome
    dev = lapply(
        unique(gtf\$seqnames), 
        function(chr) gtf %>% subset(seqnames==chr) %>%
             write.table(file=paste0(chr,".tsv"), 
             sep="\\t", row.names=F, quote=F))

    # - Write versions.yml from R
    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"'),
        paste0('    tidyr: "', packageVersion("tidyr"), '"'),
        paste0('    dplyr: "', packageVersion("dplyr"), '"'),
        paste0('    rtracklayer: "', packageVersion("rtracklayer"), '"')
    ), "versions.yml")
    """
    stub:
    """
    touch XXX.tsv
    echo "R: stub" > versions.yml
    """
}

process ISOFORM_SELECTION{
    label 'process_single'
    tag "${gtf_entries.baseName}"

    conda "conda-forge::r-argparse=2.3.1 conda-forge::r-tidyr=1.3.2 conda-forge::r-dplyr=1.1.4 conda-forge::r-data.table=1.17.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-argparse_r-data.table_r-dplyr_r-tidyr:130e888376e46c19' :
        'community.wave.seqera.io/library/r-argparse_r-data.table_r-dplyr_r-tidyr:40ad9129fa7b497b' }"

    input:
        file gtf_entries
        file script
        val distance_3end
        val distance_profile

    output:
        tuple val("${gtf_entries.baseName}"), path("${gtf_entries.baseName}.annots.gz"), emit: full_gtf
        path "stats.txt", emit: stats
        path('versions.yml'), emit: versions

    script:
    """
    Rscript ${script} \
        --input ${gtf_entries} \
        --output ${gtf_entries.baseName}.annots \
        --distance_3end ${distance_3end} \
        --distance_profile ${distance_profile}
    gzip ${gtf_entries.baseName}.annots


    # Write versions.yml for R version in bash
    echo "!${task.process}:" > versions.yml
    echo "    R: \"\$(Rscript -e 'cat(R.version.string)')\"" >> versions.yml
    """
    stub:
    """
    touch ${gtf_entries.baseName}.annots.gz
    touch stats.txt
    echo "R: stub" > versions.yml
    """
}