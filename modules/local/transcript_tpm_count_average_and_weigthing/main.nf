process TRANSCRIPT_TPM_COUNT_AVERAGE_AND_WEIGTHING{
    label 'process_single'
    tag "${genome_annotation}"

    conda "conda-forge::r-base=4.4.3 bioconda::bioconductor-rtracklayer=1.66.0 conda-forge::r-dplyr=1.2.0 conda-forge::r-tidyr=1.3.2 conda-forge::r-stringr=1.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_bioconductor-rtracklayer_r-dplyr_r-tidyr_r-stringr:4632629d2ddfef49' :
        'community.wave.seqera.io/library/r-base_bioconductor-rtracklayer_r-dplyr_r-tidyr_r-stringr:68563737921154c4' }"

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

    #calculate transcript bulk average TPM per gene (one value per transcript, as original quantification_processing.R;
    #summing over the exon rows weighted each transcript by its number of exons)
    tpm_perc = gtf %>%
        dplyr::distinct(gene_id, transcript_id, meanTPM) %>%
        dplyr::group_by(gene_id) %>%
        dplyr::mutate(TPM_perc = round(meanTPM / sum(meanTPM), 12)) %>%
        dplyr::ungroup() %>%
        dplyr::select(transcript_id, TPM_perc)
    gtf = gtf %>%
        dplyr::left_join(tpm_perc, by = "transcript_id") %>%
        dplyr::mutate(width = width - 1) %>%
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
