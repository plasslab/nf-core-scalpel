process CELL_DATA_GENERATION {
    label 'process_single'
    tag "${meta}"

    conda "conda-forge::r-base=4.4.3 conda-forge::r-data.table=1.17.8 conda-forge::r-dplyr=1.1.4 conda-forge::r-fst=0.9.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-fst:cc9686c619b6bf7a' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-fst:7d28f8f9dc2b7888' }"

    input:
        tuple val(meta), file(isoform_probabilities)

    output:
        tuple val("${meta}"), path("*.fst"), emit: cell_isoform_probabilities
        path('versions.yml'), emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
    suppressPackageStartupMessages(library(fst, quietly = TRUE))
    suppressPackageStartupMessages(library(data.table, quietly = TRUE))
    
    all_reads = lapply(
        list.files("./", pattern = "*.frag.fst"),
        function(f) {
            df = read_fst(f, as.data.table = TRUE)
        }
    ) %>% rbindlist() %>% distinct()

    #generate chunks of 500 bcs
    all_barcodes = unique(all_reads\$bc)
    bc_chunks = split(all_barcodes, ceiling(seq_along(all_barcodes)/500))

    #write cells files for each chunks
    lapply(
        names(bc_chunks),
        function(i) {
            chunk_bcs = bc_chunks[[i]]
            chunk_reads = all_reads %>% subset(bc %in% chunk_bcs)
            write_fst(chunk_reads, paste0("${meta}_","cell_chunk_", i, ".fst"))
        }
    )
    # - Write versions.yml from R
    writeLines(
        c(
            '"${task.process}":',
            paste0('    R: "', R.version.string, '"'),
            paste0('    dplyr: "', packageVersion("dplyr"), '"')
        ), "versions.yml")
    """
}
