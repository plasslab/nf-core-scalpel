

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ISOFORM_QUANTIFICATION {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    mapped_reads // channel: mapped reads from READS_PROCESSING

    main:
    ch_versions = Channel.empty()

    mapped_reads.map{ it = tuple(it[0], it[2]) }.set{ mapped_reads}

    // Calculate probability distribution for each sample on transcriptomic space
    PROBABILITY_DISTRIBUTION(
        mapped_reads,
        file("${workflow.projectDir}/subworkflows/local/isoform_quantification/bin/compute_prob.R"),
        params.gene_fraction,
        params.binsize
    )
    ch_versions = ch_versions.mix(PROBABILITY_DISTRIBUTION.out.versions)
    PROBABILITY_DISTRIBUTION.out.probabilities.join(mapped_reads).map{
        it = [it[0], it[1], it[3]]
    }.transpose().set{ reads_probs }

    // Calculate fragment probablity for each read and isoform combination
    FRAGMENT_PROBABILITY(
        reads_probs,
        file("${workflow.projectDir}/subworkflows/local/isoform_quantification/bin/frag_prob.R")
    )
    ch_versions = ch_versions.mix(FRAGMENT_PROBABILITY.out.versions)

    // Generate cell-level data for EM algorithm
    CELL_DATA_GENERATION(FRAGMENT_PROBABILITY.out.isoform_probabilities.groupTuple())
    ch_versions = ch_versions.mix(CELL_DATA_GENERATION.out.versions)

    // Create channel of sample_id, cell_chunk
    CELL_DATA_GENERATION.out.cell_isoform_probabilities.transpose().set{ cell_chunks }

    // Perform EM algorithm for each sample and chunk of cells
    EM_ALGORITHM(
        cell_chunks,
        file("${workflow.projectDir}/subworkflows/local/isoform_quantification/bin/em_algorithm.R")
    )
    ch_versions = ch_versions.mix(EM_ALGORITHM.out.versions)

    // Extract DGE count matrix for each sample
    // files: [fastq_1, fastq_2, cranger_path, bam, bai, dge_dir, barcodes]
    ch_samplesheet.map{ meta, files -> tuple(meta.id, files[5]) }.join(
        EM_ALGORITHM.out.em_results.groupTuple(by: 0)
    ).map{
        sample_id, sample_dge_path, em_results -> tuple(sample_id, file(sample_dge_path), em_results)
    }.set{ dge_inputs }

    // Generate isoform DGE matrix (sparse, 10x format) for each sample
    DGE_GENERATION(
        dge_inputs,
        file("${workflow.projectDir}/subworkflows/local/isoform_quantification/bin/generate_DGE.R")
    )
    ch_versions = ch_versions.mix(DGE_GENERATION.out.versions)

    emit:
    idge = DGE_GENERATION.out.dge_matrices    // channel: [ val(meta), path(iDGE dir) ]
    versions = ch_versions                     // channel: [ versions.yml ]

}


process PROBABILITY_DISTRIBUTION {
    label 'process_single'
    tag "${meta}"

    conda "conda-forge::r-base=4.4.3 conda-forge::r-argparse=2.3.1 conda-forge::r-data.table=1.17.8 conda-forge::r-dplyr=1.1.4 conda-forge::r-fst=0.9.8 conda-forge::r-tidyr=1.3.2 conda-forge::r-ggplot2=3.5.1 conda-forge::r-scales=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-argparse_r-data.table_r-dplyr_pruned:eccd7a6e471b877c' :
        'community.wave.seqera.io/library/r-base_r-argparse_r-data.table_r-dplyr_pruned:6aa969298cb360df' }"

    input:
        tuple val(meta), path(mapped_reads)
        path(script)
        val(gene_frac)
        val(bins)

    output:
        tuple val("${meta}"), path("${meta}.prob"), path("${meta}.pdf"), emit: probabilities
        path('versions.yml'), emit: versions

    script:
    """
    Rscript ${script} \
        --repo ./ \
        --gene_frac ${gene_frac} \
        --bins ${bins} \
        --OUTPUT_PROB ${meta}.prob \
        --OUTPUT_PDF ${meta}.pdf \
        --threads ${task.cpus}
    
    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    Rscript -e 'v = c(base = as.character(getRversion()), sapply(c("data.table", "dplyr", "fst", "ggplot2"), function(p) as.character(packageVersion(p)))); cat(paste0("    r-", tolower(names(v)), ": ", v, "\\n"), file = "versions.yml", append = TRUE, sep = "")' 
    """
    stub:
    """
    touch XXX.prob
    touch XXX.pdf
    echo "R: stub" > versions.yml
    """
}


process FRAGMENT_PROBABILITY {
    label 'process_single'
    tag "${meta}, ${mapped_reads.baseName}"

    conda "conda-forge::r-base=4.4.3 conda-forge::r-argparse=2.3.1 conda-forge::r-data.table=1.17.8 conda-forge::r-dplyr=1.1.4 conda-forge::r-fst=0.9.8 conda-forge::r-tidyr=1.3.2 conda-forge::r-stringr=1.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-argparse_r-data.table_r-dplyr_pruned:beb90925512dcf02' :
        'community.wave.seqera.io/library/r-base_r-argparse_r-data.table_r-dplyr_pruned:33d5c4529c5ead95' }"

    input:
        tuple val(meta), path(probabilities), path(mapped_reads)
        path(script)

    output:
        tuple val("${meta}"), path("${meta}_${mapped_reads.baseName}.frag.fst"), emit: isoform_probabilities
        path('versions.yml'), emit: versions

    script:
    """
    Rscript ${script} \
        --reads ${mapped_reads} \
        --probs ${probabilities} \
        --output ${meta}_${mapped_reads.baseName}.frag.fst \
        --threads ${task.cpus}

    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    Rscript -e 'v = c(base = as.character(getRversion()), sapply(c("data.table", "dplyr", "fst"), function(p) as.character(packageVersion(p)))); cat(paste0("    r-", tolower(names(v)), ": ", v, "\\n"), file = "versions.yml", append = TRUE, sep = "")' 
    """
}


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



process EM_ALGORITHM {
    label 'process_single'
    tag "${meta}, ${cell_isoform_probabilities.baseName}"

    conda "conda-forge::r-base=4.4.3 conda-forge::r-argparse=2.3.1 conda-forge::r-data.table=1.17.8 conda-forge::r-dplyr=1.1.4 conda-forge::r-fst=0.9.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-argparse_r-data.table_r-dplyr_r-fst:e79dc183dca8c207' :
        'community.wave.seqera.io/library/r-base_r-argparse_r-data.table_r-dplyr_r-fst:c22aa6005af410db' }"

    input:
        tuple val(meta), path(cell_isoform_probabilities)
        path(script)

    output:
        tuple val("${meta}"), path("${cell_isoform_probabilities.baseName}_em_results.fst"), emit: em_results
        path('versions.yml'), emit: versions

    script:
    """
    Rscript ${script} \
        --cell ${cell_isoform_probabilities} \
        --output_path ${cell_isoform_probabilities.baseName}_em_results.fst

    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    Rscript -e 'v = c(base = as.character(getRversion()), sapply(c("data.table", "dplyr", "fst"), function(p) as.character(packageVersion(p)))); cat(paste0("    r-", tolower(names(v)), ": ", v, "\\n"), file = "versions.yml", append = TRUE, sep = "")' 
    """
}



process DGE_GENERATION {
    label 'process_single'
    tag "${meta}"

    conda "conda-forge::r-base=4.4.3 conda-forge::r-argparse=2.3.1 conda-forge::r-data.table=1.17.8 conda-forge::r-dplyr=1.1.4 conda-forge::r-fst=0.9.8 conda-forge::r-matrix=1.7_6 conda-forge::r-r.utils=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-argparse_r-data.table_r-dplyr_pruned:c4165245638f3c1b' :
        'community.wave.seqera.io/library/r-base_r-argparse_r-data.table_r-dplyr_pruned:513a1d31b2feb181' }"

    input:
        tuple val(meta), path(sample_dge_path), path(em_results)
        path(script)

    output:
        tuple val("${meta}"), path("${meta}_iDGE"), emit: dge_matrices
        path('versions.yml'), emit: versions

    script:
    """
    Rscript ${script} \
        --dge ${sample_dge_path} \
        --pred ./ \
        --output ${meta}_iDGE

    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    Rscript -e 'v = c(base = as.character(getRversion()), sapply(c("data.table", "fst", "Matrix"), function(p) as.character(packageVersion(p)))); cat(paste0("    r-", tolower(names(v)), ": ", v, "\\n"), file = "versions.yml", append = TRUE, sep = "")' 
    """
    stub:
    """
    mkdir ${meta}_iDGE
    echo "R: stub" > versions.yml
    """
}