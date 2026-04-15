

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
    ch_samplesheet.map{ meta, files -> tuple(meta.id, files[2]) }.join(
        EM_ALGORITHM.out.em_results.groupTuple(by: 0)
    ).map{
        sample_id, sample_dge_path, em_results -> tuple(sample_id, file(sample_dge_path + "/outs/filtered_feature_bc_matrix/"), em_results)
    }.set{ dge_inputs }

    // Generate DGE matrix for each sample
    DGE_GENERATION(
        dge_inputs,
        file("${workflow.projectDir}/subworkflows/local/isoform_quantification/bin/generate_DGE.R")
    )

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
    
}


process PROBABILITY_DISTRIBUTION {
    label 'process_single'
    tag "${meta}"

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-argparse_r-data.table_r-dplyr_r-fst_pruned:fcc4fd7ca6b2f422' :
        'community.wave.seqera.io/library/r-argparse_r-data.table_r-dplyr_r-fst_pruned:71f85a249e6be098' }"

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
    
    # Write versions.yml for bedtools version in bash
    echo "!${task.process}:" > versions.yml
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

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-argparse_r-data.table_r-dplyr_r-fst_pruned:fcc4fd7ca6b2f422' :
        'community.wave.seqera.io/library/r-argparse_r-data.table_r-dplyr_r-fst_pruned:71f85a249e6be098' }"

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

    # Write versions.yml for bedtools version in bash
    echo "!${task.process}:" > versions.yml
    """
}


process CELL_DATA_GENERATION {
    label 'process_single'
    tag "${meta}"

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-argparse_r-data.table_r-dplyr_r-fst_pruned:fcc4fd7ca6b2f422' :
        'community.wave.seqera.io/library/r-argparse_r-data.table_r-dplyr_r-fst_pruned:71f85a249e6be098' }"

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
            paste0('    dplyr: "', packageVersion("dplyr"), '"'),
            '"'
        ), "versions.yml")
    """
}



process EM_ALGORITHM {
    label 'process_single'
    tag "${meta}, ${cell_isoform_probabilities.baseName}"

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-argparse_r-data.table_r-dplyr_r-fst_pruned:fcc4fd7ca6b2f422' :
        'community.wave.seqera.io/library/r-argparse_r-data.table_r-dplyr_r-fst_pruned:71f85a249e6be098' }"

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

    # Write versions.yml for bedtools version in bash
    echo "!${task.process}:" > versions.yml
    """
}



process DGE_GENERATION {
    label 'process_single'
    tag "${meta}"

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-argparse_r-data.table_r-dplyr_r-fst_pruned:8f91f8e3b3c44c32' :
        'community.wave.seqera.io/library/r-argparse_r-data.table_r-dplyr_r-fst_pruned:9feb5ef5dc170dba' }"

    input:
        tuple val(meta), file(sample_dge_path), path(em_results)
        path(script)

    output:
        tuple val("${meta}"), path("${meta}_dge_matrix.txt"), emit: dge_matrices
        path('versions.yml'), emit: versions

    script:
    """
    Rscript ${script} \
        --dge ${sample_dge_path} \
        --pred ./ \
        --output ${meta}_dge_matrix.txt

    # Write versions.yml from bash
    echo "!${task.process}:" > versions.yml
    """
}