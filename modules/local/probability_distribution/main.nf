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
        path("${meta}.warning.txt"), emit: warning, optional: true
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
