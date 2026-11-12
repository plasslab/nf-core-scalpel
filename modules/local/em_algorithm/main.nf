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
        val(max_iter)
        val(tol)

    output:
        tuple val("${meta}"), path("${cell_isoform_probabilities.baseName}_em_results.fst"), emit: em_results
        path('versions.yml'), emit: versions

    script:
    """
    Rscript ${script} \
        --cell ${cell_isoform_probabilities} \
        --output_path ${cell_isoform_probabilities.baseName}_em_results.fst \
        --max_iter ${max_iter} \
        --tol ${tol}

    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    Rscript -e 'v = c(base = as.character(getRversion()), sapply(c("data.table", "dplyr", "fst"), function(p) as.character(packageVersion(p)))); cat(paste0("    r-", tolower(names(v)), ": ", v, "\\n"), file = "versions.yml", append = TRUE, sep = "")' 
    """
}
