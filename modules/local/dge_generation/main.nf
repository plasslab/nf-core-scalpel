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
