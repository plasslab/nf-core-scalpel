process DIFFERENTIAL_ISOFORM_USAGE {
    label 'process_medium'
    tag "${idge_dirs.size()} samples"

    conda "conda-forge::r-base=4.4.3 conda-forge::r-seurat=5.5.1 conda-forge::r-argparse=2.3.1 conda-forge::r-data.table=1.17.8 conda-forge::r-dplyr=1.1.4 conda-forge::r-tidyr=1.3.2 conda-forge::r-tibble=3.3.1 conda-forge::r-pbapply=1.7_4 conda-forge::r-stringr=1.5.1 conda-forge::r-r.utils=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-seurat_r-argparse_r-data.table_pruned:154299b4c8b4d80a' :
        'community.wave.seqera.io/library/r-base_r-seurat_r-argparse_r-data.table_pruned:c108c4be46eec081' }"

    input:
        path(idge_dirs, stageAs: "idge/*")
        path(clusters)
        path(annotations, stageAs: "annots/*")
        val(diu_test)
        val(min_delta)
        val(permutations)
        path(script)
        path(library)

    output:
        path("DIU_table.csv"), emit: diu_table
        path("iDGE_seurat.RDS"), emit: seurat_obj
        path('versions.yml'), emit: versions

    script:
    def clusters_arg = clusters.name == 'NO_CLUSTERS' ? 'NULL' : "${clusters}"
    """
    #- one BLAS thread per worker (genes are tested in parallel on ${task.cpus} cores)...
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

    Rscript ${script} \
        --idge idge \
        --clusters ${clusters_arg} \
        --output DIU_table.csv \
        --output_rds iDGE_seurat.RDS \
        --annotations annots \
        --test ${diu_test} \
        --min_delta ${min_delta} \
        --permutations ${permutations} \
        --cpus ${task.cpus}

    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    Rscript -e 'v = c(base = as.character(getRversion()), sapply(c("Seurat", "data.table", "dplyr"), function(p) as.character(packageVersion(p)))); cat(paste0("    r-", tolower(names(v)), ": ", v, "\\n"), file = "versions.yml", append = TRUE, sep = "")' 
    """
    stub:
    """
    touch DIU_table.csv
    touch iDGE_seurat.RDS
    echo "R: stub" > versions.yml
    """
}
