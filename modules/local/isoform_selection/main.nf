process ISOFORM_SELECTION{
    label 'process_single'
    tag "${gtf_entries.baseName}"

    conda "conda-forge::r-base=4.4.3 conda-forge::r-argparse=2.3.1 conda-forge::r-tidyr=1.3.2 conda-forge::r-dplyr=1.1.4 conda-forge::r-data.table=1.17.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-argparse_r-tidyr_r-dplyr_r-data.table:91bda81fa3909649' :
        'community.wave.seqera.io/library/r-base_r-argparse_r-tidyr_r-dplyr_r-data.table:ec5de7e1233d91ec' }"

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


    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    Rscript -e 'v = c(base = as.character(getRversion()), sapply(c("data.table", "dplyr"), function(p) as.character(packageVersion(p)))); cat(paste0("    r-", tolower(names(v)), ": ", v, "\\n"), file = "versions.yml", append = TRUE, sep = "")' 
    """
    stub:
    """
    touch ${gtf_entries.baseName}.annots.gz
    touch stats.txt
    echo "R: stub" > versions.yml
    """
}
