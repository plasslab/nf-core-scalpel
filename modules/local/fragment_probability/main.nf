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
