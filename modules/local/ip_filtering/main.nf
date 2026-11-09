process IP_FILTERING {
    label 'process_single'
    tag "${meta}, ${seqname}"

    conda "conda-forge::r-base=4.4.3 conda-forge::r-argparse=2.3.1 conda-forge::r-dplyr=1.2.0 conda-forge::r-data.table=1.17.8 conda-forge::r-fst=0.9.8 conda-forge::r-stringr=1.5.1 conda-forge::r-r.utils=2.13.0 bioconda::bioconductor-genomicranges=1.58.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-argparse_r-dplyr_r-data.table_pruned:42636a19761d44d9' :
        'community.wave.seqera.io/library/r-base_r-argparse_r-dplyr_r-data.table_pruned:122f7e640136e794' }"

    input:
        tuple val(meta), val(seqname), val(mapped_reads), path(ip_reference)
        path(script)
        val(distance_threshold)

    output:
        tuple val("${meta}"), val("${seqname}"), path("${mapped_reads.baseName}.ip_filtered.fst"), emit: filtered_mapped_reads
        path('versions.yml'), emit: versions

    script:
    """
    Rscript ${script} \
        --reads ${mapped_reads} \
        --ip_reference ${ip_reference} \
        --output ${mapped_reads.baseName}.ip_filtered.fst \
        --threshold_dist ${distance_threshold} \
        --threads ${task.cpus}
    
    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    Rscript -e 'v = c(base = as.character(getRversion()), sapply(c("data.table", "dplyr", "fst", "GenomicRanges"), function(p) as.character(packageVersion(p)))); cat(paste0("    r-", tolower(names(v)), ": ", v, "\\n"), file = "versions.yml", append = TRUE, sep = "")' 
    """
    stub:
    """
    touch XXX.ip_filtered.fst
    echo "R: stub" > versions.yml
    """
}
