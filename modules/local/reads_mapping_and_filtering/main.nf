process READS_MAPPING_AND_FILTERING {
    label 'process_single'
    tag "${meta}, ${seqname}"

    conda "conda-forge::r-base=4.4.3 conda-forge::r-argparse=2.3.1 conda-forge::r-dplyr=1.2.1 conda-forge::r-data.table=1.17.8 conda-forge::r-fst=0.9.8 bioconda::bedops=2.4.42 conda-forge::r-r.utils=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-argparse_r-dplyr_r-data.table_pruned:7d88a113d2eb0730' :
        'community.wave.seqera.io/library/r-base_r-argparse_r-dplyr_r-data.table_pruned:9e005428b384d2ec' }"

    input:
        tuple val(meta), val(seqname), path(annots), val(reads)
        path(script)
        val(distance_thr)
    
    output:
        tuple val("${meta}"), val("${seqname}"), path("${meta}_${reads.baseName}.fst"), emit: mapped_reads
        path('versions.yml'), emit: versions

    script:
    """
    bam2bed --all-reads --split < ${reads} > ${reads.baseName}.tmp
    Rscript ${script} \
        --reads ${reads.baseName}.tmp \
        --annotations ${annots} \
        --output ${meta}_${reads.baseName}.fst \
        --distance_threshold ${distance_thr}
    rm -f ${reads.baseName}.tmp
    
    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    Rscript -e 'v = c(base = as.character(getRversion()), sapply(c("data.table", "dplyr", "fst"), function(p) as.character(packageVersion(p)))); cat(paste0("    r-", tolower(names(v)), ": ", v, "\\n"), file = "versions.yml", append = TRUE, sep = "")' 
    echo "    bedops: \$(bam2bed --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)" >> versions.yml
    """
    stub:
    """
    touch XXX.fst
    echo "R: stub" > versions.yml
    """
}
