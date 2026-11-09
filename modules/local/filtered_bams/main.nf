process FILTERED_BAMS {
    label 'process_low'
    tag "${meta}"

    conda "bioconda::samtools=1.25 conda-forge::r-base=4.4.3 conda-forge::r-fst=0.9.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools_r-base_r-fst:2381996845f05fb0' :
        'community.wave.seqera.io/library/samtools_r-base_r-fst:1eb0eb46db0b220d' }"

    input:
        tuple val(meta), path(bams, stageAs: "bams/*"), path(filtered_reads)

    output:
        tuple val(meta), path("${meta}_filtered.bam"), path("${meta}_filtered.bam.bai"), emit: filtered_bams
        path('versions.yml'), emit: versions

    script:
    """
    #- read ids kept by SCALPEL (after all filtering steps)...
    Rscript -e 'ids = unique(unlist(lapply(list.files(pattern = "ip_filtered.fst\$"), function(f) fst::read_fst(f, columns = "read_id")\$read_id))); writeLines(ids, "ALL_RIDS.txt")'

    #- merge chromosomes & keep filtered reads (as original SCALPEL generation_filtered_bams)...
    samtools merge -f -o tmp.bam bams/*.bam -@ ${task.cpus}
    samtools view -b -N ALL_RIDS.txt tmp.bam > ${meta}_filtered.bam
    samtools index ${meta}_filtered.bam
    rm tmp.bam

    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    echo "    samtools: \$(samtools --version | head -1 | sed 's/samtools //')" >> versions.yml
    Rscript -e 'cat("    r-fst: ", as.character(packageVersion("fst")), "\\n", file = "versions.yml", append = TRUE, sep = "")' 
    """
    stub:
    """
    touch ${meta}_filtered.bam
    touch ${meta}_filtered.bam.bai
    echo "R: stub" > versions.yml
    """
}
