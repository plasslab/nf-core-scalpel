
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow APA_CHARACTERIZATION {

    take:
    ch_idge // channel: [ val(meta), path(iDGE dir) ] from ISOFORM_QUANTIFICATION
    ch_filtered_reads // channel: [ val(meta), val(seqname), path(ip_filtered.fst) ] from READS_PROCESSING
    ch_chr_bams // channel: [ val(seqname), val(meta), path(chr.bam) ] from READS_PROCESSING

    main:
    ch_versions = channel.empty()

    // Optional clusters annotation (orig.ident, Barcodes, Clusters)
    def clusters_file = params.clusters ? file(params.clusters, checkIfExists: true) : file("${projectDir}/assets/NO_CLUSTERS")

    // Differential isoform usage between clusters (or samples if no clusters provided)
    DIFFERENTIAL_ISOFORM_USAGE(
        ch_idge.map{ it = it[1] }.collect(),
        clusters_file,
        file("${projectDir}/subworkflows/local/apa_characterization/bin/differential_isoform_usage.R"),
        file("${projectDir}/subworkflows/local/apa_characterization/bin/scalpel_library.R")
    )
    ch_versions = ch_versions.mix(DIFFERENTIAL_ISOFORM_USAGE.out.versions)

    // BAM file of the reads kept by SCALPEL for each sample (e.g. coverage visualization)
    ch_chr_bams.map{
        seqname, meta, bam -> tuple(meta, bam)
    }.groupTuple(by: 0).join(
        ch_filtered_reads.map{ meta, seqname, reads -> tuple(meta, reads) }.groupTuple(by: 0)
    ).set{ ch_bams_and_reads }
    FILTERED_BAMS(ch_bams_and_reads)
    ch_versions = ch_versions.mix(FILTERED_BAMS.out.versions)

    emit:
    diu_table = DIFFERENTIAL_ISOFORM_USAGE.out.diu_table   // path: DIU_table.csv
    seurat_obj = DIFFERENTIAL_ISOFORM_USAGE.out.seurat_obj // path: iDGE_seurat.RDS
    filtered_bams = FILTERED_BAMS.out.filtered_bams        // channel: [ val(meta), bam, bai ]
    versions = ch_versions                                 // channel: [ versions.yml ]
}


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
        path(script)
        path(library)

    output:
        path("DIU_table.csv"), emit: diu_table
        path("iDGE_seurat.RDS"), emit: seurat_obj
        path('versions.yml'), emit: versions

    script:
    def clusters_arg = clusters.name == 'NO_CLUSTERS' ? 'NULL' : "${clusters}"
    """
    Rscript ${script} \
        --idge idge \
        --clusters ${clusters_arg} \
        --output DIU_table.csv \
        --output_rds iDGE_seurat.RDS

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


process FILTERED_BAMS {
    label 'process_low'
    tag "${meta}"

    conda "bioconda::samtools=1.23.1 conda-forge::r-base=4.4.3 conda-forge::r-fst=0.9.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools_r-base_r-fst:2381996845f05fb0' :
        'community.wave.seqera.io/library/samtools_r-base_r-fst:122f35a5d549c8c8' }"

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
