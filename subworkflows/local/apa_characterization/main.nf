
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow APA_CHARACTERIZATION {

    take:
    ch_idge // channel: [ val(meta), path(iDGE dir) ] from ISOFORM_QUANTIFICATION

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

    emit:
    diu_table = DIFFERENTIAL_ISOFORM_USAGE.out.diu_table   // path: DIU_table.csv
    seurat_obj = DIFFERENTIAL_ISOFORM_USAGE.out.seurat_obj // path: iDGE_seurat.RDS
    versions = ch_versions                                 // channel: [ versions.yml ]
}


process DIFFERENTIAL_ISOFORM_USAGE {
    label 'process_medium'
    tag "${idge_dirs.size()} samples"

    conda "conda-forge::r-base=4.4.3 conda-forge::r-seurat=5.5.1 conda-forge::r-argparse=2.3.1 conda-forge::r-data.table=1.17.8 conda-forge::r-dplyr=1.1.4 conda-forge::r-tidyr=1.3.2 conda-forge::r-tibble conda-forge::r-pbapply=1.7_4 conda-forge::r-stringr=1.5.1 conda-forge::r-r.utils=2.13.0"
    container "docker.io/satijalab/seurat:5.0.0"

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

    # Write versions.yml from bash
    echo "!${task.process}:" > versions.yml
    """
    stub:
    """
    touch DIU_table.csv
    touch iDGE_seurat.RDS
    echo "R: stub" > versions.yml
    """
}
