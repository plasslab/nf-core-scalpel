process DGE_BARCODES {
    label 'process_single'
    tag "${meta.id}"

    conda "conda-forge::gawk=5.3.1 conda-forge::gzip=1.15 conda-forge::findutils=4.10.0 conda-forge::tar=1.35"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/gawk_gzip_findutils_tar:8851949c471d65d2' :
        'community.wave.seqera.io/library/gawk_gzip_findutils_tar:82f953873e50fc8f' }"

    input:
        tuple val(meta), path(dge)

    output:
        tuple val(meta), path("${meta.id}_barcodes.txt"), emit: barcodes
        path('versions.yml'), emit: versions

    script:
    """
    #- cell barcodes = DGE header (GENE, barcode_1, barcode_2, ...)
    { gzip -cdf ${dge} || true; } | head -n 1 | tr '\\t' '\\n' | tail -n +2 > ${meta.id}_barcodes.txt
    if [ ! -s ${meta.id}_barcodes.txt ]; then
        echo "ERROR -> sample '${meta.id}': no cell barcode found in the header of the DGE file '${dge}'." >&2
        exit 1
    fi
    echo "${meta.id}: \$(wc -l < ${meta.id}_barcodes.txt) cells in DGE"

    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    echo "    gawk: \$(gawk --version | head -1 | sed 's/GNU Awk //; s/,.*//')" >> versions.yml
    """
    stub:
    """
    touch ${meta.id}_barcodes.txt
    echo "R: stub" > versions.yml
    """
}
