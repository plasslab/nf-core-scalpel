process BED_CONVERSION {
    label 'process_single'
    tag "${meta.id}, ${annotated_seqname}"

    conda "bioconda::samtools=1.25"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools:1.25--f798bb3cb27e958b' :
        'community.wave.seqera.io/library/samtools:1.25--a253eb9139fd1fcd' }"

    input:
        tuple val(meta), path(bam_file), path(bai_file), path(barcodes_file), val(annotated_seqname)
        val sequencing_platform
        val min_mapq

    output:
        tuple val("${annotated_seqname}"), val("${meta.id}"), path("${annotated_seqname}.bam"), emit: bed_files, optional: true
        path('versions.yml'), emit: versions

    script:
    def barcode_tag = sequencing_platform == 'dropseq' ? 'XC' : 'CB'
    def umi_tag = sequencing_platform == 'dropseq' ? 'XM' : 'UB'
    def has_barcodes = barcodes_file.name != 'NO_BARCODES'
    def barcode_input = has_barcodes ? (barcodes_file.toString().endsWith('.gz') ? "<(gunzip -c ${barcodes_file})" : "${barcodes_file}") : ""
    def filter_cmd = has_barcodes ? "-D ${barcode_tag}:${barcode_input}" : ""
    """
    #Filter reads (cells, tags) by chromosome
    samtools view \
        -b ${bam_file} \
        -X ${bai_file} ${filter_cmd} \
        -q ${min_mapq} \
        --keep-tag "${barcode_tag},${umi_tag}" ${annotated_seqname} | samtools sort -@ ${task.cpus} > tmp.bam
    #Remove all PCR duplicates ... (as original SCALPEL bam_splitting)
    samtools markdup tmp.bam ${annotated_seqname}_tmp.bam -r --barcode-tag ${barcode_tag} --barcode-tag ${umi_tag}
    rm tmp.bam
    
    # delete empty temporary BAM files to avoid downstream errors
    samtools view ${annotated_seqname}_tmp.bam | tail -1 > check.txt
    if [ ! -s check.txt ]; then
        rm ${annotated_seqname}_tmp.bam
        rm check.txt
    else
        mv ${annotated_seqname}_tmp.bam ${annotated_seqname}.bam
    fi

    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    echo "    samtools: \$(samtools --version | head -1 | sed 's/samtools //')" >> versions.yml
    """
    stub:
    """
    touch XXX.bam
    echo "R: stub" > versions.yml
    """
}
