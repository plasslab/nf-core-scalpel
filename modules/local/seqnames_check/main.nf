process SEQNAMES_CHECK {
    label 'process_single'
    tag "${meta.id}"

    conda "bioconda::samtools=1.25"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools:1.25--f798bb3cb27e958b' :
        'community.wave.seqera.io/library/samtools:1.25--a253eb9139fd1fcd' }"

    input:
        tuple val(meta), path(bam_file), path(bai_file), path(barcodes_file)
        val(annotated_seqnames)
        val(barcode_tag)

    output:
        tuple val(meta), path(bam_file, includeInputs: true), path(bai_file, includeInputs: true), path(barcodes_file, includeInputs: true), emit: checked_bams
        path('versions.yml'), emit: versions

    script:
    """
    #- chromosome names of BAM & annotation...
    samtools view -H ${bam_file} | awk -F'\\t' '\$1=="@SQ"{ for(i=2;i<=NF;i++) if(\$i ~ /^SN:/) print substr(\$i,4) }' | sort -u > bam_seqnames.txt
    printf "%s\\n" ${annotated_seqnames.join(' ')} | sort -u > annot_seqnames.txt
    comm -12 bam_seqnames.txt annot_seqnames.txt > shared_seqnames.txt
    if [ ! -s shared_seqnames.txt ]; then
        echo "ERROR -> sample '${meta.id}': no chromosome name shared between the BAM file and the annotation (GTF)." >&2
        echo "  BAM: \$(head -5 bam_seqnames.txt | tr '\\n' ' ')..." >&2
        echo "  GTF: \$(head -5 annot_seqnames.txt | tr '\\n' ' ')..." >&2
        echo "  Use references from the same source as the BAM (e.g. Ensembl '1' vs GENCODE/UCSC 'chr1')." >&2
        exit 1
    fi
    echo "${meta.id}: \$(wc -l < shared_seqnames.txt) chromosomes shared between BAM & annotation"

    #- barcode tag (CB for 10x / STARsolo, XC for Drop-seq)...
    nb_tagged=\$(samtools view ${bam_file} 2>/dev/null | head -n 10000 | grep -c "${barcode_tag}:Z:" || true)
    if [ "\$nb_tagged" -eq 0 ]; then
        echo "ERROR -> sample '${meta.id}': no '${barcode_tag}' cell barcode tag found in the first 10,000 BAM records." >&2
        echo "  STARsolo: run with '--outSAMattributes CB UB' (and --outSAMtype BAM SortedByCoordinate)." >&2
        exit 1
    fi

    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    echo "    samtools: \$(samtools --version | head -1 | sed 's/samtools //')" >> versions.yml
    """
    stub:
    """
    echo "R: stub" > versions.yml
    """
}
