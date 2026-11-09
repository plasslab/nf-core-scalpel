process IP_SPLITTING {
    label 'process_single'
    tag "${ip_reference}"

    conda "conda-forge::gawk=5.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/gawk:5.3.1--55cafbff4f74c897' :
        'community.wave.seqera.io/library/gawk:5.3.1--d5d074ad9f29be71' }"

    input:
        path(ip_reference)
        val(annotated_seqnames)

    output:
        path("*.ip.txt"), emit: ip_chrs
        path('versions.yml'), emit: versions

    script:
    def reader = ip_reference.name.endsWith('.tar.gz') ? "tar -xzOf ${ip_reference}" : "gzip -cdf ${ip_reference}"
    """
    #- one (possibly empty) file per annotated chromosome...
    printf "%s\\n" ${annotated_seqnames.join(' ')} | sort -u > annot_seqnames.txt
    while read chr; do : > \${chr}.ip.txt; done < annot_seqnames.txt

    #- split internal priming sites (.polyA.track / ipriming_sites.bed, plain, .gz or .tar.gz)...
    ${reader} | gawk -F'\\t' 'NR==FNR{ keep[\$1]=1; next } (\$1 in keep){ print > (\$1".ip.txt") }' annot_seqnames.txt -

    if [ -z "\$(find . -name '*.ip.txt' -size +0 | head -1)" ]; then
        echo "ERROR -> internal priming reference '${ip_reference}': no chromosome name shared with the annotation (GTF)." >&2
        echo "  GTF: \$(head -5 annot_seqnames.txt | tr '\\n' ' ')..." >&2
        exit 1
    fi

    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    echo "    gawk: \$(gawk --version | head -1 | sed 's/GNU Awk //; s/,.*//')" >> versions.yml
    """
    stub:
    """
    touch XXX.ip.txt
    echo "R: stub" > versions.yml
    """
}
