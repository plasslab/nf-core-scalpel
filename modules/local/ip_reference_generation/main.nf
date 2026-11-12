process IP_REFERENCE_GENERATION {
    label 'process_medium'
    tag "${genome_fasta}"

    conda "conda-forge::r-base=4.4.3 conda-forge::r-argparse=2.3.1 conda-forge::r-dplyr=1.2.0 conda-forge::r-data.table=1.17.8 conda-forge::r-fst=0.9.8 conda-forge::r-stringr=1.5.1 conda-forge::r-r.utils=2.13.0 bioconda::bioconductor-genomicranges=1.58.0 conda-forge::gzip=1.15 conda-forge::findutils=4.10.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-argparse_r-dplyr_r-data.table_pruned:98cb79cffab48ef3' :
        'community.wave.seqera.io/library/r-base_r-argparse_r-dplyr_r-data.table_pruned:f802acc31df1be69' }"

    input:
        path(genome_fasta)
        path(genome_annotation)
        val(annotated_seqnames)
        path(script)

    output:
        path("*.ip.txt"), emit: ip_chrs
        path("ip_reference.track.gz"), emit: ip_reference
        path('versions.yml'), emit: versions

    script:
    """
    #- annotated chromosomes: GTF (start of the run, in parallel with salmon) or fixed isoform annotation...
    if [ -n "${genome_annotation}" ]; then
        gzip -cdf ${genome_annotation} | awk -F'\\t' '!/^#/ && NF >= 9 { print \$1 }' | sort -u > annot_seqnames.txt
    else
        printf "%s\\n" ${annotated_seqnames.join(' ')} | sort -u > annot_seqnames.txt
    fi

    #- one (possibly empty) file per annotated chromosome...
    while read chr; do : > \${chr}.ip.txt; done < annot_seqnames.txt

    #- single-line sequence of annotated chromosomes (plain or .gz FASTA)...
    gzip -cdf ${genome_fasta} | awk 'NR==FNR{ keep[\$1]=1; next } /^>/{ if(out != "") close(out); name=substr(\$1, 2); out = (name in keep) ? name".seq" : ""; next } out != ""{ gsub(/\\r/, ""); printf "%s", \$0 > out }' annot_seqnames.txt -

    if [ -z "\$(find . -name '*.seq' -size +0 | head -1)" ]; then
        echo "ERROR -> genome '${genome_fasta}': no chromosome name shared with the annotation (GTF)." >&2
        echo "  GTF: \$(head -5 annot_seqnames.txt | tr '\\n' ' ')..." >&2
        exit 1
    fi

    #- internal priming sites (genomic A-rich regions) per chromosome...
    ls *.seq | sed 's/\\.seq\$//' | xargs -P ${task.cpus} -I {} \
        Rscript ${script} \
            --sequence {}.seq \
            --seqname {} \
            --output {}.ip.txt
    rm -f *.seq

    #- genome-wide reference (reusable as --ip_reference)...
    cat \$(sort annot_seqnames.txt | sed 's/\$/.ip.txt/') | gzip -n > ip_reference.track.gz

    # Write versions.yml
    echo '"${task.process}":' > versions.yml
    Rscript -e 'v = c(base = as.character(getRversion()), sapply(c("argparse", "data.table"), function(p) as.character(packageVersion(p)))); cat(paste0("    r-", tolower(names(v)), ": ", v, "\\n"), file = "versions.yml", append = TRUE, sep = "")'
    """
    stub:
    """
    touch XXX.ip.txt
    echo "" | gzip > ip_reference.track.gz
    echo "R: stub" > versions.yml
    """
}
