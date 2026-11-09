
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { resolveInputPath   } from '../utils_nfcore_scalpel_pipeline'

workflow READS_PROCESSING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_all_transcripts // path: all transcripts preprocessed from ANNOTATION_PROCESSING

    main:

    ch_versions = channel.empty()
    def barcode_tag = params.sequencing_platform == 'dropseq' ? 'XC' : 'CB'

    // Parse optional barcodes whitelist (paths relative to the whitelist csv)
    if (params.barcodes_whitelist) {
        def whitelist_dir = file(params.barcodes_whitelist).parent
        ch_barcodes_whitelist = channel
            .fromPath(params.barcodes_whitelist)
            .splitCsv(header: true)
            .map { row -> tuple(row.sample, file(resolveInputPath(row.barcodes, whitelist_dir, row.sample, 'barcodes'))) }
    } else {
        ch_barcodes_whitelist = channel.empty()
    }

    // Extract BAM files, index & defaults Barcodes (resolved in PIPELINE_INITIALISATION)
    // files: [fastq_1, fastq_2, cranger_path, bam, bai, dge, barcodes]
    if( params.sequencing_platform == 'dropseq' ) {
        // Drop-seq: cells = barcodes of the DGE count matrix
        DGE_BARCODES(
            ch_samplesheet.map{ meta, files -> tuple(meta, file(files[5])) }
        )
        ch_versions = ch_versions.mix(DGE_BARCODES.out.versions)
        ch_samplesheet.map{
            meta, files -> tuple(meta.id, meta, files[3], files[4])
        }.join(
            DGE_BARCODES.out.barcodes.map{ meta, barcodes -> tuple(meta.id, barcodes) }
        ).map{
            sample_id, meta, bam, bai, barcodes -> tuple(meta, bam, bai, barcodes)
        }.set{ ch_bam_files }
    } else {
        // 10x Chromium: cells = filtered barcodes (Cell Ranger / STARsolo)
        ch_samplesheet.map{
            meta, files -> tuple(meta, files[3], files[4], file(files[6] ?: "${projectDir}/assets/NO_BARCODES"))
        }.set{ ch_bam_files }
    }

    // Join barcodes whitelist with BAM files if provided (override default barcodes)
    ch_bam_files.map{
        meta, bam, bai, barcodes -> tuple(meta.id, meta, bam, bai, barcodes)
    }.join(ch_barcodes_whitelist, by: 0, remainder: true).filter{
        it[1] != null
    }.map{
        sample_id, meta, bam, bai, barcodes, custom_barcodes -> tuple(meta, bam, bai, custom_barcodes ?: barcodes)
    }.branch{
        indexed: it[2]
        not_indexed: true
    }.set{ ch_bam_with_barcodes }

    // Index BAM files without index (e.g. STARsolo)
    SAMTOOLS_INDEX(
        ch_bam_with_barcodes.not_indexed.map{ meta, bam, bai, barcodes -> tuple(meta, file(bam)) }
    )
    ch_bam_with_barcodes.not_indexed.map{
        meta, bam, bai, barcodes -> tuple(meta.id, meta, bam, barcodes)
    }.join(
        SAMTOOLS_INDEX.out.index.map{ meta, bai -> tuple(meta.id, bai) }
    ).map{
        sample_id, meta, bam, barcodes, bai -> tuple(meta, bam, bai.toString(), barcodes)
    }.mix(ch_bam_with_barcodes.indexed).set{ ch_indexed_bams }

    // Check chromosome names (BAM vs annotation) & barcode tags before processing (#20, #26, #36)
    SEQNAMES_CHECK(
        ch_indexed_bams,
        ch_all_transcripts.map{ it = it[0] }.collect(),
        barcode_tag
    )
    ch_versions = ch_versions.mix(SEQNAMES_CHECK.out.versions)

    // Convert BAM to BED files
    BED_CONVERSION(
        SEQNAMES_CHECK.out.checked_bams.combine(ch_all_transcripts.map{ it=it[0] }),
        params.sequencing_platform
    )
    ch_versions = ch_versions.mix(BED_CONVERSION.out.versions)
    

    // Map reads in BED files to transcripts from annotation
    ch_all_transcripts.combine(
        BED_CONVERSION.out.bed_files, by: 0
    ).map{
        it = tuple(it[2], it[0], it[1], it[3])
    }.set{ reads_and_annots }


    READS_MAPPING_AND_FILTERING(
        reads_and_annots,
        file("${projectDir}/subworkflows/local/reads_processing/bin/reads_merging_v3.R"),
        params.distance_profile
    )
    ch_versions = ch_versions.mix(READS_MAPPING_AND_FILTERING.out.versions)


    // Split internal priming reference by chromosome (once, instead of per task)
    IP_SPLITTING(
        file(params.ip_reference, checkIfExists: true),
        ch_all_transcripts.map{ it = it[0] }.collect()
    )
    ch_versions = ch_versions.mix(IP_SPLITTING.out.versions)
    IP_SPLITTING.out.ip_chrs.flatten().map{
        it = tuple(it.name - '.ip.txt', it)
    }.set{ ch_ip_chrs }

    // Filter reads associated to internal priming positions
    READS_MAPPING_AND_FILTERING.out.mapped_reads.map{
        meta, seqname, mapped_reads -> tuple(seqname, meta, mapped_reads)
    }.combine(ch_ip_chrs, by: 0).map{
        seqname, meta, mapped_reads, ip_chr -> tuple(meta, seqname, mapped_reads, ip_chr)
    }.set{ reads_and_ips }

    IP_FILTERING(
        reads_and_ips,
        file("${projectDir}/subworkflows/local/reads_processing/bin/ip_filtering.R"),
        params.distance_ip
    )
    ch_versions = ch_versions.mix(IP_FILTERING.out.versions)


    emit:
    reads = IP_FILTERING.out.filtered_mapped_reads
    chr_bams = BED_CONVERSION.out.bed_files     // channel: [ val(seqname), val(meta), path(chr.bam) ]
    versions = ch_versions                     // channel: [ versions.yml ]
}



process DGE_BARCODES {
    label 'process_single'
    tag "${meta.id}"

    conda "conda-forge::gawk=5.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/gawk:5.3.1--55cafbff4f74c897' :
        'community.wave.seqera.io/library/gawk:5.3.1--d5d074ad9f29be71' }"

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


process SEQNAMES_CHECK {
    label 'process_single'
    tag "${meta.id}"

    conda "bioconda::samtools=1.23.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools:1.23.1--ddf9171cf37418c2' :
        'community.wave.seqera.io/library/samtools:1.23.1--51a9c8a523bed249' }"

    input:
        tuple val(meta), val(bam_file), val(bai_file), path(barcodes_file)
        val(annotated_seqnames)
        val(barcode_tag)

    output:
        tuple val(meta), val(bam_file), val(bai_file), path(barcodes_file), emit: checked_bams
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


process BED_CONVERSION {
    label 'process_single'
    tag "${meta.id}, ${annotated_seqname}"

    conda "bioconda::samtools=1.23.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools:1.23.1--ddf9171cf37418c2' :
        'community.wave.seqera.io/library/samtools:1.23.1--51a9c8a523bed249' }"

    input:
        tuple val(meta), val(bam_file), val(bai_file), path(barcodes_file), val(annotated_seqname)
        val sequencing_platform

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


