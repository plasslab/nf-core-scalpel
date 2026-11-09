
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scalpel_pipeline'
include { ANNOTATION_PROCESSING  } from '../subworkflows/local/annotation_processing/main.nf'
include { READS_PROCESSING       } from '../subworkflows/local/reads_processing/main.nf'
include { ISOFORM_QUANTIFICATION } from '../subworkflows/local/isoform_quantification/main.nf'
include { APA_CHARACTERIZATION   } from '../subworkflows/local/apa_characterization/main.nf'


workflow SCALPEL {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    genome_fasta // channel: genome FASTA file read in from --genome
    transcriptome_fasta // channel: transcriptome FASTA file read in from --transcriptome
    genome_annotation // channel: genome annotation GTF file read in from --gtf
    main:

    ch_versions = channel.empty()
    
    if( params.isoform_annotation ) {
        // Fixed annotation from a previous run (<outdir>/isoform/*.annots.gz), e.g. adding samples (#31)
        channel.fromPath("${params.isoform_annotation}/*.annots.gz", checkIfExists: true).map{
            it = tuple(it.name - '.annots.gz', it)
        }.set{ ch_all_transcripts }
    } else {
        // Annotation processing (pseudobulk quantification & isoform selection) from all samples
        ANNOTATION_PROCESSING(
            ch_samplesheet,
            genome_fasta,
            transcriptome_fasta,
            genome_annotation)
        ANNOTATION_PROCESSING.out.all_transcripts.set{ ch_all_transcripts }
    }

    READS_PROCESSING(
        ch_samplesheet,
        ch_all_transcripts
    )

    ISOFORM_QUANTIFICATION(
        ch_samplesheet,
        // group by sample_id [Index 0]
        READS_PROCESSING.out.reads.groupTuple(by: 0)
    )

    // Differential isoform usage (DIU_table.csv & iDGE_seurat.RDS)
    APA_CHARACTERIZATION(
        ISOFORM_QUANTIFICATION.out.idge
    )

    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'scalpel_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    versions = ch_versions
}