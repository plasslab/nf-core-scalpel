//
// Subworkflow with functionality specific to the nf-core/scalpel pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    samplesheet       //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/scalpel ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { "    https://doi.org/${it.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/scalpel/blob/master/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Always print pipeline-specific options regardless of whether they differ from defaults
    //
    log.info """
Pipeline options
  sequencing_platform : ${params.sequencing_platform}
  libtype             : ${params.libtype}
  distance_profile    : ${params.distance_profile}
  distance_3end       : ${params.distance_3end}
  distance_ip         : ${params.distance_ip}
  gene_fraction       : ${params.gene_fraction}
  binsize             : ${params.binsize}
  barcodes_whitelist  : ${params.barcodes_whitelist ?: 'not provided'}
------------------------------------------------------
""".stripIndent()

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.samplesheet
    //

    Channel
        .fromList(samplesheetToList(params.samplesheet, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2, cranger_path, dge ->
                def base_dir = file(params.samplesheet).parent
                def sample_path = resolveInputPath(cranger_path, base_dir, meta.id, 'cranger_path')
                if( params.sequencing_platform == 'dropseq' && !dge ) {
                    error("Samplesheet error -> sample '${meta.id}': column 'dge' (Drop-seq DGE count matrix) is required with --sequencing_platform dropseq.")
                }
                def dge_path = dge ? resolveInputPath(dge, base_dir, meta.id, 'dge') : ''
                // files: [fastq_1, fastq_2, cranger_path, bam, bai, dge_dir, barcodes]
                [ meta.id, meta, [
                    resolveInputPath(fastq_1, base_dir, meta.id, 'fastq_1'),
                    resolveInputPath(fastq_2, base_dir, meta.id, 'fastq_2'),
                    sample_path
                ] + resolveSampleLayout(meta.id, sample_path, params.sequencing_platform, dge_path) ]
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Resolve samplesheet paths (relative paths -> samplesheet directory) & check existence
//
def resolveInputPath(input_path, base_dir, sample_id, column) {
    def is_absolute = input_path ==~ /^(\/|[a-zA-Z][a-zA-Z0-9+.-]*:\/\/).*/
    def resolved_path = is_absolute ? file(input_path) : base_dir.resolve(input_path)
    if( !resolved_path.exists() ) {
        error("Samplesheet error -> sample '${sample_id}', column '${column}': '${input_path}' not found (resolved to '${resolved_path}'). Relative paths are resolved from the samplesheet directory.")
    }
    return resolved_path.toString()
}

//
// Resolve sample layout (Cell Ranger, Cell Ranger ARC multiome, STARsolo, Drop-seq BAM + DGE)
// returns [bam, bai, dge, barcodes] ('' when not available)
//
def resolveSampleLayout(sample_id, sample_path, platform, dge_path) {
    def sample_dir = file(sample_path)
    def pick = { candidates -> candidates.collect{ sample_dir.resolve(it) }.find{ it.exists() } }
    def find_bai = { bam_path -> ["${bam_path}.bai", "${bam_path}".replaceAll('\\.bam$', '.bai')].findAll{ it != "${bam_path}" }.collect{ file(it) }.find{ it.exists() } }

    // Drop-seq: BAM file (XC/XM tags) + DGE count matrix (Drop-seq tools); cells = DGE barcodes
    if( platform == 'dropseq' ) {
        if( sample_dir.isDirectory() || !sample_path.endsWith('.bam') ) {
            error("Samplesheet error -> sample '${sample_id}': with --sequencing_platform dropseq, column 'cranger_path' must be the BAM file (got '${sample_path}').")
        }
        def bai_dropseq = find_bai(sample_path)
        return [ sample_path, bai_dropseq ? bai_dropseq.toString() : '', dge_path, '' ]
    }

    // 10x Chromium: Cell Ranger (GEX / ARC multiome) or STARsolo output directory
    def bam_candidates = [
        'outs/possorted_genome_bam.bam',            // Cell Ranger count
        'outs/gex_possorted_bam.bam',               // Cell Ranger ARC (multiome GEX)
        'possorted_genome_bam.bam',
        'Aligned.sortedByCoord.out.bam'             // STARsolo
    ]
    def dge_candidates = [
        'outs/filtered_feature_bc_matrix',          // Cell Ranger (GEX / ARC)
        'filtered_feature_bc_matrix',
        'Solo.out/Gene/filtered',                   // STARsolo
        'Solo.out/GeneFull/filtered'
    ]
    def bam = pick(bam_candidates)
    def dge_dir = pick(dge_candidates)
    if( !bam ) {
        error("Samplesheet error -> sample '${sample_id}': no BAM file found in '${sample_path}'. Expected one of: ${bam_candidates.join(', ')}.")
    }
    if( !dge_dir ) {
        error("Samplesheet error -> sample '${sample_id}': no filtered count matrix found in '${sample_path}'. Expected one of: ${dge_candidates.join(', ')}.")
    }
    def barcodes = ['barcodes.tsv.gz', 'barcodes.tsv'].collect{ dge_dir.resolve(it) }.find{ it.exists() }
    if( !barcodes ) {
        error("Samplesheet error -> sample '${sample_id}': no barcodes.tsv(.gz) in '${dge_dir}'.")
    }
    def bai = find_bai(bam)
    return [ bam.toString(), bai ? bai.toString() : '', dge_dir.toString(), barcodes.toString() ]
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""


    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Send notification via hook URL
//
def imNotification(summary_params, hook_url) {
    def summary = [:]
    summary_params
        .keySet()
        .sort()
        .each { group ->
            summary << summary_params[group]
        }

    def misc_fields = [:]
    misc_fields['start']          = workflow.start
    misc_fields['complete']       = workflow.complete
    misc_fields['scriptfile']     = workflow.scriptFile
    misc_fields['scriptid']       = workflow.scriptId
    if (workflow.repository) {
        misc_fields['repository'] = workflow.repository
    }
    if (workflow.commitId) {
        misc_fields['commitid']   = workflow.commitId
    }
    if (workflow.revision) {
        misc_fields['revision']   = workflow.revision
    }
    misc_fields['nxf_version']    = workflow.nextflow.version
    misc_fields['nxf_build']      = workflow.nextflow.build
    misc_fields['nxf_timestamp']  = workflow.nextflow.timestamp

    def msg_fields = [:]
    msg_fields['version']      = getWorkflowVersion()
    msg_fields['runName']      = workflow.runName
    msg_fields['success']      = workflow.success
    msg_fields['dateComplete'] = workflow.complete
    msg_fields['duration']     = workflow.duration
    msg_fields['exitStatus']   = workflow.exitStatus
    msg_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    msg_fields['errorReport']  = (workflow.errorReport ?: 'None')
    msg_fields['commandLine']  = workflow.commandLine.replaceFirst(/ +--hook_url +[^ ]+/, "")
    msg_fields['projectDir']   = workflow.projectDir
    msg_fields['summary']      = summary << misc_fields

    // Render the JSON template
    def engine       = new groovy.text.GStringTemplateEngine()
    // Different JSON depending on the service provider
    // Defaults to "Adaptive Cards" (https://adaptivecards.io), except Slack which has its own format
    def json_path     = hook_url.contains("hooks.slack.com") ? "slackreport.json" : "adaptivecard.json"
    def hf            = new File("${workflow.projectDir}/assets/${json_path}")
    def json_template = engine.createTemplate(hf).make(msg_fields)
    def json_message  = json_template.toString()

    // POST
    def post = new URL(hook_url).openConnection()
    post.setRequestMethod("POST")
    post.setDoOutput(true)
    post.setRequestProperty("Content-Type", "application/json")
    post.getOutputStream().write(json_message.getBytes("UTF-8"))
    def postRC = post.getResponseCode()
    if (!postRC.equals(200)) {
        log.warn(post.getErrorStream().getText())
    }
}
