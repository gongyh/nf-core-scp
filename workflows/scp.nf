/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowScp.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.hmm ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check optional parameters
ch_hmm = params.hmm ? file(params.hmm, checkIfExists: true) : file("$projectDir/assets/dummy_file.txt", checkIfExists: true)
ch_proteins = params.proteins ? Channel.fromPath(params.proteins) : []
ch_prodigal_tf = params.prodigal_tf ? Channel.fromPath(params.prodigal_tf) : []
ch_gtdb = params.gtdb ? file(params.gtdb, checkIfExists: true) : Channel.empty()
ch_extra_genomes = params.extra_genomes ? Channel.fromPath(params.extra_genomes) : Channel.empty()

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : file("$projectDir/assets/multiqc_custom_config.yaml", checkIfExists: true)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { CHECKM } from '../modules/local/checkm' addParams( options: modules['checkm'] )
include { GTDBTK } from '../modules/local/gtdbtk' addParams( options: modules['gtdbtk'] )
include { RETRIEVE_GENOMES } from '../modules/local/retrieve_genomes' addParams( options: modules['retrieve_genomes'] )
include { ROARY } from '../modules/local/roary' addParams( options: modules['roary'] )
include { ROARY2FRIPAN } from '../modules/local/roary2fripan' addParams( options: modules['roary2fripan'] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'  addParams( options: modules['trimgalore'] )
include { SPADES } from '../modules/nf-core/modules/spades/main'  addParams( options: modules['spades'] )
// tailored modules
include { PROKKA } from '../modules/local/prokka'  addParams( options: modules['prokka'] )
include { MULTIQC } from '../modules/local/multiqc' addParams( options: multiqc_options )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SCP {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    //
    // MODULE: Run TrimGalore
    //
    TRIMGALORE (
        INPUT_CHECK.out.reads
    )
    ch_software_versions = ch_software_versions.mix(TRIMGALORE.out.version.first().ifEmpty(null))

    //
    // MODULE: Run SPAdes
    //
    SPADES (
        TRIMGALORE.out.reads,
        ch_hmm
    )
    ch_software_versions = ch_software_versions.mix(SPADES.out.version.first().ifEmpty(null))

    ch_extra_genomes.map { it ->
        def meta = [:]; meta.id = it.name.replaceFirst(~/\.[^\.]+$/, ''); meta.single_end = false
        [meta, it]
    }
    .set {ch_extra_genomes_meta}

    //
    // MODULE: Run CheckM
    //
    CHECKM (
        SPADES.out.scaffolds.mix(ch_extra_genomes_meta)
    )
    ch_software_versions = ch_software_versions.mix(CHECKM.out.version.first().ifEmpty(null))

    CHECKM.out.completeness
        .map { meta, genome, completeness_tsv ->
            [meta, genome, WorkflowScp.getCompletenessContamination(completeness_tsv)]
        }
        .branch { meta, genome, checkm ->
            pass: checkm[0]
                return [meta, genome]
        }
        .set { ch_hq_genomes }

    //
    // MODULE: Run GTDB-Tk
    //
    GTDBTK (
        ch_hq_genomes.pass,
        ch_gtdb
    )
    ch_software_versions = ch_software_versions.mix(GTDBTK.out.version.first().ifEmpty(null))

    GTDBTK.out.scaffolds
        .map { meta, scaffolds, taxa -> [meta, WorkflowScp.getGenus(taxa), scaffolds] }
        .multiMap { meta, genus, scaffolds ->
            genome:
                def meta2 = meta.clone(); meta2.genus = genus
                [meta2, scaffolds]
            taxa:
                def meta3 = meta.clone(); meta3.id = ""; meta3.genus = genus
                meta3
        }
        .set {ch_gtdbtk_genomes}

    ch_gtdbtk_genomes.taxa.unique().set{ch_genus_meta}

    //
    // MODULE: Run retrieve_genomes
    //
    RETRIEVE_GENOMES (
        ch_genus_meta,
        ch_gtdb
    )

    RETRIEVE_GENOMES.out.scaffolds
        .map { meta, scaffolds ->
            def refs = []
            for (scf in scaffolds) {
                def meta2 = meta.clone()
                meta2.id = scf.name.replaceFirst(~/\.[^\.]+$/, '')
                refs << [meta2, scf]
            }
            refs
        }
        .flatten()
        .collate(2)
        .set {ch_gtdb_genomes}

    //
    // MODULE: Run Prokka
    //
    PROKKA (
        ch_gtdbtk_genomes.genome.mix(ch_gtdb_genomes),
        ch_proteins,
        ch_prodigal_tf
    )
    ch_software_versions = ch_software_versions.mix(PROKKA.out.version.first().ifEmpty(null))

    //
    // MODULE: Run Roary
    //
    ROARY (
        PROKKA.out.gff
            .map { meta, gff ->
                meta.id = ""; meta.single_end = false
                [meta, gff]
            }
            .groupTuple()
    )
    ch_software_versions = ch_software_versions.mix(ROARY.out.version.ifEmpty(null))

    //
    // MODULE: Run roary2fripan
    //
    ROARY2FRIPAN (
        ROARY.out.roary_pa
    )
    ch_software_versions = ch_software_versions.mix(ROARY2FRIPAN.out.version.ifEmpty(null))

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowScp.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    MULTIQC (
        ch_multiqc_config,
        ch_multiqc_custom_config,
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
        GET_SOFTWARE_VERSIONS.out.yaml.collect(),
        FASTQC.out.zip.collect{it[1]}.ifEmpty([]),
        TRIMGALORE.out.zip.collect{it[1]}.ifEmpty([]),
        TRIMGALORE.out.log.collect{it[1]}.ifEmpty([]),
        CHECKM.out.mqc_tsv.collect{it[1]}.ifEmpty([]),
        PROKKA.out.txt.collect{it[1]}.ifEmpty([])
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
