// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CHECKM {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::checkm-genome==1.1.3--py_1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/checkm-genome:1.1.3--py_1"
    } else {
        container "quay.io/biocontainers/checkm-genome:1.1.3--py_1"
    }

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("genomes/${genome}"), path("${prefix}_completeness.tsv"),   emit: completeness
    path  '*.version.txt'   , emit: version
    path  "mqc_${prefix}_completeness.tsv", emit: mqc_tsv

    script:
    def software    = getSoftwareName(task.process)
    prefix          = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    mkdir -p checkm_out genomes
    cp $genome genomes/

    checkm lineage_wf \\
        --extension fa \\
        --threads $task.cpus \\
        $options.args \\
        --tab_table \\
        --file ${prefix}_completeness.tsv \\
        genomes checkm_out

    cat ${prefix}_completeness.tsv | grep -v Completeness | cut -f1,2,12,13,14 > mqc_${prefix}_completeness.tsv

    echo \$(checkm -h 2>&1) | grep CheckM | sed 's/^.*CheckM v//; s/ .*\$//' > ${software}.version.txt
    """
}
