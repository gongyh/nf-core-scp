// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CHECKM {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::checkm-genome==1.1.3--py_1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/checkm-genome:1.1.3--py_1"
    } else {
        container "quay.io/biocontainers/checkm-genome:1.1.3--py_1"
    }

    input:
    path genome_files

    output:
    path "genome_completeness.tsv", emit: completeness
    path  '*.version.txt'         , emit: version

    script:
    def software    = getSoftwareName(task.process)
    """
    mkdir -p checkm_out
    checkm lineage_wf \\
        --extension fa \\
        --threads $task.cpus \\
        $options.args \\
        --tab_table \\
        --file genome_completeness.tsv \\
        . checkm_out

    echo \$(checkm -h 2>&1) | grep CheckM | sed 's/^.*CheckM v//; s/ .*\$//' > ${software}.version.txt
    """
}
