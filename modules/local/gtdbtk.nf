// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GTDBTK {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::gtdbtk==1.5.1--pyhdfd78af_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gtdbtk:1.5.1--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/gtdbtk:1.5.1--pyhdfd78af_0"
    }

    input:
    path genome
    path gtdb

    output:
    path "${prefix}"
    path  '*.version.txt'         , emit: version

    script:
    def software    = getSoftwareName(task.process)
    prefix          = options.suffix ?: software
    """
    export GTDBTK_DATA_PATH=$gtdb
    gtdbtk classify_wf \\
        $options.args \\
        --genome_dir . \\
        --extension fa \\
        --out_dir $prefix \\
        --cpus $task.cpus

    echo \$(gtdbtk -h 2>&1) | grep 'GTDB-Tk' | sed 's/^.*GTDB-Tk v//; s/ .*\$//' > ${software}.version.txt
    """
}
