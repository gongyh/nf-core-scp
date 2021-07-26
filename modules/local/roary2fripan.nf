// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ROARY2FRIPAN {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::roary2fripan.py==0.1--hdfd78af_2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/roary2fripan.py:0.1--hdfd78af_2"
    } else {
        container "quay.io/biocontainers/roary2fripan.py:0.1--hdfd78af_2"
    }

    input:
    path roary_pa

    output:
    path "${prefix}.*"
    path  '*.version.txt' , emit: version

    script:
    def software    = getSoftwareName(task.process)
    prefix          = options.suffix ?: software
    """
    roary2fripan.py --input $roary_pa ${prefix}

    echo \$(roary2fripan.py --version 2>&1 | grep roary2fripan) | sed 's/roary2fripan.py v//' > ${software}.version.txt
    """
}
