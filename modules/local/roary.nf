// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ROARY {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::roary==3.13.0--pl526h516909a_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/roary:3.13.0--pl526h516909a_0"
    } else {
        container "quay.io/biocontainers/roary:3.13.0--pl526h516909a_0"
    }

    input:
    path prokka_gff3

    output:
    path "${prefix}/gene_presence_absence.csv" , emit: roary_pa
    path "${prefix}"
    path  '*.version.txt'                      , emit: version

    when:
    prokka_gff3.size() > 1

    script:
    def software    = getSoftwareName(task.process)
    prefix          = options.suffix ?: software
    """
    echo $prokka_gff3
    roary -p $task.cpus -f $prefix *.gff

    echo \$(roary -w 2>/dev/null) > ${software}.version.txt
    """
}