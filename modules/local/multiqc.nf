// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::multiqc=1.11--pyhdfd78af_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
    }

    input:
    path multiqc_config
    path multiqc_custom_config
    path workflow_summary
    path software_versions
    path ('fastqc/*')
    path ('trimgalore/fastqc/*')
    path ('cutadapt/*')
    path ('prokka/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "*.version.txt"       , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    multiqc -f --config $multiqc_custom_config $options.args .
    multiqc --version | sed -e "s/multiqc, version //g" > ${software}.version.txt
    """
}
