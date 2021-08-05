// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RETRIEVE_GENOMES {
    tag "$meta.genus"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    val(meta)
    path gtdb

    output:
    tuple val(meta), path("${prefix}/*.fa")      , emit: scaffolds
    path  '*.version.txt'                        , emit: version

    script:
    def software    = getSoftwareName(task.process)
    prefix          = options.suffix ? "${meta.genus}${options.suffix}" : "${meta.genus}"
    """
    realDB=$gtdb
    if [[ -f $gtdb ]]; then
        mkdir -p db && gzip -cd $gtdb | tar xvf /dev/stdin -C ./db
        realDB=\$PWD/db
    fi
    cat \$realDB/taxonomy/gtdb_taxonomy.tsv | grep \";${meta.genus};\" | cut -f1 | cut -c 4- > refs.id
    mkdir -p ${prefix}
    for id in `cat refs.id`; do
        x=(`grep \${id}_genomic.fna.gz \$realDB/fastani/genome_paths.tsv`)
        gzip -cd \$realDB/fastani/\${x[1]}\${x[0]} > ${prefix}/\${id}.fa
    done

    touch ${software}.version.txt
    """
}
