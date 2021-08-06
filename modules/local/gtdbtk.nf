// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GTDBTK {
    tag "$meta.id"
    label 'process_high'
    maxForks 1
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gtdbtk==1.5.1--pyhdfd78af_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gtdbtk:1.5.1--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/gtdbtk:1.5.1--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(genome)
    path gtdb

    output:
    path "${prefix}"
    tuple val(meta), path("genome/${prefix}.fa"), path('taxa.txt')  , emit: scaffolds
    path  '*.version.txt'                  , emit: version

    script:
    def software    = getSoftwareName(task.process)
    prefix          = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    export GTDBTK_DATA_PATH=$gtdb

    mkdir -p genome
    cp $genome genome/${prefix}.fa
    echo 'g__' >  taxa.txt

    if [[ -f $gtdb ]]; then
        mkdir -p $prefix
        echo 'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli' > taxa.txt
        touch ${software}.version.txt
    else
        gtdbtk classify_wf \\
            $options.args \\
            --genome_dir genome \\
            --extension fa \\
            --out_dir $prefix \\
            --cpus $task.cpus

        if [ -f ${prefix}/*.summary.tsv ]; then
            cut -f2 ${prefix}/*.summary.tsv | grep -v classification > taxa.txt
        fi

        echo \$(gtdbtk -h 2>&1) | grep 'GTDB-Tk' | sed 's/^.*GTDB-Tk v//; s/ .*\$//' > ${software}.version.txt
    fi
    """
}
