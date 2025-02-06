process CONVERT_TO_SAMPLE_COORDS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'jlalli/g2gtools:0.2.9' }"

    input:
    tuple val(meta), path(vci)
    path (files_to_convert)
    val (suffix)

    output:
    tuple val(meta), path("*.tab")   , emit: tab
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.refcoords"
    """
    g2gtools convert -c ${vci} \\
                     -f ${suffix} \\
                     -o ${prefix}.${suffix} \\
                     -i ${files_to_convert}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version 2>&1)
    END_VERSIONS
    """
}
