process CONVERT_TO_SECOND_SAMPLE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'docker.io/jlalli/g2gtools:3.1-792b2da' }"

    input:
    tuple val(meta), path(file_to_convert)
    path (vci_one)
    path (vci_two)
    val (suffix)

    output:
    tuple val(meta), path("*_remapped_to_*")   , emit: gtf
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    g2gtools convert -c ${vci_one} \\
                    -i ${file_to_convert} \\
                    -f ${suffix} \\
                    --reverse
                    -o tmp.refcoords

    g2gtools convert -c ${vci_two} \\
                    -i tmp.refcoords \\
                    -f ${suffix} \\
                    -o ${prefix}_remapped_to_${vci_two.baseName()}.${suffix}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version | tail -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.refcoords"
    """
    touch ${prefix}_remapped_to_${vci_two.baseName()}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version | tail -n 1)
    END_VERSIONS
    """
}
