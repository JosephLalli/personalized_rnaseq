process GTF2DB {
    tag "$meta.id"
    label 'io_heavy'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'docker.io/jlalli/g2gtools:3.1-792b2da' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.db")   , emit: db
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    g2gtools gtf2db ${args} --gtf ${gtf} --db ${prefix}.db


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version | tail -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version | tail -n 1)
    END_VERSIONS
    """
}
