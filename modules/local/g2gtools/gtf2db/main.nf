process GTF2DB {
    tag "$meta.id"
    label 'io_heavy'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'jlalli/g2gtools:0.3.1' }"

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
    def previously_generated_file_path = params.force_resume ? task.publishDir.path[0] : ""
    """
    if [[ -f ${previously_generated_file_path}/${prefix}.db ]]
    then
        ln -s ${previously_generated_file_path}/${prefix}.db .
    else
        g2gtools gtf2db ${args} -i ${gtf} -o ${prefix}.db
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version 2>&1)
    END_VERSIONS
    """
}
