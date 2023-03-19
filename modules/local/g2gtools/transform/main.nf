process TRANSFORM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'jlalli/g2gtools:0.3.1' }"

    input:
    tuple val(meta), path(vci), path(vci_tbi), path(patched_fasta), path(patched_fasta_fai), path(patched_fasta_gzi)
    
    output:
    tuple val(meta), path("*.fa.gz"), path("*.fa.gz.fai"), path("*.fa.gz.gzi"), emit: fasta
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def previously_generated_file_path = params.force_resume ? task.publishDir.path[0] : ""
    """
    if [[ -s ${previously_generated_file_path}/${prefix}.fa.gz ]]
    then
        ln -s ${previously_generated_file_path}/${prefix}.fa.gz .
        ln -s ${previously_generated_file_path}/${prefix}.fa.gz.fai .
        ln -s ${previously_generated_file_path}/${prefix}.fa.gz.gzi .
    else
        g2gtools transform ${args} -c ${vci} \\
                        -i ${patched_fasta} \\
                        -p ${task.cpus} \\
                        -o ${prefix}.fa \\
                        --bgzip
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version 2>&1)
    END_VERSIONS
    """
}
