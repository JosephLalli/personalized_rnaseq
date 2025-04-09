process PATCH {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'docker.io/jlalli/g2gtools:3.1-792b2da' }"

    input:
    tuple val(meta), path(vci), path(vci_tbi), path (ref_fasta), path(fai), path(gzi)

    output:
    tuple val(meta), path("*.fa.gz"), path("*.fa.gz.fai"), path("*.fa.gz.gzi"), emit: reference_with_snps
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    g2gtools patch ${args} --vci ${vci} \\
                --fasta ${ref_fasta} \\
                --out ${prefix}.patched.fa \\
                --bgzip && \\
    samtools faidx ${prefix}.patched.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version | tail -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.patched.fa.gz
    touch  ${prefix}.patched.fa.gz.fai
    touch  ${prefix}.patched.fa.gz.gzi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version | tail -n 1)
    END_VERSIONS
    """
}

