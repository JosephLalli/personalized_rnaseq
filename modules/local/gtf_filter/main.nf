process GTF_FILTER {
    tag "$fasta"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(fasta), path(gtf)

    output:
    tuple val(meta), path ("*.filtered.gtf"), emit: genome_gtf
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // filter_gtf.py is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    filter_gtf.py \\
        --gtf $gtf \\
        --fasta $fasta \\
        --prefix ${fasta.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta.baseName}.filtered.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
