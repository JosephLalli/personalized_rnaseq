process SALMON_TX2GENE {
    tag "$ref_gtf"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    path ( "salmon/*" )
    path ( ref_gtf )

    output:
    path "*.tsv"       , emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_tx2gene.py \\
        --gtf $ref_gtf \\
        --salmon salmon \\
        --id $params.gtf_group_features \\
        --extra $params.gtf_extra_attributes \\
        -o salmon_tx2gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
