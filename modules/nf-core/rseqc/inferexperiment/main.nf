process RSEQC_INFEREXPERIMENT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(bam), path (bed)

    output:
    tuple val(meta), path("*.infer_experiment.txt"), emit: txt
    path  "versions.yml"                           , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    infer_experiment.py \\
        -i $bam \\
        -r $bed \\
        $args \\
        > ${prefix}.infer_experiment.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(infer_experiment.py --version | sed -e "s/infer_experiment.py //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.infer_experiment.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(infer_experiment.py --version | sed -e "s/infer_experiment.py //g")
    END_VERSIONS
    """
}
