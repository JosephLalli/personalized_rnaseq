process REGTOOLS_BAMTOJUNC {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::regtools=0.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/regtools:0.6.1--hd03093a_0':
        'quay.io/biocontainers/regtools:0.6.1--hd03093a_0'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path ("*.junc"), emit: junctions
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }

    """
    regtools junctions extract \\
        -s $strandedness \\
        -a 8 \\
        $args \\
        -o ${prefix}.junc \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        regtools: \$(echo \$(regtools 2>&1) | sed 's/.*Version:\s*//; s/Usage.*//')

    END_VERSIONS
    """
}

        // -m ${params.leafcutter_min_intron_length} \\
        // -M ${params.leafcutter_max_intron_length} \\
