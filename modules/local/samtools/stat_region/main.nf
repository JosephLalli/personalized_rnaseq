process SAMTOOLS_STATS_REGION {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(bedfile), path (fasta)
    val (region_name)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${input}"
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    cat ${bedfile} | awk 'BEGIN{OFS="\\t"} FNR==0{print; next} {for (i=2;i<=NF;i++)\$i=\$i+1}1' > regions_file.txt

    samtools \\
        stats \\
        --threads ${task.cpus-1} \\
        -t regions_file.txt \\
        $args \\
        ${reference} \\
        ${input} \\
        > ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${input}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
