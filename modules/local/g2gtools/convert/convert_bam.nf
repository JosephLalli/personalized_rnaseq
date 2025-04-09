process CONVERT_BAM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'docker.io/jlalli/g2gtools:3.1-792b2da' }"

    input:
    tuple val(meta), path(vci), path(vci_tbi), path (bam), path(bai)

    output:
    tuple val(meta), path("*bam")   , emit: bam
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // and now, for the hackiest of hacks to get around the fact that
    // some of my test files have no variants in chrMT, leading to no
    // chrMT gtf records being converted in a mistaken belief that the
    // contig has been deleted.
    // TODO: figure out a better way!
    def previously_generated_file_path = params.force_resume ? task.publishDir.path[0] : ""
    """
    g2gtools convert -c ${vci} \\
                    -i ${bam} \\
                    -f bam \\
                    --reverse \\
                    -o ${prefix}_reference_coords.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version | tail -n 1)
    END_VERSIONS
    """
}
