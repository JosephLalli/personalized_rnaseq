process CONVERT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'docker.io/jlalli/g2gtools:3.1-792b2da' }"

    input:
    tuple val(meta), path(vci), path(vci_tbi), path (infile)
    val (suffix)

    output:
    tuple val(meta), path("*$suffix")   , emit: converted_file
    path  "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_specific"
    // and now, for the hackiest of hacks to get around the fact that
    // some of my test files have no variants in chrMT, leading to no
    // chrMT gtf records being converted in a mistaken belief that the
    // contig has been deleted.
    // TODO: figure out a better way!
    def append_chrMT = ""
    if (meta.ploidy && (meta.ploidy == 'haploid') && (suffix == 'gtf')){
        append_chrMT = "cat *.unmapped | grep ^chrMT >> ${prefix}.gtf || true"
    }

    """
    g2gtools convert --vci ${vci} \\
                     --in ${infile} \\
                     --file-format $suffix \\
                     $args \\
                     --out ${prefix}.$suffix

    ${append_chrMT}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version | tail -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_specific"
    """
    touch ${prefix}.$suffix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version | tail -n 1)
    END_VERSIONS
    """
}
