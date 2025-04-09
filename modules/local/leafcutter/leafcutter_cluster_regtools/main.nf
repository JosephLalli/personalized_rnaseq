process LEAFCUTTER_CLUSTERINTRONS {
    tag "cluster_juncfiles"
    label 'process_medium'

    conda (params.enable_conda ? "python=3.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1':
        'quay.io/biocontainers/python:3.9--1'}"

    input:
    path (junction_files)
    path (junction_files_txt)

    output:
    path ("leafcutter_perind*.gz"), emit: leafcutter_perind_counts
    path ("versions.yml")         , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    leafcutter_cluster_regtools \\
        --juncfiles $junction_files_txt \\
        $args \\
        -o leafcutter \\
    && \\
    zcat leafcutter_perind_numers.counts.gz | \\
    sed 's/^/phenotype_id /' | \\
    sed 's/.sorted//g' | \\
    sed -e 's/ /\t/g' | \\
    gzip -c > leafcutter_perind_numers.counts.formatted.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leafcutter: \$(echo \$(leafcutter_cluster_regtools --version 2>&1))
    END_VERSIONS
    """

    // """
    // ls
    // echo *.junc
    // echo *.junc | tr ' ' '\n' > list_of_junc_files.txt && \\
    // leafcutter_cluster_regtools \\
    //     --juncfiles list_of_junc_files.txt \\
    //     -o leafcutter \\
    // \\
    // && \\
    // gcat leafcutter_perind_numers.counts.gz | \\
    // sed 's/^/phenotype_id /' | \\
    // sed 's/.sorted//g' | \\
    // sed -e 's/ /\t/g' | \\
    // gzip -c > leafcutter_perind_numers.counts.formatted.gz

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     leafcutter: \$(echo \$(leafcutter_cluster_regtools --version 2>&1))
    // END_VERSIONS
    // """
}
