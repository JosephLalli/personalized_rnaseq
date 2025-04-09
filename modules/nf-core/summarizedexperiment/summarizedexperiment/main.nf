process SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0' :
        'biocontainers/bioconductor-summarizedexperiment:1.32.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(matrix_files), path(rowdata), path(coldata)

    output:
    tuple val(meta), path("*.rds")              , emit: rds
    tuple val(meta), path("*.R_sessionInfo.log"), emit: log
    path "versions.yml"                         , emit: versions


    script:
    template 'summarizedexperiment.r'

    stub:
    """
    touch ${meta.id}.SummarizedExperiment.rds
    touch ${meta.id}.R_sessionInfo.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-summarizedexperiment: \$(Rscript -e "library(SummarizedExperiment); cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """
}
