process EXTRACT {
    tag "$meta.id"
    label 'io_heavy'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'jlalli/g2gtools:0.2.9' }"

    input:
    tuple val(meta), path(personal_fasta), path(personal_fasta_fai), path(personal_fasta_gzi), path(db)

    output:
    // tuple val(meta), path("*.exons.fasta.gz")   , emit: exons
    tuple val(meta), path("*.transcripts.fa.gz")   , emit: transcripts
    // tuple val(meta), path("*.genes.fasta.gz")   , emit: genes
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def previously_generated_file_path = params.force_resume ? task.publishDir.path[0] : ""
    """
    if [[ -s ${previously_generated_file_path}/${prefix}.transcripts.fa.gz ]]
    then
        ln -s ${previously_generated_file_path}/${prefix}.transcripts.fa.gz .
    else
        g2gtools extract -i ${personal_fasta} \\
                        -db ${db} \\
                        --transcripts ${args} \\
                        > ${prefix}.transcripts.fa
    
        bgzip ${prefix}.transcripts.fa
    fi
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version 2>&1)
    END_VERSIONS
    """
}

    // g2gtools extract -i ${personal_fasta} \\
    //                  -db ${db} \\
    //                  --exons ${args} | gzip -c > ${prefix}.exons.fasta.gz
    
    // g2gtools extract -i ${personal_fasta} \\
    //                  -db ${db} \\
    //                  --genes ${args} | gzip -c > ${prefix}.genes.fasta.gz
