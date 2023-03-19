process SALMON_INDEX {
    tag "${meta.id}"
    label "process_medium"

    publishDir "${params.outdir}/salmon_indexes", 
                mode: params.publish_dir_mode, 
                saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
                overwrite: false,
                failOnError: true,
                enabled: {!file("${params.outdir}/salmon_indexes/${meta.id}_salmon_indexes").exists() }


    conda (params.enable_conda ? 'bioconda::salmon=1.5.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.5.2--h84f40af_0' :
        'quay.io/biocontainers/salmon:1.9.0--h7e5ed60_1' }"
    
    // if (file(task.publishDir.path[0].toString() + '/' +  meta.id + "_salmon_indexes").exists()){
    //     publishDir '', enabled: false 
    // }




    input:
    tuple val(meta), path (genome_fasta), path (genome_fasta_fai), path (genome_fasta_gzi), path (transcript_fasta), path (transcript_fasta_fai), path (transcript_fasta_gzi)

    output:
    tuple val(meta), path ("*_salmon_indexes")      , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cat_genome  = genome_fasta.toString().endsWith('.gz') ? "<(gunzip -c $genome_fasta)" : "$genome_fasta"
    def gentrome    = genome_fasta.toString().endsWith('.gz') ? "gentrome.fa.gz"             : "gentrome.fa"
    def previously_generated_file_path = params.force_resume ? task.publishDir.path[0] : ""    
    def is_active = !file("${params.outdir}/salmon_indexes/${meta.id}_salmon_indexes").exists() 
    """
    if [[ -d ${previously_generated_file_path}/${meta.id}_salmon_indexes ]]
    then
        ln -s ${previously_generated_file_path}/${meta.id}_salmon_indexes .
    else
        echo ${is_active}
        grep '^>' $cat_genome | cut -d ' ' -f 1 > decoys.txt
        sed -i.bak -e 's/>//g' decoys.txt
        cat $transcript_fasta $genome_fasta > $gentrome

        salmon \\
            index \\
            --threads $task.cpus \\
            -t $gentrome \\
            -d decoys.txt \\
            $args \\
            -i ${meta.id}_salmon_index
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
