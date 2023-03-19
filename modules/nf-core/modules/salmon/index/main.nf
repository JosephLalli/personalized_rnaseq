process SALMON_INDEX {
    tag "${meta.id}"
    label "process_medium"

    conda (params.enable_conda ? 'bioconda::salmon=1.9.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.9.0--h7e5ed60_1' :
        'quay.io/biocontainers/salmon:1.9.0--h7e5ed60_1' }"

    input:
    tuple val(meta), path (genome_fasta), path (transcript_fasta)

    output:
    tuple val(meta), path ("${meta.id}_index")     , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cat_genome  = genome_fasta.toString().endsWith('.gz') ? "<(gunzip -c $genome_fasta)" : "$genome_fasta"
    def gentrome    = genome_fasta.toString().endsWith('.gz') ? "gentrome.fa.gz"             : "gentrome.fa"
    

    """
    grep '^>' $cat_genome | cut -d ' ' -f 1 > decoys.txt
    sed -i.bak -e 's/>//g' decoys.txt
    cat $transcript_fasta $genome_fasta > $gentrome

    salmon \\
        index \\
        --threads $task.cpus \\
        -t $gentrome \\
        -d decoys.txt \\
        $args \\
        -i ${meta.id}_index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
