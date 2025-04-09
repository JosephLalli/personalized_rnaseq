process SALMON_INDEX {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.3--h6dccd9a_2' :
        'biocontainers/salmon:1.10.3--h6dccd9a_2' }"

    input:
    tuple val(meta), path (genome_fasta), path (transcript_fasta)

    output:
    tuple val(meta), path ("${meta.id}_index"), emit: index
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def decoys = ''
    def fasta = transcript_fasta

    if (genome_fasta){
        if (genome_fasta.toString().endsWith('.gz')) {
            genome_fasta = "<(gunzip -c $genome_fasta)"
        }
        decoys='-d decoys.txt'
        fasta='gentrome.fa'
    }
    if (transcript_fasta.toString().endsWith('.gz')) {
        transcript_fasta = "<(gunzip -c $transcript_fasta)"
    }
    
    """
    if [ -n '$genome_fasta' ]; then
        grep '^>' $genome_fasta | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 | sed 's/>//g' > decoys.txt
        cat $transcript_fasta $genome_fasta > $fasta
    fi

    salmon \\
        index \\
        --threads $task.cpus \\
        -t $fasta \\
        $decoys \\
        $args \\
        -i ${meta.id}_index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """

    stub:
    """
    mkdir ${meta.id}_index
    touch ${meta.id}_index/complete_ref_lens.bin
    touch ${meta.id}_index/ctable.bin
    touch ${meta.id}_index/ctg_offsets.bin
    touch ${meta.id}_index/duplicate_clusters.tsv
    touch ${meta.id}_index/info.json
    touch ${meta.id}_index/mphf.bin
    touch ${meta.id}_index/pos.bin
    touch ${meta.id}_index/pre_indexing.log
    touch ${meta.id}_index/rank.bin
    touch ${meta.id}_index/refAccumLengths.bin
    touch ${meta.id}_index/ref_indexing.log
    touch ${meta.id}_index/reflengths.bin
    touch ${meta.id}_index/refseq.bin
    touch ${meta.id}_index/seq.bin
    touch ${meta.id}_index/versionInfo.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
