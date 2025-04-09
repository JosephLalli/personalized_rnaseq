process STAR_GENOMEGENERATE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'nf-core/htslib_samtools_star_gawk:311d422a50e6d829'

    input:
    tuple val(meta), path (fasta), path (fai), path (gtf), path(vcf)

    output:
    tuple val(meta), path ("${meta.id}"), emit: index
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def args_list   = args.tokenize()
    def memory      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def include_gtf = gtf ? "--sjdbGTFfile $gtf" : ''
    if (args_list.contains('--genomeSAindexNbases')) {
        """
        mkdir ${meta.id}
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir ${meta.id}/ \\
            --genomeFastaFiles $fasta \\
            $include_gtf \\
            --runThreadN $task.cpus \\
            $memory \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    } else {
        """
        NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' ${fai}`

        mkdir ${meta.id}
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir ${meta.id}/ \\
            --genomeFastaFiles $fasta \\
            $include_gtf \\
            --runThreadN $task.cpus \\
            --genomeSAindexNbases \$NUM_BASES \\
            $memory \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    }

    stub:
    if (gtf) {
        """
        mkdir ${meta.id}
        touch ${meta.id}/Genome
        touch ${meta.id}/Log.out
        touch ${meta.id}/SA
        touch ${meta.id}/SAindex
        touch ${meta.id}/chrLength.txt
        touch ${meta.id}/chrName.txt
        touch ${meta.id}/chrNameLength.txt
        touch ${meta.id}/chrStart.txt
        touch ${meta.id}/exonGeTrInfo.tab
        touch ${meta.id}/exonInfo.tab
        touch ${meta.id}/geneInfo.tab
        touch ${meta.id}/genomeParameters.txt
        touch ${meta.id}/sjdbInfo.txt
        touch ${meta.id}/sjdbList.fromGTF.out.tab
        touch ${meta.id}/sjdbList.out.tab
        touch ${meta.id}/transcriptInfo.tab

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    } else {
        """
        mkdir ${meta.id}
        touch ${meta.id}/Genome
        touch ${meta.id}/Log.out
        touch ${meta.id}/SA
        touch ${meta.id}/SAindex
        touch ${meta.id}/chrLength.txt
        touch ${meta.id}/chrName.txt
        touch ${meta.id}/chrNameLength.txt
        touch ${meta.id}/chrStart.txt
        touch ${meta.id}/genomeParameters.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """
    }
}
