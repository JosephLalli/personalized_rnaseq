process MERGE_DIPLOID_HAPLOID_REFS {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_0' }"

    input:
    tuple val(meta),
        path (transcript_fastas, arity:'1..*'),
        path (genome_fastas, arity:'1..*'),
        path (genome_fasta_fais, arity:'1..*'),
        path (genome_fasta_gzis, arity:'1..*'),
        path (gtfs, arity:'1..*'),
        path (vcis, arity:'1..*'),
        path (vci_tbis, arity:'1..*')

// Todo: ensure that input files are sorted in diploid,haploid order
// Especially important for vcis

    output:
        tuple val(meta), path ("*.genome.fa.gz"), path('*.genome.fa.gz.fai'), emit: merged_refseq
        tuple val(meta), path ("*.transcriptome.fa.gz"), path('*.transcriptome.fa.gz.fai'), emit: merged_transcriptome
        tuple val(meta), path ("*.gtf"), emit: merged_gtf
        tuple val(meta), path ("*.vci.gz"), path("*.vci.gz.tbi"), emit: merged_vci
        path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    diploid_vci = vcis[0]
    haploid_vci = vcis[1]
    if (genome_fastas.size() > 1) {
        """
        cat ${genome_fastas.join(' ')} > ${prefix}.genome.fa.gz && \\
        samtools faidx ${prefix}.genome.fa.gz & \\
        \\
        cat ${transcript_fastas.join(' ')} > ${prefix}.transcriptome.fa.gz && \\
        samtools faidx ${prefix}.transcriptome.fa.gz & \\
        \\
        cat ${gtfs.join(' ')} > ${prefix}.gtf & \\
        \\
        zcat ${vcis.join(' ')} | grep '##CONTIG' > contigs.txt && \\
        contig_insertion_line=\$(zcat $diploid_vci | grep -n \"#CHROM\" | cut -f1 -d:) && \\
        contig_insertion_line=\$((\$contig_insertion_line-1)) && \\
        zcat $diploid_vci | sed \"\$contig_insertion_line r contigs.txt\" > ${prefix}.vci && \\
        zcat $haploid_vci | grep -v '#' >> ${prefix}.vci && \\
        \\
        bgzip ${prefix}.vci && tabix -p vcf ${prefix}.vci.gz & \\
        wait

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        mv ${genome_fastas[0]} ${prefix}.genome.fa.gz
        mv ${genome_fasta_fais[0]} ${prefix}.genome.fa.gz.fai
        mv ${genome_fasta_gzis[0]} ${prefix}.genome.fa.gz.gzi
        mv ${transcript_fastas[0]} ${prefix}.transcriptome.fa.gz
        samtools faidx ${prefix}.transcriptome.fa.gz
        mv ${gtfs[0]} ${prefix}.gtf
        mv ${vcis[0]} ${prefix}.vci.gz
        mv ${vci_tbis[0]} ${prefix}.vci.gz.tbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genome.fa.gz
    touch ${prefix}.genome.fa.gz.fai
    touch ${prefix}.genome.fa.gz.gzi
    touch ${prefix}.transcriptome.fa.gz
    touch ${prefix}.transcriptome.fa.gz.fai
    touch ${prefix}.transcriptome.fa.gz.gzi
    touch ${prefix}.gtf
    touch ${prefix}.vci.gz
    touch ${prefix}.vci.gz.tbi
    touch contigs.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
