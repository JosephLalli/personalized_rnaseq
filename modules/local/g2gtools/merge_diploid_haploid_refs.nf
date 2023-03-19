process MERGE_DIPLOID_HAPLOID_REFS {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_0' }"

    input:
    tuple val(meta), 
        path (haploid_transcript_fasta), path (diploid_transcript_fasta), 
        path (haploid_genome_fasta), path (diploid_genome_fasta),
        path (haploid_genome_fasta_fai), path (diploid_genome_fasta_fai),
        path (haploid_genome_fasta_gzi), path (diploid_genome_fasta_gzi),
        path (haploid_gtf), path (diploid_gtf),
        path (haploid_vci), path (diploid_vci),
        path (haploid_vci_tbi), path (diploid_vci_tbi)

    output:
    tuple val(meta), path ("*.genome.fa.gz"), path('*.genome.fa.gz.fai'), path('*.genome.fa.gz.gzi'), emit: merged_refseq
    tuple val(meta), path ("*.transcriptome.fa.gz"), path('*.transcriptome.fa.gz.fai'), path('*.transcriptome.fa.gz.gzi'), emit: merged_transcriptome
    tuple val(meta), path ("*.gtf"), emit: merged_gtf
    tuple val(meta), path ("*.vci.gz"), path("*.vci.gz.tbi"), emit: merged_vci
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // zstdcat -T ${task.cpus} $diploid_genome_fasta $haploid_genome_fasta | bgzip -c -@ ${task.cpus} > ${prefix}.genome.fa.gz
    // zstdcat -T ${task.cpus} $diploid_transcript_fasta $haploid_transcript_fasta | bgzip -c -@ ${task.cpus} > ${prefix}.transcriptome.fa.gz
    def previously_generated_file_path = params.force_resume ? task.publishDir.path[0] : ""
    """
    cat $diploid_genome_fasta $haploid_genome_fasta > ${prefix}.genome.fa.gz && \\
    samtools faidx ${prefix}.genome.fa.gz & \\
    \\
    cat $diploid_transcript_fasta $haploid_transcript_fasta > ${prefix}.transcriptome.fa.gz && \\
    samtools faidx ${prefix}.transcriptome.fa.gz & \\
    \\
    cat $diploid_gtf $haploid_gtf > ${prefix}.gtf & \\
    \\
    zcat $haploid_vci | grep '##CONTIG' > contigs.txt && \\
    contig_insertion_line=\$(zcat $diploid_vci | grep -n \"#CHROM\" | cut -f1 -d:) && \\
    contig_insertion_line=\$((\$contig_insertion_line-1)) && \\
    zcat $diploid_vci | sed \"\$contig_insertion_line r contigs.txt\" > ${prefix}.vci && \\
    zcat $haploid_vci | grep -v '#' >> ${prefix}.vci && \\
    \\
    bgzip ${prefix}.vci && tabix -p vcf ${prefix}.vci.gz & \\
    wait

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1)
    END_VERSIONS
    """
}
