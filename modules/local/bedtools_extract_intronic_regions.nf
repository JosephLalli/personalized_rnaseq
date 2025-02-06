process EXTRACT_INTRONIC_REGIONS {
    tag "${meta.id}"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path (gtf), path(fai)

    output:
    tuple val(meta), path ('*intron.bed')     , emit: intronic_regions
    tuple val(meta), path ('*exon_regions.bed')       , emit: exonic_regions
    tuple val(meta), path ('*intergenic.bed') , emit: intergenic_regions
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $fai | cut -f 1-2 | sort > sorted_sizes.txt
    cat $gtf | cut -f 1-8 | grep gene | sed "s/^chrM\\t/chrMT\\t/" | cut -f 1,4-8 > ${prefix}.gene_regions.unsorted.bed
    cat $gtf | cut -f 1-8 | grep exon | sed "s/^chrM\\t/chrMT\\t/" | cut -f 1,4-8 > ${prefix}.exon_regions.unsorted.bed

    bedtools sort -i ${prefix}.gene_regions.unsorted.bed > ${prefix}.gene_regions.bed
    bedtools sort -i ${prefix}.exon_regions.unsorted.bed | cut -f 1-3 > ${prefix}.exon_regions.bed

    bedtools complement -i ${prefix}.gene_regions.bed -g sorted_sizes.txt > ${prefix}.intergenic.bed

    bedtools complement -i <(cat ${prefix}.exon_regions.bed ${prefix}.intergenic.bed | sort -k1,1 -k2,2n) -g sorted_sizes.txt > ${prefix}.intron.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
