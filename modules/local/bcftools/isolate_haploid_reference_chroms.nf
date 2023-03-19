process SPLIT_OUT_HAPLOID_CHROMS {
    tag "sex_specific_references"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'nanozoo/samtools:1.14--d8fb865'}" // has both python and samtools

    input:
    path (ref_fasta)
    path (ref_fai)

    output:
    tuple path ("no_Y.fa.gz"), path('no_Y.fa.gz.fai'), path('no_Y.fa.gz.gzi'), 
          path ("chrMT.fa.gz"), path('chrMT.fa.gz.fai'), path('chrMT.fa.gz.gzi'), 
          path ("male_haploid.fa.gz"), path('male_haploid.fa.gz.fai'), path('male_haploid.fa.gz.gzi'),
          path ("autosomes.fa.gz"), path("autosomes.fa.gz.fai"), path("autosomes.fa.gz.gzi"), emit: fastas

    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ploidy = task.ext.ploidy ?: "--diploid"
    def chrX = params.chrX ?: 'chrX'
    def chrY = params.chrY ?: 'chrY'
    def chrMT = params.chrMT ?: 'chrMT'
    def chrs_to_sep_out = "${chrX} ${chrY} ${chrMT}"
    def bcftools_cpus = Math.round((task.cpus-1)/3)
    """
    if [ -f ${params.diploid_haploid_ref_folder}/autosomes.fa.gz ]
    then
        ln -s ${params.diploid_haploid_ref_folder}/*.*.gz.* .
    
    else
        split_out_chrs.py -i $ref_fasta --chrs $chrX $chrY $chrMT

        cat others.fa ${chrX}.fa              | bgzip -c -@ ${task.cpus} > no_Y.fa.gz & \\
        cat ${chrMT}.fa                       | bgzip -c -@ ${task.cpus} > chrMT.fa.gz & \\
        cat others.fa                         | bgzip -c -@ ${task.cpus} > autosomes.fa.gz & \\
        cat ${chrX}.fa ${chrY}.fa ${chrMT}.fa | bgzip -c -@ ${task.cpus} > male_haploid.fa.gz & \\
        wait

        samtools faidx no_Y.fa.gz & \\
        samtools faidx chrMT.fa.gz & \\
        samtools faidx autosomes.fa.gz & \\
        samtools faidx male_haploid.fa.gz & \\
        wait
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1)
    END_VERSIONS
    """
}
