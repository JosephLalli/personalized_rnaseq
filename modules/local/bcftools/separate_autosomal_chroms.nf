process SEPARATE_AUTOSOMAL_CHROMS {
    tag "sex_specific_references"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'jlalli/g2gtools:0.2.9' }"

    input:
    path (vcf)
    path (vcf_index)
    path (ref_fasta)
    path (ref_fasta_fai)
    path (ref_gtf)

    output:
    tuple path ('no_Y.fa.gz'), path ('no_Y.fa.gz.fai'), path ('no_Y.fa.gz.gzi'), 
          path ("no_Y.vcf.gz"), path('no_Y.vcf.gz.tbi'),
          path ("female_diploid.gtf"), emit: diploid_female_refs
    tuple path ('chrMT.fa.gz'), path ('chrMT.fa.gz.fai'), path ('chrMT.fa.gz.gzi'), path ("chrMT.vcf.gz"), path('chrMT.vcf.gz.tbi'), path ("female_haploid.gtf"), emit: haploid_female_refs
    tuple path ('autosomes.fa.gz'), path ('autosomes.fa.gz.fai'), path ('autosomes.fa.gz.gzi'), path ("autosomes.vcf.gz"), path('autosomes.vcf.gz.tbi'), path ("male_diploid.gtf"), emit: diploid_male_refs
    tuple path ('male_haploid.fa.gz'), path ('male_haploid.fa.gz.fai'), path ('male_haploid.fa.gz.gzi'), path ("male_haploid.vcf.gz"), path('male_haploid.vcf.gz.tbi'), path ("male_haploid.gtf"), emit: haploid_male_refs
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def chrX = params.chrX ?: 'chrX'
    def chrY = params.chrY ?: 'chrY'
    def chrMT = params.chrMT ?: 'chrMT'
    def chrs_to_sep_out = "${chrX} ${chrY} ${chrMT}"
    def bcftools_cpus = Math.round((task.cpus-1)/3)
    def diploid_haploid_ref_folder = params.diploid_haploid_ref_folder
    // println(ref_fasta_fai.getClass())
    // def all_chrs = file(ref_fasta_fai.toString()).readLines().each{ it.split{'\t'}[0] }
    // // println(all_chrs)
    // def female_diploid_chrs = 'banana' //all_chrs.remove(chrY).remove(chrMT)
    // def male_diploid_chrs = 'banana' //female_diploid_chrs.remove(chrX)
    // println(male_diploid_chrs)
    // 
    // 
    """
    if [ -f ${diploid_haploid_ref_folder}/autosomes.fa.gz ]
    then
        ln -s ${diploid_haploid_ref_folder}/*.gz* .
        ln -s ${diploid_haploid_ref_folder}/*.gtf .
    else

        samtools faidx $ref_fasta \$(cut -f 1 < $ref_fasta_fai | tr '\\n' ' ' | sed 's/$chrX//' | sed 's/$chrY//' | sed 's/$chrMT//') | bgzip -c -@ ${task.cpus} > autosomes.fa.gz & \\
        samtools faidx $ref_fasta \$(cut -f 1 < $ref_fasta_fai | tr '\\n' ' ' | sed 's/$chrY//' | sed 's/$chrMT//') | bgzip -c -@ ${task.cpus} > no_Y.fa.gz & \\
        samtools faidx $ref_fasta $chrMT | bgzip -c > chrMT.fa.gz & \\
        samtools faidx $ref_fasta $chrs_to_sep_out | bgzip -c > male_haploid.fa.gz & \\
        wait

        cat ${ref_gtf} | grep -Ev "(^$chrY|^$chrMT)"        > female_diploid.gtf & \\
        cat ${ref_gtf} | grep -E "^$chrMT"                  > female_haploid.gtf & \\
        cat ${ref_gtf} | grep -Ev "(^$chrX|^$chrY|^$chrMT)" > male_diploid.gtf & \\
        cat ${ref_gtf} | grep -E  "(^$chrX|^$chrY|^$chrMT)" > male_haploid_unaltered.gtf & \\
        \\
        bcftools view --threads ${bcftools_cpus} -O z -t ^$chrY,$chrMT $vcf > no_Y.vcf.gz  & \\
        bcftools view --threads ${bcftools_cpus} -O z -t $chrMT $vcf > chrMT.vcf.gz  & \\
        bcftools view --threads ${bcftools_cpus} -O z -t $chrX,$chrY,$chrMT $vcf > male_haploid.vcf.gz & \\
        bcftools view --threads ${bcftools_cpus} -O z -t ^$chrX,$chrY,$chrMT $vcf > autosomes.vcf.gz  & \\
        wait
        
        append_suffixes_to_sex_chroms.py -i male_haploid_unaltered.gtf -o male_haploid.gtf -X $chrX -Y $chrY & \\
        samtools faidx no_Y.fa.gz & \\
        samtools faidx chrMT.fa.gz & \\
        samtools faidx autosomes.fa.gz & \\
        samtools faidx male_haploid.fa.gz & \\
        bcftools index -t no_Y.vcf.gz & \\
        bcftools index -t chrMT.vcf.gz & \\
        bcftools index -t male_haploid.vcf.gz & \\
        bcftools index -t autosomes.vcf.gz & \\
        wait
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1)
    END_VERSIONS
    """
}
