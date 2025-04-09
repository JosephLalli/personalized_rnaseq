process SEPARATE_AUTOSOMAL_CHROMS {
    tag "sex_specific_references"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'community.wave.seqera.io/library/bcftools_bedtools_samtools:750401edc217be14' }"

    input:
    tuple val (meta), path (vcf)
    tuple val (meta2), path (vcf_index)
    tuple val (meta3), path (ref_fasta)
    tuple val (meta4), path (ref_fasta_fai)
    tuple val (meta5), path (ref_gtf)
    tuple val (meta6), path (ref_par_bedfile)
    val reference_contigs
    val par_regions

    output:
    tuple path ('female_diploid.fa.gz'), path ('female_diploid.fa.gz.fai'), path ('female_diploid.fa.gz.gzi'),
          path ("female_diploid.vcf.gz"), path('female_diploid.vcf.gz.tbi'),
          path ("female_diploid.gtf"), emit: diploid_female_refs, optional: true
    tuple path ('female_haploid.fa.gz'), path ('female_haploid.fa.gz.fai'), path ('female_haploid.fa.gz.gzi'), path ("female_haploid.vcf.gz"), path('female_haploid.vcf.gz.tbi'), path ("female_haploid.gtf"), emit: haploid_female_refs, optional: true
    tuple path ('male_diploid.fa.gz'), path ('male_diploid.fa.gz.fai'), path ('male_diploid.fa.gz.gzi'), path ("male_diploid.vcf.gz"), path('male_diploid.vcf.gz.tbi'), path ("male_diploid.gtf"), emit: diploid_male_refs, optional: true
    tuple path ('male_haploid.fa.gz'), path ('male_haploid.fa.gz.fai'), path ('male_haploid.fa.gz.gzi'), path ("male_haploid.vcf.gz"), path('male_haploid.vcf.gz.tbi'), path ("male_haploid.gtf"), emit: haploid_male_refs, optional: true
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Original chromosome definitions
    def chrX = params.chrX ?: 'chrX'
    def chrY = params.chrY ?: 'chrY'
    def chrM = params.chrMT ?: 'chrM'

    // Define potential haploid chromosomes
    def potential_male_hap_chrs = [chrX, chrY, chrM]
    def potential_female_hap_chrs = [chrM]  // chrX is haploid in males but diploid in females

    def bcftools_male_haploid_region = []
    def bcftools_male_diploid_region = []

    // Filter to only include chromosomes present in the reference
    def female_haploid_chrs = potential_female_hap_chrs.intersect(reference_contigs)
    def bcftools_female_haploid_region = female_haploid_chrs.join(',')
    def male_haploid_chrs = potential_male_hap_chrs.intersect(reference_contigs)
    if (reference_contigs.contains('chrX')){
        male_haploid_chrs = (male_haploid_chrs+[chrX]).unique()
        def PAR1_end = par_regions[0].split(':')[1].split('-')[1]
        def PAR2_start = par_regions[1].split(':')[1].split('-')[0]
        def haploid_chrX = "${chrX}:${PAR1_end}-${PAR2_start}"
        bcftools_male_haploid_region = ((male_haploid_chrs - [chrX]) + haploid_chrX).join(',')
        // bcftools_male_haploid_region = bcftools_male_haploid_region.addAll(male_haploid_chrs)
        // bcftools_male_haploid_region = bcftools_male_haploid_region.remove(chrX)
        // bcftools_male_haploid_region = bcftools_male_haploid_region.addAll(par_regions.collect { "^${it}" })
        // bcftools_male_haploid_region = bcftools_male_haploid_region.join(',')
    }

    // Define diploid chromosomes
    def female_diploid_chrs = reference_contigs.findAll { it2 -> !(it2 in female_haploid_chrs) }
    def bcftools_female_diploid_region = female_diploid_chrs.join(',')
    def male_diploid_chrs = reference_contigs.findAll { it2 -> !(it2 in male_haploid_chrs) }
    if (reference_contigs.contains('chrX')){
        male_diploid_chrs = (male_diploid_chrs+[chrX]).unique()
        bcftools_male_diploid_region = ((male_diploid_chrs - [chrX]) + par_regions).join(',')
        // bcftools_male_diploid_region = bcftools_male_diploid_region.addAll(male_diploid_chrs)
        // bcftools_male_diploid_region = bcftools_male_diploid_region.remove(chrX)
        // bcftools_male_diploid_region = bcftools_male_diploid_region.addAll(par_regions.collect { "${it}" })
        // bcftools_male_diploid_region = bcftools_male_diploid_region.join(',')
    }

    def bcftools_cpus = Math.round((task.cpus-1)/3)
    def diploid_haploid_ref_folder = params.diploid_haploid_ref_folder ?: ""
    """
    if [ -d "${diploid_haploid_ref_folder}" ]
    then
        ln -s ${diploid_haploid_ref_folder}/*.gz* .
        ln -s ${diploid_haploid_ref_folder}/*.gtf .
    else
        # Process female diploid chromosomes
        if [ ! -z "${female_diploid_chrs.join(' ')}" ]; then
            samtools faidx $ref_fasta ${female_diploid_chrs.join(' ')} | bgzip -c -@ ${task.cpus} > female_diploid.fa.gz && \\
                samtools faidx female_diploid.fa.gz &
            cat ${ref_gtf} | grep -Ev "(^$chrY|^$chrM)" > female_diploid.gtf &
            bcftools view --threads ${bcftools_cpus} -O z -t ${bcftools_female_diploid_region} $vcf > female_diploid.vcf.gz && \\
                bcftools index -t female_diploid.vcf.gz &
        fi

        # Process female haploid chromosomes
        if [ ! -z "${female_haploid_chrs.join(' ')}" ]; then
            samtools faidx $ref_fasta ${female_haploid_chrs.join(' ')} | bgzip -c -@ ${task.cpus} > female_haploid.fa.gz && \\
                samtools faidx female_haploid.fa.gz &
            cat $ref_gtf | grep -E "^$chrM" > female_haploid.gtf &
            bcftools view --threads ${bcftools_cpus} -O z -t ${bcftools_female_haploid_region} $vcf > female_haploid.vcf.gz  && \\
                bcftools index -t female_haploid.vcf.gz &
        fi

        # Process male diploid chromosomes
        if [ ! -z "${male_diploid_chrs.join(' ')}" ]; then
            samtools faidx $ref_fasta ${male_diploid_chrs.join(' ')} | bgzip -c -@ ${task.cpus} > male_diploid.fa.gz && \\
                samtools faidx male_diploid.fa.gz &
            grep -Ev "(^$chrX|^$chrY|^$chrM)" $ref_gtf > male_diploid.gtf && \\
            bedtools intersect -a $ref_gtf -b $ref_par_bedfile -v | grep "^$chrX" >> male_diploid.gtf &
            bcftools view --threads ${bcftools_cpus} -O z -t ${bcftools_male_diploid_region} $vcf > male_diploid.vcf.gz && \\
                bcftools index -t male_diploid.vcf.gz &
        fi

        # Process male haploid chromosomes
        if [ ! -z "${male_haploid_chrs.join(' ')}" ]; then
            samtools faidx $ref_fasta ${male_haploid_chrs.join(' ')} | bgzip -c > male_haploid.fa.gz && \\
                samtools faidx male_haploid.fa.gz &
            bedtools intersect -a $ref_gtf -b $ref_par_bedfile > male_haploid.gtf && \\
            grep -E "(^$chrY|^$chrM)" $ref_gtf >> male_haploid.gtf &
            bcftools view --threads ${bcftools_cpus} -O z -t ${bcftools_male_haploid_region} $vcf > male_haploid.vcf.gz && \\
                bcftools index -t male_haploid.vcf.gz &
        fi
        wait
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n 1 | cut -d' ' -f2)
        htslib: \$(bcftools --version 2>&1 | head -n 2 | tail -n 1 | cut -d' ' -f3)
    END_VERSIONS
    """

    stub:
    """
    touch male_diploids.fa.gz
    touch male_diploids.fa.gz.fai
    touch male_diploids.fa.gz.gzi
    touch male_diploids.vcf.gz
    touch male_diploids.vcf.gz.tbi
    touch no_Y.fa.gz
    touch no_Y.fa.gz.fai
    touch no_Y.fa.gz.gzi
    touch no_Y.vcf.gz
    touch no_Y.vcf.gz.tbi
    touch chrM.fa.gz
    touch chrM.fa.gz.fai
    touch chrM.fa.gz.gzi
    touch chrM.vcf.gz
    touch chrM.vcf.gz.tbi
    touch female_diploid.gtf
    touch female_haploid.gtf
    touch male_haploid.fa.gz
    touch male_haploid.fa.gz.fai
    touch male_haploid.fa.gz.gzi
    touch male_haploid.vcf.gz
    touch male_haploid.vcf.gz.tbi
    touch male_haploid.gtf
    touch male_diploid.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n 1 | cut -d' ' -f2)
        htslib: \$(bcftools --version 2>&1 | head -n 2 | tail -n 1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
