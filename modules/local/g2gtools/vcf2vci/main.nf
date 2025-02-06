process VCF2VCI {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'jlalli/g2gtools:0.2.9' }"

    input:
    tuple val (meta), path (ref_fasta), path (ref_fai), path(ref_gzi), path (vcf), path (vcf_tbi), path(gtf)

    output:
    tuple val(meta), path("*.vci.gz"), path("*.vci.gz.tbi"), emit: vci
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ploidy_arg = meta.ploidy == 'diploid' ? "--diploid" : ""
    def previously_generated_file_path = params.force_resume ? task.publishDir.path[0] : ""
    """
    if [[ -s ${previously_generated_file_path}/${prefix}.vci.gz ]]
    then
        ln -s ${previously_generated_file_path}/${prefix}.vci.gz .
        ln -s ${previously_generated_file_path}/${prefix}.vci.gz.tbi .
    else
        g2gtools vcf2vci -o ${prefix}.vci \\
                        -s ${meta.sample} \\
                        -f ${ref_fasta} \\
                        ${ploidy_arg} \\
                        -p ${task.cpus}\\
                        ${args} \\
                        -i ${vcf}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        g2gtools: \$(g2gtools --version 2>&1)
    END_VERSIONS
    """
}
// (cat HSB100-haploid.vci | cut -f 1 | uniq | grep ^[A-Za-z] | sort ) > in_vci.txt
// (cat HSB100-haploid.vci | cut -f 1 | uniq | grep "##CONTIG" | sed 's/##CONTIG=//' | sed 's/:.*$//' | sort ) > supposed_to_be_in_vci.txt

// for contig in \$(comm -13 <(cat in_vci.txt) <(cat supposed_to_be_in_vci.txt))
// do
//     echo "$contig   1   .   .   .   1" >> HSB100-haploid.vci
// done