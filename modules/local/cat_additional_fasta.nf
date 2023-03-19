process CAT_ADDITIONAL_FASTA {
    tag "${meta.id}"

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path (fasta), path (gtf), path (add_fasta)
    val (biotype)

    output:
    tuple val(meta), path ("*.fasta"), emit: fasta
    tuple val(meta), path ("*.gtf")  , emit: gtf
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome_name  = params.genome ? params.genome : fasta.getBaseName()
    def biotype_name = biotype ? "-b $biotype" : ''
    def add_name     = add_fasta.getBaseName()
    def name             = "${genome_name}_${add_name}"
    """
    fasta2gtf.py \\
        -o ${add_fasta.baseName}.gtf \\
        $biotype_name \\
        $add_fasta

    cat $fasta $add_fasta > ${name}.fasta
    cat $gtf ${add_fasta.baseName}.gtf > ${name}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}