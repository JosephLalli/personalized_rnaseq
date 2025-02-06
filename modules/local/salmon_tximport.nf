process SALMON_TXIMPORT {
    label "process_medium"

    conda (params.enable_conda ? "bioconda::bioconductor-tximeta=1.8.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta:1.8.0--r40_0' :
        'quay.io/biocontainers/bioconductor-tximeta:1.8.0--r40_0' }"

    input:
    path ("salmon/*")
    path (tx2gene)

    output:
    path "*gene_tpm.tsv"                   , emit: tpm_gene
    path "*gene_counts.tsv"                , emit: counts_gene
    path "*gene_length.tsv"                , emit: lengths_gene
    path "*gene_tpm_length_scaled.tsv"     , emit: tpm_gene_length_scaled
    path "*gene_counts_length_scaled.tsv"  , emit: counts_gene_length_scaled
    path "*gene_length_length_scaled.tsv"  , emit: lengths_gene_length_scaled
    path "*gene_tpm_scaled.tsv"            , emit: tpm_gene_scaled
    path "*gene_counts_scaled.tsv"         , emit: counts_gene_scaled
    path "*gene_length_scaled.tsv"         , emit: lengths_gene_scaled
    path "*transcript_tpm.tsv"             , emit: tpm_transcript
    path "*transcript_counts.tsv"          , emit: counts_transcript
    path "*transcript_length.tsv"          , emit: length_transcript
    path "*transcript_tpm_dtuscaled.tsv"   , emit: tpm_transcript_dtu
    path "*transcript_counts_dtuscaled.tsv", emit: counts_transcript_dtu
    path "*transcript_length_dtuscaled.tsv", emit: length_transcript_dtu
    path "*gene_tpm_length_scaled_median_bootstraps.tsv", emit: tpm_gene_length_scaled_median_bootstraps
    path "*gene_counts_length_scaled_median_bootstraps.tsv", emit: counts_gene_length_scaled_median_bootstraps
    path "*gene_length_length_scaled_median_bootstraps.tsv", emit: length_gene_length_scaled_median_bootstraps
    path "*transcript_tpm_length_scaled_median_bootstraps.tsv", emit: tpm_transcript_length_scaled_median_bootstraps
    path "*transcript_counts_length_scaled_median_bootstraps.tsv", emit: counts_transcript_length_scaled_median_bootstraps
    path "*transcript_length_length_scaled_median_bootstraps.tsv", emit: length_transcript_length_scaled_median_bootstraps

    path "*transcript_tpm.tximport_se.rds", emit: merged_transcript_rds
    path "*transcript_tpm_dtuscaled.tximport_se.rds", emit: merged_transcript_dtuscaled_rds
    path "*transcript_tpm_length_scaled.median_bootstrap.tximport_se.rds", emit: merged_transcript_length_scaled_median_bootstrap_rds
    path "*gene_tpm.tximport_se.rds", emit: merged_gene_rds
    path "*gene_tpm_length_scaled.tximport_se.rds", emit: merged_gene_rds_length_scaled
    path "*gene_tpm_scaled.tximport_se.rds", emit: merged_gene_rds_scaled
    path "*gene_tpm_length_scaled.median_bootstrap.tximport_se.rds", emit: merged_gene_length_scaled_median_bootstrap_rds


    path "*.tsv", emit: tsv
    path "*se.rds", emit: rds

    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_tximport.r \\
        NULL \\
        salmon \\
        salmon.merged

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """
}
