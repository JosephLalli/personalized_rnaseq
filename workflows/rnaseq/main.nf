/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { DESEQ2_QC as DESEQ2_QC_STAR_SALMON } from '../../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_RSEM        } from '../../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_PSEUDO      } from '../../modules/local/deseq2_qc'
include { MULTIQC_CUSTOM_BIOTYPE             } from '../../modules/local/multiqc_custom_biotype'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { ALIGN_STAR                            } from '../../subworkflows/local/align_star'
include { QUANTIFY_RSEM                         } from '../../subworkflows/local/quantify_rsem'
include { BAM_DEDUP_UMI as BAM_DEDUP_UMI_STAR   } from '../../subworkflows/nf-core/bam_dedup_umi'
include { BAM_DEDUP_UMI as BAM_DEDUP_UMI_HISAT2 } from '../../subworkflows/nf-core/bam_dedup_umi'

include { checkSamplesAfterGrouping      } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getStarPercentMapped           } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { biotypeInGtf                   } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getInferexperimentStrandedness } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { methodsDescriptionText         } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { checkMaxContigSize             } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PREPARE_GENOME as PREPARE_PERSONALIZED_GENOMES } from '../../subworkflows/local/prepare_genome'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { DUPRADAR                   } from '../../modules/nf-core/dupradar'
include { PRESEQ_LCEXTRAP            } from '../../modules/nf-core/preseq/lcextrap'
include { QUALIMAP_RNASEQ            } from '../../modules/nf-core/qualimap/rnaseq'
include { STRINGTIE_STRINGTIE        } from '../../modules/nf-core/stringtie/stringtie'
include { SUBREAD_FEATURECOUNTS      } from '../../modules/nf-core/subread/featurecounts'
include { KRAKEN2_KRAKEN2 as KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN } from '../../modules/nf-core/bracken/bracken/main'
include { MULTIQC                    } from '../../modules/nf-core/multiqc'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_FW          } from '../../modules/nf-core/bedtools/genomecov'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_REV         } from '../../modules/nf-core/bedtools/genomecov'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { samplesheetToList                } from 'plugin/nf-schema'
include { multiqcTsvFromList               } from '../../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
include { paramsSummaryMultiqc             } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQ_ALIGN_HISAT2               } from '../../subworkflows/nf-core/fastq_align_hisat2'
include { BAM_MARKDUPLICATES_PICARD        } from '../../subworkflows/nf-core/bam_markduplicates_picard'
include { BAM_RSEQC                        } from '../../subworkflows/nf-core/bam_rseqc'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { QUANTIFY_PSEUDO_ALIGNMENT as QUANTIFY_STAR_SALMON } from '../../subworkflows/nf-core/quantify_pseudo_alignment'
include { QUANTIFY_PSEUDO_ALIGNMENT                         } from '../../subworkflows/nf-core/quantify_pseudo_alignment'
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS              } from '../../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {

    take:
    params
    ch_samplesheet       // channel: [ meta, path/to/samplesheet.csv ]
    ch_fastq             // channel: [ meta, path/to/reads ]
    ch_versions          // channel: [ path(versions.yml) ]
    ch_fasta             // channel: [ meta, path(genome.fasta) ]
    ch_gtf               // channel: [ meta, path(genome.gtf) ]
    ch_fai               // channel: [ meta, path(genome.fai) ]
    ch_chrom_sizes       // channel: [ meta, path(genome.sizes) ]
    ch_gene_bed          // channel: [ meta, path(gene.bed) ]
    ch_transcript_fasta  // channel: [ meta, path(transcript.fasta) ]
    ch_star_index        // channel: [ meta, path(star/index/) ]
    ch_rsem_index        // channel: [ meta, path(rsem/index/) ]
    ch_hisat2_index      // channel: [ meta, path(hisat2/index/) ]
    ch_salmon_index      // channel: [ meta, path(salmon/index/) ]
    ch_kallisto_index    // channel: [ meta, path(kallisto/index/) ] ]
    ch_bbsplit_index     // channel: [ meta, path(bbsplit/index/) ]
    ch_ribo_db           // channel: [ meta, path(sortmerna_fasta_list) ]
    ch_sortmerna_index   // channel: [ meta, path(sortmerna/index/) ]
    ch_splicesites       // channel: [ meta, path(genome.splicesites.txt) ]

    main:

    ch_multiqc_files = Channel.empty()
    ch_trim_status = Channel.empty()
    ch_map_status = Channel.empty()
    ch_strand_status = Channel.empty()

    // Header files for MultiQC
    ch_pca_header_multiqc           = file("$projectDir/workflows/rnaseq/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
    sample_status_header_multiqc    = file("$projectDir/workflows/rnaseq/assets/multiqc/sample_status_header.txt", checkIfExists: true)
    ch_clustering_header_multiqc    = file("$projectDir/workflows/rnaseq/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
    ch_biotypes_header_multiqc      = file("$projectDir/workflows/rnaseq/assets/multiqc/biotypes_header.txt", checkIfExists: true)
    ch_dummy_file                   = Channel.value(ch_pca_header_multiqc)

    //
    // Run RNA-seq FASTQ preprocessing subworkflow
    //

    // The subworkflow only has to do Salmon indexing if it discovers 'auto'
    // samples, and if we haven't already made one elsewhere
    salmon_index_available = params.salmon_index || (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon')

    // add sample to meta to ensure compatibility with downstream practices
    // ch_fastq = ch_fastq.map{meta, reads -> [meta + ['sample':meta.id], reads] }

    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
        ch_fastq,
        ch_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        ch_sortmerna_index,
        ch_bbsplit_index,
        ch_ribo_db,
        params.skip_bbsplit,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming,
        params.skip_umi_extract,
        !salmon_index_available,
        !params.sortmerna_index && params.remove_ribo_rna,
        params.trimmer,
        params.min_trimmed_reads,
        params.save_trimmed,
        params.remove_ribo_rna,
        params.with_umi,
        params.umi_discard_read,
        params.stranded_threshold,
        params.unstranded_threshold,
        params.skip_linting,
        false
    )

    ch_multiqc_files                  = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)
    ch_versions                       = ch_versions.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.versions)
    ch_strand_inferred_filtered_fastq = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads
    ch_trim_read_count                = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_read_count

    if ( workflow.stubRun ) {
        params.min_trimmed_reads = -1
    }

    ch_trim_status = ch_trim_read_count
        .map {
            meta, num_reads ->
                return [ meta.id, num_reads >= params.min_trimmed_reads.toFloat() ]
        }


    ch_meta = ch_strand_inferred_filtered_fastq.map{ it -> it[0] + ['sample': it[0].id]}
    ch_samplesheet = ch_meta.combine(ch_samplesheet)
    ch_dummy_file = ch_meta.combine(ch_dummy_file)
    if (params.use_personalized_references) {
        PREPARE_PERSONALIZED_GENOMES(
            params,
            ch_meta,
            params.fasta,
            params.fai,
            params.gtf,
            params.gff,
            params.additional_fasta,
            params.transcript_fasta,
            params.vcf,
            params.vcf_index,
            params.gene_bed,
            params.par_bed,
            params.splicesites,
            params.bbsplit_fasta_list,
            params.ribo_database_manifest,
            params.star_index,
            params.rsem_index,
            params.salmon_index,
            params.kallisto_index,
            params.hisat2_index,
            params.bbsplit_index,
            params.sortmerna_index,
            params.gencode,
            params.featurecounts_group_type,
            params.aligner,
            params.pseudo_aligner,
            params.skip_gtf_filter,
            params.skip_bbsplit,
            !params.remove_ribo_rna,
            params.skip_alignment,
            params.skip_pseudo_alignment
        )

        ch_fasta = PREPARE_PERSONALIZED_GENOMES.out.fasta
        ch_fai = PREPARE_PERSONALIZED_GENOMES.out.fai
        ch_gtf = PREPARE_PERSONALIZED_GENOMES.out.gtf
        ch_chrom_sizes = PREPARE_PERSONALIZED_GENOMES.out.chrom_sizes
        ch_gene_bed = PREPARE_PERSONALIZED_GENOMES.out.gene_bed
        ch_transcript_fasta = PREPARE_PERSONALIZED_GENOMES.out.transcript_fasta
        ch_star_index = PREPARE_PERSONALIZED_GENOMES.out.star_index
        ch_hisat2_index = PREPARE_PERSONALIZED_GENOMES.out.hisat2_index
        ch_salmon_index = PREPARE_PERSONALIZED_GENOMES.out.salmon_index
        ch_kallisto_index = PREPARE_PERSONALIZED_GENOMES.out.kallisto_index
        ch_bbsplit_index = PREPARE_PERSONALIZED_GENOMES.out.bbsplit_index
        ch_ribo_db = PREPARE_PERSONALIZED_GENOMES.out.rrna_fastas
        ch_sortmerna_index = PREPARE_PERSONALIZED_GENOMES.out.sortmerna_index
        ch_splicesites = PREPARE_PERSONALIZED_GENOMES.out.splicesites

        ch_versions = PREPARE_PERSONALIZED_GENOMES.out.versions

        if (!params.skip_alignment && !params.bam_csi_index) {
            PREPARE_PERSONALIZED_GENOMES
                .out
                .fai
                .map { _meta, fai -> checkMaxContigSize(fai) }
        }

    } else {
        ch_fasta = ch_meta.combine(ch_fasta.map{ it -> it[1] })
        ch_fai = ch_meta.combine(ch_fai.map{ it -> it[1] })
        ch_gtf = ch_meta.combine(ch_gtf.map{ it -> it[1] })
        ch_chrom_sizes = ch_meta.combine(ch_chrom_sizes.map{ it -> it[1] })
        ch_gene_bed = ch_meta.combine(ch_gene_bed.map{ it -> it[1] })
        ch_transcript_fasta = ch_meta.combine(ch_transcript_fasta.map{ it -> it[1] })
        ch_star_index = ch_meta.combine(ch_star_index.map{ it -> it[1] })
        ch_hisat2_index = ch_meta.combine(ch_hisat2_index.map{ it -> it[1] })
        ch_salmon_index = ch_meta.combine(ch_salmon_index.map{ it -> it[1] })
        ch_kallisto_index = ch_meta.combine(ch_kallisto_index.map{ it -> it[1] })
        ch_bbsplit_index = ch_meta.combine(ch_bbsplit_index.map{ it -> it[1] })
        ch_ribo_db = ch_meta.combine(ch_ribo_db.map{ it -> it[1] })
        ch_sortmerna_index = ch_meta.combine(ch_sortmerna_index.map{ it -> it[1] })
        ch_splicesites = ch_meta.combine(ch_splicesites.map{ it -> it[1] })
    }

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //
    ch_genome_bam          = Channel.empty()
    ch_genome_bam_index    = Channel.empty()
    ch_star_log            = Channel.empty()
    ch_unaligned_sequences = Channel.empty()
    ch_transcriptome_bam   = Channel.empty()

    if (!params.skip_alignment && params.aligner == 'star_salmon') {
        // Check if an AWS iGenome has been provided to use the appropriate version of STAR
        def is_aws_igenome = false
        if (params.fasta && params.gtf) {
            if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
                is_aws_igenome = true
            }
        }

        ALIGN_STAR (
            ch_strand_inferred_filtered_fastq,
            ch_star_index,
            ch_gtf,
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            is_aws_igenome,
            ch_fasta
        )
        ch_genome_bam          = ALIGN_STAR.out.bam
        ch_genome_bam_index    = ALIGN_STAR.out.bai
        ch_transcriptome_bam   = ALIGN_STAR.out.bam_transcript
        ch_star_log            = ALIGN_STAR.out.log_final
        ch_unaligned_sequences = ALIGN_STAR.out.fastq
        ch_multiqc_files = ch_multiqc_files.mix(ch_star_log.collect{ it -> it[1] })

        if (params.bam_csi_index) {
            ch_genome_bam_index = ALIGN_STAR.out.csi
        }
        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)

        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //
        if (params.with_umi) {

            BAM_DEDUP_UMI_STAR(
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                ch_fasta,
                params.umi_dedup_tool,
                params.umitools_dedup_stats,
                params.bam_csi_index,
                ch_transcriptome_bam,
                ch_transcript_fasta
            )

            ch_genome_bam        = BAM_DEDUP_UMI_STAR.out.bam
            ch_transcriptome_bam = BAM_DEDUP_UMI_STAR.out.transcriptome_bam
            ch_genome_bam_index  = BAM_DEDUP_UMI_STAR.out.bai
            ch_versions          = ch_versions.mix(BAM_DEDUP_UMI_STAR.out.versions)

            ch_multiqc_files = ch_multiqc_files
                .mix(BAM_DEDUP_UMI_STAR.out.multiqc_files)

        } else {
            // The deduplicated stats should take priority for MultiQC, but use
            // them straight out of the aligner otherwise

            ch_multiqc_files = ch_multiqc_files
                .mix(ALIGN_STAR.out.stats.collect{ it -> it[1] })
                .mix(ALIGN_STAR.out.flagstat.collect{ it -> it[1] })
                .mix(ALIGN_STAR.out.idxstats.collect{ it -> it[1] })
        }

        //
        // SUBWORKFLOW: Count reads from BAM alignments using Salmon
        //
        QUANTIFY_STAR_SALMON (
            ch_samplesheet,
            ch_transcriptome_bam,
            ch_dummy_file,
            ch_transcript_fasta,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            'salmon',
            true,
            '',
            // params.salmon_quant_libtype ?: '',
            params.kallisto_quant_fraglen,
            params.kallisto_quant_fraglen_sd
        )
        ch_versions = ch_versions.mix(QUANTIFY_STAR_SALMON.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_STAR_SALMON (
                QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled.map { it -> it[1] },
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_STAR_SALMON.out.pca_multiqc.collect())
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_STAR_SALMON.out.dists_multiqc.collect())
            ch_versions = ch_versions.mix(DESEQ2_QC_STAR_SALMON.out.versions)
        }
    }

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with RSEM
    //
    if (!params.skip_alignment && params.aligner == 'star_rsem') {
        QUANTIFY_RSEM (
            ch_strand_inferred_filtered_fastq,
            ch_rsem_index,
            ch_fasta
        )
        ch_genome_bam       = QUANTIFY_RSEM.out.bam
        ch_genome_bam_index = QUANTIFY_RSEM.out.bai
        ch_star_log         = QUANTIFY_RSEM.out.logs
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.stats.collect{ it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.flagstat.collect{ it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.idxstats.collect{ it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(ch_star_log.collect{ it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_RSEM.out.stat.collect{ it -> it[1] })

        if (params.bam_csi_index) {
            ch_genome_bam_index = QUANTIFY_RSEM.out.csi
        }
        ch_versions = ch_versions.mix(QUANTIFY_RSEM.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_RSEM (
                QUANTIFY_RSEM.out.merged_counts_gene,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_RSEM.out.pca_multiqc.collect())
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_RSEM.out.dists_multiqc.collect())
            ch_versions = ch_versions.mix(DESEQ2_QC_RSEM.out.versions)
        }
    }

    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        FASTQ_ALIGN_HISAT2 (
            ch_strand_inferred_filtered_fastq,
            ch_hisat2_index,
            ch_splicesites,
            ch_fasta
        )
        ch_genome_bam          = FASTQ_ALIGN_HISAT2.out.bam
        ch_genome_bam_index    = FASTQ_ALIGN_HISAT2.out.bai
        ch_unaligned_sequences = FASTQ_ALIGN_HISAT2.out.fastq
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.summary.collect{ it -> it[1] })

        if (params.bam_csi_index) {
            ch_genome_bam_index = FASTQ_ALIGN_HISAT2.out.csi
        }
        ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //

        if (params.with_umi) {

            BAM_DEDUP_UMI_HISAT2(
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                ch_fasta,
                params.umi_dedup_tool,
                params.umitools_dedup_stats,
                params.bam_csi_index,
                ch_transcriptome_bam,
                ch_transcript_fasta
            )

            ch_genome_bam        = BAM_DEDUP_UMI_HISAT2.out.bam
            ch_genome_bam_index  = BAM_DEDUP_UMI_HISAT2.out.bai
            ch_versions          = ch_versions.mix(BAM_DEDUP_UMI_HISAT2.out.versions)

            ch_multiqc_files = ch_multiqc_files
                .mix(BAM_DEDUP_UMI_HISAT2.out.multiqc_files)
        } else {

            // The deduplicated stats should take priority for MultiQC, but use
            // them straight out of the aligner otherwise
            ch_multiqc_files = ch_multiqc_files
                .mix(FASTQ_ALIGN_HISAT2.out.stats.collect{ it -> it[1] })
                .mix(FASTQ_ALIGN_HISAT2.out.flagstat.collect{ it -> it[1] })
                .mix(FASTQ_ALIGN_HISAT2.out.idxstats.collect{ it -> it[1] })
        }
    }

    //
    // Filter channels to get samples that passed STAR minimum mapping percentage
    //
    if (!params.skip_alignment && params.aligner.contains('star')) {
        ch_star_log
            .map { meta, align_log -> [ meta ] + getStarPercentMapped(params, align_log) }
            .set { ch_percent_mapped }

        // Save status for workflow summary
        ch_map_status = ch_percent_mapped
            .map {
                meta, _mapped, pass ->
                    return [ meta.id, pass ]
            }

        ch_genome_bam
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, _mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam }

        ch_genome_bam_index
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, _mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam_index }

        ch_percent_mapped
            .branch { meta, mapped, pass ->
                pass: pass
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    return [ "$meta.id\t$mapped" ]
            }
            .set { ch_pass_fail_mapped }

        ch_pass_fail_mapped
            .fail
            .collect()
            .map {
                tsv_data ->
                    def header = ["Sample", "STAR uniquely mapped reads (%)"]
                    sample_status_header_multiqc.text + multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_mapping_multiqc }
        ch_multiqc_files = ch_multiqc_files.mix(ch_fail_mapping_multiqc.collectFile(name: 'fail_mapped_samples_mqc.tsv'))
    }

    //
    // MODULE: Run Preseq
    //
    if (!params.skip_alignment && !params.skip_qc && !params.skip_preseq) {
        PRESEQ_LCEXTRAP (
            ch_genome_bam
        )
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.collect{ it -> it[1] })
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
    }

    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    if (!params.skip_alignment && !params.skip_markduplicates && !params.with_umi) {
        BAM_MARKDUPLICATES_PICARD (
            ch_genome_bam,
            ch_fasta,
            ch_fai
        )
        ch_genome_bam       = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.bai
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect{ it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{ it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect{ it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect{ it -> it[1] })

        if (params.bam_csi_index) {
            ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.csi
        }
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
    }

    //
    // MODULE: STRINGTIE
    //
    if (!params.skip_alignment && !params.skip_stringtie) {
        STRINGTIE_STRINGTIE (
            ch_genome_bam.join(ch_gtf)
        )
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())
    }

    //
    // MODULE: Feature biotype QC using featureCounts
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    if (!params.skip_alignment && !params.skip_qc && !params.skip_biotype_qc && biotype) {

        ch_gtf
            .map { meta, gtf -> [ meta, biotypeInGtf(gtf, biotype) ]}
            .set { biotype_in_gtf }

        // Prevent any samples from running if GTF file doesn't have a valid biotype
        ch_genome_bam
            .join(ch_gtf)
            .join(biotype_in_gtf)
            .filter { it -> it[-1] }
            .map { it -> it[0..<it.size()-1] }
            .set { ch_featurecounts }

        SUBREAD_FEATURECOUNTS (
            ch_featurecounts
        )
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

        MULTIQC_CUSTOM_BIOTYPE (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_biotypes_header_multiqc
        )
        ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_BIOTYPE.out.tsv.collect{ it -> it[1] })
        ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())
    }

    //
    // MODULE: Genome-wide coverage with BEDTools
    //
    if (!params.skip_alignment && !params.skip_bigwig) {

        ch_genomecov_input = ch_genome_bam.map { meta, bam -> [ meta, bam, 1 ] }

        BEDTOOLS_GENOMECOV_FW (
            ch_genomecov_input,
            [],
            'bedGraph',
            true
        )
        BEDTOOLS_GENOMECOV_REV (
            ch_genomecov_input,
            [],
            'bedGraph',
            true
        )

        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_FW.out.versions.first())

        //
        // SUBWORKFLOW: Convert bedGraph to bigWig
        //
        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD (
            BEDTOOLS_GENOMECOV_FW.out.genomecov,
            ch_chrom_sizes
        )
        ch_versions = ch_versions.mix(BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD.out.versions)

        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE (
            BEDTOOLS_GENOMECOV_REV.out.genomecov,
            ch_chrom_sizes
        )
    }

    //
    // MODULE: Downstream QC steps
    //
    if (!params.skip_alignment && !params.skip_qc) {
        if (!params.skip_qualimap) {
            QUALIMAP_RNASEQ (
                ch_genome_bam.join(ch_gtf)
            )
            ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_RNASEQ.out.results.collect{ it -> it[1] })
            ch_versions = ch_versions.mix(QUALIMAP_RNASEQ.out.versions.first())
        }

        if (!params.skip_dupradar) {
            DUPRADAR (
                ch_genome_bam.join(ch_gtf)
            )
            ch_multiqc_files = ch_multiqc_files.mix(DUPRADAR.out.multiqc.collect{ it -> it[1] })
            ch_versions = ch_versions.mix(DUPRADAR.out.versions.first())
        }

        // Get RSeqC modules to run
        def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
        if (params.bam_csi_index) {
            rseqc_modules = rseqc_modules.findAll { !['read_distribution', 'inner_distance', 'tin'].contains(it) }
        }

        if (!params.skip_rseqc && rseqc_modules.size() > 0) {
            BAM_RSEQC (
                ch_genome_bam.join(ch_genome_bam_index).join(ch_gene_bed),
                rseqc_modules
            )
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.bamstat_txt.collect{ it -> it[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.inferexperiment_txt.collect{ it -> it[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.innerdistance_freq.collect{ it -> it[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.junctionannotation_log.collect{ it -> it[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.junctionsaturation_rscript.collect{ it -> it[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readdistribution_txt.collect{ it -> it[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readduplication_pos_xls.collect{ it -> it[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.tin_txt.collect{ it -> it[1] })
            ch_versions = ch_versions.mix(BAM_RSEQC.out.versions)

            // Compare predicted supplied or Salmon-predicted strand with what we get from RSeQC
            ch_strand_comparison = BAM_RSEQC.out.inferexperiment_txt
                .map {
                    meta, strand_log ->
                        def rseqc_inferred_strand = getInferexperimentStrandedness(strand_log, params.stranded_threshold, params.unstranded_threshold)
                        def rseqc_strandedness = rseqc_inferred_strand.inferred_strandedness

                        def status = 'fail'
                        def multiqc_lines = []
                        if (meta.salmon_strand_analysis) {
                            def salmon_strandedness = meta.salmon_strand_analysis.inferred_strandedness

                            if (salmon_strandedness == rseqc_strandedness && rseqc_strandedness != 'undetermined') {
                                status = 'pass'
                            }
                            multiqc_lines = [
                                "$meta.id \tSalmon\t$status\tauto\t${meta.salmon_strand_analysis.values().join('\t')}",
                                "$meta.id\tRSeQC\t$status\tauto\t${rseqc_inferred_strand.values().join('\t')}"
                            ]
                        }
                        else {
                            if (meta.strandedness == rseqc_strandedness) {
                                status = 'pass'
                            }
                            multiqc_lines = [ "$meta.id\tRSeQC\t$status\t$meta.strandedness\t${rseqc_inferred_strand.values().join('\t')}" ]
                        }
                        return [ meta, status, multiqc_lines ]
                }
                .multiMap {
                    meta, status, multiqc_lines ->
                        status: [ meta.id, status == 'pass' ]
                        multiqc_lines: multiqc_lines
                }

            // Store the statuses for output
            ch_strand_status = ch_strand_comparison.status

            // Take the lines formatted for MultiQC and output
            ch_strand_comparison.multiqc_lines
                .flatten()
                .collect()
                .map {
                    tsv_data ->
                        def header = [
                            "Sample",
                            "Strand inference method",
                            "Status",
                            "Provided strandedness",
                            "Inferred strandedness",
                            "Sense (%)",
                            "Antisense (%)",
                            "Unstranded (%)"
                        ]
                        sample_status_header_multiqc.text + multiqcTsvFromList(tsv_data, header)
                }
                .set { ch_fail_strand_multiqc }

            ch_multiqc_files = ch_multiqc_files.mix(ch_fail_strand_multiqc.collectFile(name: 'fail_strand_check_mqc.tsv'))
        }

        if (params.contaminant_screening in ['kraken2', 'kraken2_bracken'] ) {
            KRAKEN2 (
                ch_unaligned_sequences,
                params.kraken_db,
                params.save_kraken_assignments,
                params.save_kraken_unassigned
            )
            ch_kraken_reports = KRAKEN2.out.report
            ch_versions = ch_versions.mix(KRAKEN2.out.versions)

            if (params.contaminant_screening == 'kraken2') {
                ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2.out.report.collect{ it -> it[1] })
            } else if (params.contaminant_screening == 'kraken2_bracken') {
                BRACKEN (
                    ch_kraken_reports,
                    params.kraken_db
                )
                ch_versions = ch_versions.mix(BRACKEN.out.versions)
                ch_multiqc_files = ch_multiqc_files.mix(BRACKEN.out.txt.collect{ it -> it[1] })
            }
        }
    }
    //
    // SUBWORKFLOW: Pseudoalignment and quantification with Salmon
    //
    if (!params.skip_pseudo_alignment && params.pseudo_aligner) {

        if (params.pseudo_aligner == 'salmon') {
            ch_pseudo_index = ch_salmon_index
        } else {
            ch_pseudo_index = ch_kallisto_index
        }
        QUANTIFY_PSEUDO_ALIGNMENT (
            ch_samplesheet,
            ch_strand_inferred_filtered_fastq,
            ch_pseudo_index,
            ch_dummy_file,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            params.pseudo_aligner,
            false,
            params.salmon_quant_libtype ?: '',
            params.kallisto_quant_fraglen,
            params.kallisto_quant_fraglen_sd
        )
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_PSEUDO_ALIGNMENT.out.multiqc.collect{ it -> it[1] })
        ch_versions = ch_versions.mix(QUANTIFY_PSEUDO_ALIGNMENT.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_PSEUDO (
                QUANTIFY_PSEUDO_ALIGNMENT.out.counts_gene_length_scaled.map { it -> it[1] },
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_PSEUDO.out.pca_multiqc.collect())
            ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC_PSEUDO.out.dists_multiqc.collect())
            ch_versions = ch_versions.mix(DESEQ2_QC_PSEUDO.out.versions)
        }
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_rnaseq_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = Channel.empty()

    if (!params.skip_multiqc) {

        // Load MultiQC configuration files
        ch_multiqc_config        = Channel.fromPath("$projectDir/workflows/rnaseq/assets/multiqc/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo)   : Channel.empty()

        // Prepare the workflow summary
        ch_workflow_summary = Channel.value(
            paramsSummaryMultiqc(
                paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
            )
        ).collectFile(name: 'workflow_summary_mqc.yaml')

        // Prepare the methods section
        ch_methods_description = Channel.value(
            methodsDescriptionText(
                params.multiqc_methods_description
                    ? file(params.multiqc_methods_description)
                    : file("$projectDir/workflows/rnaseq/assets/multiqc/methods_description_template.yml", checkIfExists: true)
            )
        ).collectFile(name: 'methods_description_mqc.yaml')

        // Add summary, versions, and methods to the MultiQC input file list
        ch_multiqc_files = ch_multiqc_files
            .mix(ch_workflow_summary)
            .mix(ch_collated_versions)
            .mix(ch_methods_description)

        // Provide MultiQC with rename patterns to ensure it uses sample names
        // for single-techrep samples not processed by CAT_FASTQ, and trims out
        // _raw or _trimmed

        ch_name_replacements = ch_fastq
            .map{ meta, reads ->
                def name1 = file(reads[0][0]).simpleName + "\t" + meta.id + '_1'
                def fastqcnames = meta.id + "_raw\t" + meta.id + "\n" + meta.id + "_trimmed\t" + meta.id
                if (reads[0][1] ){
                    def name2 = file(reads[0][1]).simpleName + "\t" + meta.id + '_2'
                    def fastqcnames1 = meta.id + "_raw_1\t" + meta.id + "_1\n" + meta.id + "_trimmed_1\t" + meta.id + "_1"
                    def fastqcnames2 = meta.id + "_raw_2\t" + meta.id + "_2\n" + meta.id + "_trimmed_2\t" + meta.id + "_2"
                    return [ name1, name2, fastqcnames1, fastqcnames2 ]
                } else{
                    return [ name1, fastqcnames ]
                }
            }
            .flatten()
            .collectFile(name: 'name_replacement.txt', newLine: true)

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            ch_name_replacements,
            []
        )
        ch_multiqc_report = MULTIQC.out.report
    }

    emit:
    trim_status    = ch_trim_status    // channel: [id, boolean]
    map_status     = ch_map_status     // channel: [id, boolean]
    strand_status  = ch_strand_status  // channel: [id, boolean]
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
