#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnaseq
    Website: https://nf-co.re/rnaseq
    Slack  : https://nfcore.slack.com/channels/rnaseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// params.fasta            = getGenomeAttribute('fasta')
// params.additional_fasta = getGenomeAttribute('additional_fasta')
// params.transcript_fasta = getGenomeAttribute('transcript_fasta')
// params.gff              = getGenomeAttribute('gff')
// params.gtf              = getGenomeAttribute('gtf')
// params.gene_bed         = getGenomeAttribute('bed12')
// params.bbsplit_index    = getGenomeAttribute('bbsplit')
// params.sortmerna_index  = getGenomeAttribute('sortmerna')
// params.star_index       = getGenomeAttribute('star')
// params.rsem_index       = getGenomeAttribute('rsem')
// params.hisat2_index     = getGenomeAttribute('hisat2')
// params.salmon_index     = getGenomeAttribute('salmon')
// params.kallisto_index   = getGenomeAttribute('kallisto')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList                } from 'plugin/nf-schema'
include { checkSamplesAfterGrouping      } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { RNASEQ                   } from './workflows/rnaseq'
include { PREPARE_GENOME           } from './subworkflows/local/prepare_genome'
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS              } from './subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
include { PIPELINE_INITIALISATION  } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PIPELINE_COMPLETION      } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { checkMaxContigSize       } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow NFCORE_RNASEQ {
    take:
    params

    main:

    ch_versions = Channel.empty()
    def ref_meta = ['id': 'reference']
    ch_ref_meta = Channel.value(ref_meta)

    ch_samplesheet = Channel.value(file(params.input, checkIfExists: true))


    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            checkSamplesAfterGrouping(samplesheet)
        }
        .set { ch_fastq }


    // Todo: change order to filter reads and determine strandedness before prepare genome
    // The meta "strandedness" needs to be consistent across all reference channels
    //
    // SUBWORKFLOW: Prepare reference genome files
    //
    PREPARE_GENOME (
        params,
        ch_ref_meta,
        params.fasta,
        params.fai,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        [],
        [],
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
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME
            .out
            .fai
            .map { _meta, fai -> checkMaxContigSize(fai) }
    }

    //
    // WORKFLOW: Run nf-core/rnaseq workflow
    //
    RNASEQ (
        params,
        ch_samplesheet,
        ch_fastq,
        ch_versions,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.gene_bed,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.star_index,
        PREPARE_GENOME.out.rsem_index,
        PREPARE_GENOME.out.hisat2_index,
        PREPARE_GENOME.out.salmon_index,
        PREPARE_GENOME.out.kallisto_index,
        PREPARE_GENOME.out.bbsplit_index,
        PREPARE_GENOME.out.rrna_fastas,
        PREPARE_GENOME.out.sortmerna_index,
        PREPARE_GENOME.out.splicesites,
    )
    ch_versions = ch_versions.mix(RNASEQ.out.versions)

    emit:
    trim_status    = RNASEQ.out.trim_status    // channel: [id, boolean]
    map_status     = RNASEQ.out.map_status     // channel: [id, boolean]
    strand_status  = RNASEQ.out.strand_status  // channel: [id, boolean]
    multiqc_report = RNASEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions               // channel: [version1, version2, ...]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RNASEQ (params)

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_RNASEQ.out.multiqc_report,
        NFCORE_RNASEQ.out.trim_status,
        NFCORE_RNASEQ.out.map_status,
        NFCORE_RNASEQ.out.strand_status
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

// def getGenomeAttribute(attribute) {
//     if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
//         if (params.genomes[ params.genome ].containsKey(attribute)) {
//             return params.genomes[ params.genome ][ attribute ]
//         }
//     }
//     return null
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
