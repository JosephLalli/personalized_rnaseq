//
// Pseudoalignment and quantification with Salmon or Kallisto
//

include { SALMON_QUANT     } from '../../../modules/nf-core/salmon/quant'
include { KALLISTO_QUANT   } from '../../../modules/nf-core/kallisto/quant'
include { CUSTOM_TX2GENE   } from '../../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT } from '../../../modules/nf-core/tximeta/tximport'

include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE               } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE_LENGTH_SCALED } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE_SCALED        } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_TRANSCRIPT         } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'

workflow QUANTIFY_PSEUDO_ALIGNMENT {
    take:
    samplesheet               // channel: [ val(meta), /path/to/samplsheet ]
    ch_reads                     // channel: [ val(meta), [ reads ] ]
    index                     // channel: [ val(meta), /path/to/index/ ]
    transcript_fasta          // channel: [ val(meta), /path/to/transcript.fasta ]
    gtf                       // channel: [ val(meta), /path/to/genome.gtf ]
    gtf_id_attribute          //     val: GTF gene ID attribute
    gtf_extra_attribute       //     val: GTF alternative gene attribute (e.g. gene_name)
    pseudo_aligner            //     val: kallisto or salmon
    alignment_mode            //    bool: Run Salmon in alignment mode
    lib_type                  //     val: String to override Salmon library type
    kallisto_quant_fraglen    //     val: Estimated fragment length required by Kallisto in single-end mode
    kallisto_quant_fraglen_sd //     val: Estimated standard error for fragment length required by Kallisto in single-end mode

    main:
    ch_versions = Channel.empty()    

    ch_reads = ch_reads.map{meta, reads -> [meta + ['sample':meta.id], reads] }

    //
    // Quantify and merge counts across samples
    //
    // NOTE: MultiQC needs Salmon outputs, but Kallisto logs
    if (pseudo_aligner == 'salmon') {
        SALMON_QUANT (
            ch_reads.join(index).join(gtf).join(transcript_fasta.map{meta, fasta -> [meta, fasta, []]}),
            alignment_mode
        )
        ch_pseudo_results = SALMON_QUANT.out.results
        ch_pseudo_multiqc = ch_pseudo_results
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
    } else {
        KALLISTO_QUANT (
            ch_reads.join(index).join(gtf).join(transcript_fasta).map{meta, gtf_file, transcript_fasta_file -> [meta, gtf_file, transcript_fasta_file, []]},
            kallisto_quant_fraglen,
            kallisto_quant_fraglen_sd
        )
        ch_pseudo_results = KALLISTO_QUANT.out.results
        ch_pseudo_multiqc = KALLISTO_QUANT.out.log
        ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions.first())
    }

    // CUSTOM_TX2GENE (
    //     gtf.join(ch_pseudo_results),
    //     pseudo_aligner,
    //     gtf_id_attribute,
    //     gtf_extra_attribute
    // )
    // ch_versions = ch_versions.mix(CUSTOM_TX2GENE.out.versions)

    // reference_tx2gene = CUSTOM_TX2GENE.out.tx2gene.map{ it -> it[1] }//.view{ it -> "CUSTOM_TX2GENE.out.tx2gene: $it" }
    // // really not sure if tximeta_tximport should be working with all samples at once, or should be working with personalized tx2gene values
    
    // TXIMETA_TXIMPORT (
    //     ch_pseudo_results.view{it -> "TXIMETA_TXIMPORT: $it"}.collect{ it -> it[1] }.view{it -> "all_samples:\n$it"}map { it -> [ ['id': 'all_samples'], it ]
    //     }.concat(reference_tx2gene),
    //     pseudo_aligner
    // )
    // ch_versions = ch_versions.mix(TXIMETA_TXIMPORT.out.versions)

    // the following code assumes that SE Gene requires the reference tx2gene and samplesheet. If instead these items need to be personalized, then the code won't work.
    // TXIMETA_TXIMPORT.out.counts_gene.view{it -> "TXIMPORT.counts_gene:${it}"}.concat(TXIMETA_TXIMPORT.out.tpm_gene).view{it -> "TXIMPORT.tpm_gene:${it}"}.groupTuple().view{it -> "grouped:${it}"}.concat(reference_tx2gene.view{it -> "TXIMPORT.concatted w/reference_tx2gene:${it}"}).concat(samplesheet).view{it -> "SE_GENE input:\n${it}"}

    // SE_GENE (
    //     TXIMETA_TXIMPORT.out.counts_gene.view{it -> "TXIMPORT.counts_gene:${it}"}.concat(TXIMETA_TXIMPORT.out.tpm_gene).view{it -> "TXIMPORT.tpm_gene:${it}"}.groupTuple().view{it -> "grouped:${it}"}
    //                                     .concat(reference_tx2gene.view{it -> "TXIMPORT.concatted w/reference_tx2gene:${it}"})
    //                                     .concat(samplesheet),
    // )
    // ch_versions = ch_versions.mix(SE_GENE.out.versions)

    // SE_GENE_LENGTH_SCALED (
    //     TXIMETA_TXIMPORT.out.counts_gene_length_scaled.concat(TXIMETA_TXIMPORT.out.tpm_gene).groupTuple()
    //                                                   .concat(reference_tx2gene)
    //                                                   .concat(samplesheet),
    // )

    // SE_GENE_SCALED (
    //     TXIMETA_TXIMPORT.out.counts_gene_scaled.concat(TXIMETA_TXIMPORT.out.tpm_gene).groupTuple()
    //                                            .concat(reference_tx2gene)
    //                                            .concat(samplesheet),
    // )

    // SE_TRANSCRIPT (
    //     TXIMETA_TXIMPORT.out.counts_transcript.concat(TXIMETA_TXIMPORT.out.tpm_transcript).groupTuple()
    //                                           .concat(reference_tx2gene)
    //                                           .concat(samplesheet),
    // )

    emit:
    results                       = ch_pseudo_results                              // channel: [ val(meta), results_dir ]
    multiqc                       = ch_pseudo_multiqc                              // channel: [ val(meta), files_for_multiqc ]

    // tpm_gene                      = TXIMETA_TXIMPORT.out.tpm_gene                  //    path: *gene_tpm.tsv
    // counts_gene                   = TXIMETA_TXIMPORT.out.counts_gene               //    path: *gene_counts.tsv
    // lengths_gene                  = TXIMETA_TXIMPORT.out.lengths_gene              //    path: *gene_lengths.tsv
    // counts_gene_length_scaled     = TXIMETA_TXIMPORT.out.counts_gene_length_scaled //    path: *gene_counts_length_scaled.tsv
    // counts_gene_scaled            = TXIMETA_TXIMPORT.out.counts_gene_scaled        //    path: *gene_counts_scaled.tsv
    // tpm_transcript                = TXIMETA_TXIMPORT.out.tpm_transcript            //    path: *gene_tpm.tsv
    // counts_transcript             = TXIMETA_TXIMPORT.out.counts_transcript         //    path: *transcript_counts.tsv
    // lengths_transcript            = TXIMETA_TXIMPORT.out.lengths_transcript        //    path: *transcript_lengths.tsv

    // merged_gene_rds               = Channel.value('SE_GENE.out.rds')                                //    path: *.rds
    // merged_gene_rds_length_scaled = Channel.value('SE_GENE_LENGTH_SCALED.out.rds')                  //    path: *.rds
    // merged_gene_rds_scaled        = Channel.value('SE_GENE_SCALED.out.rds')                         //    path: *.rds
    // merged_transcript_rds         = Channel.value('SE_TRANSCRIPT.out.rds')                          //    path: *.rds

    versions                      = ch_versions                                    // channel: [ versions.yml ]
}
