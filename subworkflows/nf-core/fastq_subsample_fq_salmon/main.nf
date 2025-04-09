//
// Sub-sample FastQ files and pseudo-align with Salmon
//      can be used to infer strandedness of library
//

include { SALMON_INDEX } from '../../../modules/nf-core/salmon/index/main'
include { FQ_SUBSAMPLE } from '../../../modules/nf-core/fq/subsample/main'
include { SALMON_QUANT } from '../../../modules/nf-core/salmon/quant/main'

workflow FASTQ_SUBSAMPLE_FQ_SALMON {
    take:
    ch_reads            // channel: [ val(meta), [ reads ] ]
    ch_genome_fasta     // channel: /path/to/genome.fasta
    ch_transcript_fasta // channel: /path/to/transcript.fasta
    ch_gtf              // channel: /path/to/genome.gtf
    ch_index            // channel: /path/to/salmon/index/
    make_index          // boolean: Whether to create salmon index before running salmon quant

    main:

    ch_versions = Channel.empty()

    //
    // Create Salmon index if required
    //

    if (make_index) {
        ch_index = SALMON_INDEX ( ch_genome_fasta.join( ch_transcript_fasta )).index
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
    }
    
    //
    // Sub-sample FastQ files with fq
    //
    FQ_SUBSAMPLE ( ch_reads )
    ch_versions = ch_versions.mix(FQ_SUBSAMPLE.out.versions.first())

    //
    // Pseudo-alignment with Salmon
    //
    def lib_type = Channel.value('A')
    def alignment_mode = false

    quant_input = ch_reads.map { it -> it[0] }.combine(ch_index.join(ch_gtf).join(ch_transcript_fasta).map{ it -> it.drop(1) })
    // FQ_SUBSAMPLE.out.fastq.view{it -> "FQ_SUBSAMPLE.out.fastq: ${it}"}
    // ch_index.view{it -> "ch_index: ${it}"}
    // ch_index.join(ch_gtf).view{it -> "ch_index.join(ch_gtf): ${it}"}
    // ch_index.join(ch_gtf).join(ch_transcript_fasta).view{it -> "ch_index.join(ch_gtf).join(ch_transcript_fasta): ${it}"}
    // quant_input
    SALMON_QUANT ( FQ_SUBSAMPLE.out.fastq.join(quant_input)
                                         .combine(lib_type), alignment_mode)
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    if ( workflow.stubRun ) {
        def dummy_report = "$projectDir/subworkflows/nf-core/fastq_subsample_fq_salmon/assets/dummy_salmon_counts.json"
        ch_meta = FQ_SUBSAMPLE.out.fastq.map{ it -> it[0] }
        ch_lib_format_counts = ch_meta.combine(Channel.value(file(dummy_report)))
    } else {
        ch_lib_format_counts = SALMON_QUANT.out.lib_format_counts
    }

    emit:
    index             = ch_index                           // channel: [ index ]

    reads             = FQ_SUBSAMPLE.out.fastq             // channel: [ val(meta), fastq ]

    results           = SALMON_QUANT.out.results           // channel: [ val(meta), results_dir ]
    json_info         = SALMON_QUANT.out.json_info         // channel: [ val(meta), json_info
    lib_format_counts = ch_lib_format_counts               // channel: [ val(meta), json_info

    versions          = ch_versions                        // channel: [ versions.yml ]
}
