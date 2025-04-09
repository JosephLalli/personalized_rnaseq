//
// Run bedClip and bedGraphToBigWig
//

include { UCSC_BEDCLIP          } from '../../../modules/nf-core/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG {
    take:
    bedgraph // channel: [ val(meta), [ bedgraph ] ]
    sizes    // channel: [ val(meta), [ chrom_sizes ] ]

    main:

    ch_versions = Channel.empty()

    //
    // Clip bedGraph file
    //
    UCSC_BEDCLIP ( bedgraph.join(sizes) )
    ch_versions = ch_versions.mix(UCSC_BEDCLIP.out.versions.first())

    //
    // Convert bedGraph to bigWig
    //
    UCSC_BEDGRAPHTOBIGWIG ( UCSC_BEDCLIP.out.bedgraph.join(sizes) )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions.first())

    emit:
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), [ bigwig ] ]
    bedgraph = UCSC_BEDCLIP.out.bedgraph        // channel: [ val(meta), [ bedgraph ] ]

    versions = ch_versions                      // channel: [ versions.yml ]
}
