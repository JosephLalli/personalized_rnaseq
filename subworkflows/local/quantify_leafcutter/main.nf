//
// Intron/splice site usage quantification with leafcutter
//

include { REGTOOLS_BAMTOJUNC               } from '../../../modules/local/regtools/junctions/extract'
include { LEAFCUTTER_CLUSTERINTRONS        } from '../../../modules/local/leafcutter/leafcutter_cluster_regtools'
include { CONVERT as G2GTOOLS_CONVERT_JUNC } from '../../../modules/local/g2gtools/convert'

workflow QUANTIFY_LEAFCUTTER {
    take:
    bam_sorted_indexed  // channel: [ val(meta), [ bam ], [ bai ] ]
    vci                 // channel: [ val(meta), vci.gz, vci.gz.tbi ]

    main:
    ch_versions = Channel.empty()


    // extract junctions from bams
    REGTOOLS_BAMTOJUNC( bam_sorted_indexed )

    // if aligning to personalized genomes, remap back to reference genomes
    if (params.use_personalized_references) {
        G2GTOOLS_CONVERT_JUNC( vci.join( REGTOOLS_BAMTOJUNC.out.junctions ), 'junc' )
        junc_files = G2GTOOLS_CONVERT_JUNC.out.converted_file

        ch_versions = ch_versions.mix(G2GTOOLS_CONVERT_JUNC.out.versions)
    } else {
        junc_files = REGTOOLS_BAMTOJUNC.out.junctions
    }

    //
    // once all are done, cluster introns
    //
    junc_files_ch = junc_files.collect{ it -> it[1] }
    junc_files_list_ch = junc_files.map({ it -> it[1].toString()}).collectFile(name:'input_files.txt', newLine:true)

    LEAFCUTTER_CLUSTERINTRONS( junc_files_ch, junc_files_list_ch )

    ch_versions = ch_versions.mix(REGTOOLS_BAMTOJUNC.out.versions.first())
    ch_versions = ch_versions.mix(LEAFCUTTER_CLUSTERINTRONS.out.versions)

    emit:
    leafcutter_intron_counts = LEAFCUTTER_CLUSTERINTRONS.out.leafcutter_perind_counts       // channel: [ path(countsfile) ]
    versions                 = ch_versions
}
