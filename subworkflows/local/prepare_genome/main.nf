//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../../modules/nf-core/gunzip'

include { UNTAR as UNTAR_BBSPLIT_INDEX      } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SORTMERNA_INDEX    } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_RSEM_INDEX         } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_KALLISTO_INDEX     } from '../../../modules/nf-core/untar'

include { CUSTOM_CATADDITIONALFASTA         } from '../../../modules/nf-core/custom/catadditionalfasta'
include { CUSTOM_GETCHROMSIZES              } from '../../../modules/nf-core/custom/getchromsizes'
include { SAMTOOLS_FAIDX                    } from '../../../modules/nf-core/samtools/faidx/main'
include { GFFREAD                           } from '../../../modules/nf-core/gffread'
include { BBMAP_BBSPLIT                     } from '../../../modules/nf-core/bbmap/bbsplit'
include { SORTMERNA as SORTMERNA_INDEX      } from '../../../modules/nf-core/sortmerna'
include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate'
include { HISAT2_EXTRACTSPLICESITES         } from '../../../modules/nf-core/hisat2/extractsplicesites'
include { HISAT2_BUILD                      } from '../../../modules/nf-core/hisat2/build'
include { SALMON_INDEX                      } from '../../../modules/nf-core/salmon/index'
include { KALLISTO_INDEX                    } from '../../../modules/nf-core/kallisto/index'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../../modules/nf-core/rsem/preparereference'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../../modules/nf-core/rsem/preparereference'

include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../../modules/local/preprocess_transcripts_fasta_gencode'
include { GTF2BED                              } from '../../../modules/local/gtf2bed'
include { GTF_FILTER                           } from '../../../modules/local/gtf_filter'
include { STAR_GENOMEGENERATE_IGENOMES         } from '../../../modules/local/star_genomegenerate_igenomes'

include { BCFTOOLS_INDEX                       } from '../../../modules/nf-core/bcftools/index'
include { PERSONALIZE_REFERENCES               } from '../personalize_references'

workflow PREPARE_GENOME {
    take:
    params
    sample_metadata          //       val: [ meta ]
    ref_fasta                //      file: /path/to/genome.fasta
    ref_fai                  //      file: /path/to/genome.fasta.fai
    ref_gtf                  //      file: /path/to/genome.gtf
    ref_gff                  //      file: /path/to/genome.gff
    additional_fasta         //      file: /path/to/additional.fasta
    transcript_fasta         //      file: /path/to/transcript.fasta
    vcf                      //      file: /path/to/phased_sample_variant_calls_to_create_personalized_references.vcf.gz
    vcf_index                //      file: /path/to/phased_sample_variant_calls_to_create_personalized_references.vcf.gz.tbi
    gene_bed                 //      file: /path/to/gene.bed
    par_bed                  //      file: /path/to/chrX_par.bed
    splicesites              //      file: /path/to/splicesites.txt
    bbsplit_fasta_list       //      file: /path/to/bbsplit_fasta_list.txt
    sortmerna_fasta_list     //      file: /path/to/sortmerna_fasta_list.txt
    star_index               // directory: /path/to/star/index/
    rsem_index               // directory: /path/to/rsem/index/
    salmon_index             // directory: /path/to/salmon/index/
    kallisto_index           // directory: /path/to/kallisto/index/
    hisat2_index             // directory: /path/to/hisat2/index/
    bbsplit_index            // directory: /path/to/bbsplit/index/
    sortmerna_index          // directory: /path/to/sortmerna/index/
    gencode                  //   boolean: whether the genome is from GENCODE
    featurecounts_group_type //    string: The attribute type used to group feature types in the GTF file when generating the biotype plot with featureCounts
    aligner                  //    string: Specifies the alignment algorithm to use - available options are 'star_salmon', 'star_rsem' and 'hisat2'
    pseudo_aligner           //    string: Specifies the pseudo aligner to use - available options are 'salmon'. Runs in addition to '--aligner'
    skip_gtf_filter          //   boolean: Skip filtering of GTF for valid scaffolds and/ or transcript IDs
    skip_bbsplit             //   boolean: Skip BBSplit for removal of non-reference genome reads
    skip_sortmerna           //   boolean: Skip sortmerna for removal of reads mapping to sequences in sortmerna_fasta_list
    skip_alignment           //   boolean: Skip all of the alignment-based processes within the pipeline
    skip_pseudo_alignment    //   boolean: Skip all of the pseudoalignment-based processes within the pipeline
    // personalize_reference    //   boolean: Personalize the reference genome using the provided VCF file

    main:
    ch_fasta    = Channel.empty()
    ch_fai      = Channel.empty()
    ch_gtf      = Channel.empty()
    ch_versions = Channel.empty()

    def ref_meta = ['id': 'reference']


    //
    // Uncompress genome fasta file if required
    //
    if (ref_fasta.endsWith('.gz')) {
        ch_ref_fasta    = GUNZIP_FASTA ( [ ref_meta , file(ref_fasta, checkIfExists: true) ] ).gunzip
        ch_ref_fai      = SAMTOOLS_FAIDX ( ch_ref_fasta, Channel.of([[], []]), false ).fai
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_ref_fasta = Channel.value([ ref_meta, file(ref_fasta, checkIfExists: true) ])
        if (ref_fai) {
            ch_ref_fai = Channel.value([ ref_meta, file(ref_fai, checkIfExists: true) ])
        } else {
            ch_ref_fai = SAMTOOLS_FAIDX ( ch_ref_fasta, Channel.of([[], []]), false).fai
        }
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (ref_gtf || ref_gff) {
        if (ref_gtf) {
            ch_ref_gtf = Channel.value([ ref_meta, file(ref_gtf, checkIfExists: true) ])

            if (ref_gtf.endsWith('.gz')) {
                ch_ref_gtf      = GUNZIP_GTF ( ch_ref_gtf ).gunzip
                ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
            }
        } else if (ref_gff) {
            ch_ref_gff = Channel.value([ ref_meta, file(ref_gff, checkIfExists: true) ])
            if (ref_gff.endsWith('.gz')) {
                ch_ref_gff      = GUNZIP_GFF ( ch_ref_gff ).gunzip
                ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
            }
            ch_ref_gtf      = GFFREAD ( ch_ref_gff, [] ).gtf
            ch_versions = ch_versions.mix(GFFREAD.out.versions)
        }

        // Determine whether to filter the GTF or not
        def filter_gtf =
            ((
                // Condition 1: Alignment is required and aligner is set
                !skip_alignment && aligner
            ) ||
            (
                // Condition 2: Pseudoalignment is required and pseudoaligner is set
                !skip_pseudo_alignment && pseudo_aligner
            ) ||
            (
                // Condition 3: Transcript FASTA file is not provided
                !transcript_fasta
            )) &&
            (
                // Condition 4: --skip_gtf_filter is not provided
                !skip_gtf_filter
            )
        if (filter_gtf) {
            GTF_FILTER ( ch_ref_fasta.join(ch_ref_gtf) )
            ch_gtf = GTF_FILTER.out.genome_gtf
            ch_versions = ch_versions.mix(GTF_FILTER.out.versions)
        }
    }

    //
    // Uncompress additional fasta file and concatenate with reference fasta and gtf files
    //
    def biotype = gencode ? "gene_type" : featurecounts_group_type
    if (additional_fasta) {
        if (additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( [ ref_meta, file(additional_fasta, checkIfExists: true) ] ).gunzip
            ch_versions  = ch_versions.mix(GUNZIP_ADDITIONAL_FASTA.out.versions)
        } else {
            ch_add_fasta = Channel.value(file(additional_fasta, checkIfExists: true)).map{ it -> [ref_meta, it] }
        }
        CUSTOM_CATADDITIONALFASTA (
            ch_ref_fasta.join(ch_ref_gtf).join(ch_add_fasta.map { it -> [ ref_meta, it ] } ),
            biotype
        )
        ch_ref_fasta    = CUSTOM_CATADDITIONALFASTA.out.fasta
        ch_ref_gtf      = CUSTOM_CATADDITIONALFASTA.out.gtf
        ch_versions = ch_versions.mix(CUSTOM_CATADDITIONALFASTA.out.versions)
    }

    if (vcf){
        ch_vcf = Channel.value(file(vcf, checkIfExists: true)).map{ it -> [ref_meta, it] }
        if ( ! vcf_index ) {
            BCFTOOLS_INDEX(ch_vcf)
            ch_vcf_index = BCFTOOLS_INDEX.out.csi
            ch_vcf_index = ch_vcf_index.mix(BCFTOOLS_INDEX.out.tbi)
        } else {
            ch_vcf_index = Channel.value(file(vcf_index, checkIfExists: true)).map{ it -> [ref_meta, it] }
        }
    } else {
        ch_vcf = Channel.empty()
        ch_vcf_index = Channel.empty()
    }

    //
    // If mapping to personalized transcriptome, generate personal transcriptomes
    // ch_fasta, ch_fai, ch_gtf contain [ sample_metadata, personalized fasta/fai/gtf ]
    // All downstream resources will be generated from those personalized files
    // If mapping to reference transcriptome, ch_fasta/fai/gtf will contain [ ref_meta, reference fasta/fai/gtf ]
    // And all downstream resources will be generated from those references
    //
    if ( params.use_personalized_references && vcf ) {
        println ("Personalizing references using VCF file")
        ch_par_bed = Channel.value(file(par_bed, checkIfExists: true)).map{ it -> [ref_meta, it] }
        PERSONALIZE_REFERENCES(
                sample_metadata.filter{ it -> it.id != 'reference'}, // meta information for each sample
                ch_ref_fasta,    // path: reference fasta
                ch_ref_fai,      // path: reference fasta index
                ch_ref_gtf,      // path: reference gtf
                ch_vcf,          // path: phased vcf of sample genetic variation
                ch_vcf_index,    // path: index of phased vcf
                ch_par_bed       // path: reference_bed_file
        )

        ch_fasta            = PERSONALIZE_REFERENCES.out.ch_fasta
        ch_fai              = PERSONALIZE_REFERENCES.out.ch_fai
        ch_gtf              = PERSONALIZE_REFERENCES.out.ch_gtf
        ch_transcript_fasta = PERSONALIZE_REFERENCES.out.ch_transcript_fasta
        ch_vci              = PERSONALIZE_REFERENCES.out.ch_vci
        ch_versions         = ch_versions.mix(PERSONALIZE_REFERENCES.out.versions)
    } else {
        ch_fasta = ch_ref_fasta
        ch_fai = ch_ref_fai
        ch_gtf = ch_ref_gtf
        ch_vci = Channel.empty()
    }

    // Uncompress gene BED annotation file or create from GTF if required
    //
    if (gene_bed) {
        if (gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ ref_meta, file(gene_bed, checkIfExists: true) ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = Channel.value(file(gene_bed, checkIfExists: true)).map{ it -> [ref_meta, it] }
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }


    //
    // Uncompress transcript fasta file / create if required
    //
    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ ref_meta, file(transcript_fasta, checkIfExists: true) ] ).gunzip
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = Channel.value(file(transcript_fasta, checkIfExists: true)).map{ it -> [ref_meta, it] }
        }
        if (gencode) {
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE ( ch_transcript_fasta )
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
            ch_versions         = ch_versions.mix(PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.versions)
        }
    } else if ( ! ( params.use_personalized_references && vcf ) ) {  // personalize references subworkflow generates this file
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta.join(ch_gtf) ).transcript_fasta
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.join(ch_fai))
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Get list of indices that need to be created
    //
    def prepare_tool_indices = []
    if (!skip_bbsplit) { prepare_tool_indices << 'bbsplit' }
    if (!skip_sortmerna) { prepare_tool_indices << 'sortmerna' }
    if (!skip_alignment) { prepare_tool_indices << aligner }
    if (!skip_pseudo_alignment && pseudo_aligner) { prepare_tool_indices << pseudo_aligner }

    //
    // Uncompress BBSplit index or generate from scratch if required
    // Todo: I do not think the BBMAP_BBSPLIT call makes sense, and I do not know enough about the program to know what to do here
    // Need to implement
    //

    ch_bbsplit_index = Channel.empty()
    if ('bbsplit' in prepare_tool_indices) {
        if (bbsplit_index) {
            if (bbsplit_index.endsWith('.tar.gz')) {
                ch_bbsplit_index = UNTAR_BBSPLIT_INDEX ( [ ref_meta, bbsplit_index ] ).untar
                ch_versions      = ch_versions.mix(UNTAR_BBSPLIT_INDEX.out.versions)
            } else {
                ch_bbsplit_index = Channel.value(file(bbsplit_index, checkIfExists: true)).map{ it -> [ref_meta, it] }
            }
        } else {
            ch_bbsplit_index = Channel.empty()
            bbsplit_fasta_list
    //         Channel
    //             .from(file(bbsplit_fasta_list))
    //             .splitCsv() // Read in 2 column csv file: short_name,path_to_fasta
    //             .flatMap { id, bbsplit_fasta -> [ [ 'id', id ], [ 'fasta', file(bbsplit_fasta, checkIfExists: true) ] ] } // Flatten entries to be able to groupTuple by a common key
    //             .groupTuple()
    //             .map { it -> it[1] } // Get rid of keys and keep grouped values
    //             .collect { [ it ] } // Collect entries as a list to pass as "tuple val(short_names), path(path_to_fasta)" to module
    //             .set { ch_bbsplit_fasta_list }

    //         ch_bbsplit_index = BBMAP_BBSPLIT ( [ ref_meta, [], [] ], ch_fasta, ch_bbsplit_fasta_list, true ).index
    //         ch_versions      = ch_versions.mix(BBMAP_BBSPLIT.out.versions)
        }
    }


    //
    // Uncompress sortmerna index or generate from scratch if required
    //
    ch_sortmerna_index = Channel.empty()
    ch_rrna_fastas = Channel.empty()

    if ('sortmerna' in prepare_tool_indices) {
        ribo_db = file(sortmerna_fasta_list, checkIfExists: true)

        // SortMeRNA needs the rRNAs even if we're providing the index
        ch_rrna_fastas = Channel.from(ribo_db.readLines())
            .map { row -> file(row, checkIfExists: true) }

        if (sortmerna_index) {
            if (sortmerna_index.endsWith('.tar.gz')) {
                ch_sortmerna_index = UNTAR_SORTMERNA_INDEX ( [ [ ref_meta ], sortmerna_index ] ).untar
                ch_versions = ch_versions.mix(UNTAR_SORTMERNA_INDEX.out.versions)
            } else {
                ch_sortmerna_index = Channel.value([[ ref_meta ], file(sortmerna_index, checkIfExists: true)])
            }
        } else {
            ch_sortmerna_index = Channel.empty()
            // SORTMERNA_INDEX (
            //     Channel.of([ [],[] ]),
            //     ch_rrna_fastas.collect().map { [ 'rrna_refs', it ] },
            //     Channel.of([ [],[] ])
            // )
            // ch_sortmerna_index = SORTMERNA_INDEX.out.index.first()
            // ch_versions = ch_versions.mix(SORTMERNA_INDEX.out.versions)
        }
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if ('star_salmon' in prepare_tool_indices) {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ ref_meta, star_index ] ).untar
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.value(file(star_index, checkIfExists: true)).map{ it -> [ref_meta, it] }
            }
        } else {
            // Check if an AWS iGenome has been provided to use the appropriate version of STAR
            def is_aws_igenome = false
            if (ref_fasta && ref_fasta) {
                if ((file(ref_fasta).getName() - '.gz' == 'genome.fa') && (file(ref_fasta).getName() - '.gz' == 'genes.gtf')) {
                    is_aws_igenome = true
                }
            }
            if (is_aws_igenome) {
                ch_star_index = STAR_GENOMEGENERATE_IGENOMES ( ch_fasta.join(ch_gtf).join(ch_vcf) ).index
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE_IGENOMES.out.versions)
            } else {
                ch_star_index = STAR_GENOMEGENERATE ( ch_fasta.join(ch_fai).join(ch_gtf).join(ch_vcf) ).index
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
            }
        }
    }

    //
    // Uncompress RSEM index or generate from scratch if required
    //
    ch_rsem_index = Channel.empty()
    if ('star_rsem' in prepare_tool_indices) {
        if (rsem_index) {
            if (rsem_index.endsWith('.tar.gz')) {
                ch_rsem_index = UNTAR_RSEM_INDEX ( [ ref_meta, rsem_index ] ).untar.map { it -> it[1] }.map{ it -> [ref_meta, it] }
                ch_versions   = ch_versions.mix(UNTAR_RSEM_INDEX.out.versions)
            } else {
                ch_rsem_index = file(rsem_index).map{ it -> [ref_meta, it] }
            }
        } else {
            ch_rsem_index = RSEM_PREPAREREFERENCE_GENOME ( ch_fasta.join(ch_gtf) ).index
            ch_versions   = ch_versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
        }
    }

    //
    // Uncompress HISAT2 index or generate from scratch if required
    //
    ch_splicesites  = Channel.empty()
    ch_hisat2_index = Channel.empty()
    if ('hisat2' in prepare_tool_indices) {
        if (!splicesites) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf ).txt
            ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        } else {
            ch_splicesites = file(splicesites, checkIfExists: true).map{ it -> [ref_meta, it] }
        }
        if (hisat2_index) {
            if (hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ( [ ref_meta, hisat2_index ] ).untar.map { it -> it[1] }.map{ it -> [ref_meta, it] }
                ch_versions     = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
            } else {
                ch_hisat2_index = file(hisat2_index, checkIfExists: true).map{ it -> [ref_meta, it] }
            }
        } else {
            ch_hisat2_index = HISAT2_BUILD ( ch_fasta.join(ch_gtf).join(ch_splicesites) ).index
            ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
    }

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index = Channel.empty()
    if (salmon_index) {
        if (salmon_index.endsWith('.tar.gz')) {
            ch_salmon_index = UNTAR_SALMON_INDEX ( [ ref_meta, salmon_index ] ).untar
            ch_versions     = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
        } else {
            ch_salmon_index = Channel.value(file(salmon_index, checkIfExists: true)).map{ it -> [ref_meta, it] }
        }
    } else {
        if ('salmon' in prepare_tool_indices) {
            ch_salmon_index = SALMON_INDEX ( ch_fasta.map{ it -> [ it[0], it[1] ]}.join(ch_transcript_fasta.map{ it -> [ it[0], it[1] ]}) ).index
            ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
        }
    }

    //
    // Uncompress Kallisto index or generate from scratch if required
    //
    ch_kallisto_index = Channel.empty()
    if (kallisto_index) {
        if (kallisto_index.endsWith('.tar.gz')) {
            ch_kallisto_index = UNTAR_KALLISTO_INDEX ( [ ref_meta, kallisto_index ] ).untar
            ch_versions     = ch_versions.mix(UNTAR_KALLISTO_INDEX.out.versions)
        } else {
            ch_kallisto_index = Channel.value([ ref_meta, file(kallisto_index, checkIfExists: true)])
        }
    } else {
        if ('kallisto' in prepare_tool_indices) {
            ch_kallisto_index = KALLISTO_INDEX ( ch_transcript_fasta ).index
            ch_versions     = ch_versions.mix(KALLISTO_INDEX.out.versions)
        }
    }

    // if not personalizing references, drop the reference tuple and instead
    // create each channel as [ [sample_metadata], refernece_index_object ]
    // If personalizing references, this work will have already been done
    // Net result is that no matter what, this pipeline emits
    // One value per sample, along with the index to be used for that sample

    if ( ! ( params.use_personalized_references && vcf ) ) {
        ch_fasta = sample_metadata.combine(ch_fasta.map { it -> it[1] })
        ch_gtf = sample_metadata.combine(ch_gtf.map { it -> it[1] })
        ch_fai = sample_metadata.combine(ch_fai.map { it -> it[1] })
        ch_gene_bed = sample_metadata.combine(ch_gene_bed.map { it -> it[1] })
        ch_transcript_fasta = sample_metadata.combine(ch_transcript_fasta.map { it -> it[1] })
        ch_vci = sample_metadata.combine(ch_vci.map { it -> it[1] })
        ch_chrom_sizes = sample_metadata.combine(ch_chrom_sizes.map { it -> it[1] })
        ch_splicesites = sample_metadata.combine(ch_splicesites.map { it -> it[1] })
        ch_bbsplit_index = sample_metadata.combine(ch_bbsplit_index.map { it -> it[1] })
        ch_rrna_fastas = sample_metadata.combine(ch_rrna_fastas.map { it -> it[1] })
        ch_sortmerna_index = sample_metadata.combine(ch_sortmerna_index.map { it -> it[1] })
        ch_star_index = sample_metadata.combine(ch_star_index.map { it -> it[1] })
        ch_rsem_index = sample_metadata.combine(ch_rsem_index.map { it -> it[1] })
        ch_hisat2_index = sample_metadata.combine(ch_hisat2_index.map { it -> it[1] })
        ch_salmon_index = sample_metadata.combine(ch_salmon_index.map { it -> it[1] })
        ch_kallisto_index = sample_metadata.combine(ch_kallisto_index.map { it -> it[1] })
    }

    emit:
    fasta            = ch_fasta                  // channel: [ meta, path(genome.fasta) ]
    gtf              = ch_gtf                    // channel: [ meta, path(genome.gtf) ]
    fai              = ch_fai                    // channel: [ meta, path(genome.fai) ]
    gene_bed         = ch_gene_bed               // channel: [ meta, path(gene.bed) ]
    transcript_fasta = ch_transcript_fasta       // channel: [ meta, path(transcript.fasta) ]
    vci              = ch_vci                    // channel: [ meta, path(gene.vci) ]
    chrom_sizes      = ch_chrom_sizes            // channel: [ meta, path(genome.sizes) ]
    splicesites      = ch_splicesites            // channel: [ meta, path(genome.splicesites.txt) ]
    bbsplit_index    = ch_bbsplit_index          // channel: [ meta, path(bbsplit/index/) ]
    rrna_fastas      = ch_rrna_fastas            // channel: [ meta, path(sortmerna_fasta_list) ]
    sortmerna_index  = ch_sortmerna_index        // channel: [ meta, path(sortmerna/index/) ]
    star_index       = ch_star_index             // channel: [ meta, path(star/index/) ]
    rsem_index       = ch_rsem_index             // channel: [ meta, path(rsem/index/) ]
    hisat2_index     = ch_hisat2_index           // channel: [ meta, path(hisat2/index/) ]
    salmon_index     = ch_salmon_index           // channel: [ meta, path(salmon/index/) ]
    kallisto_index   = ch_kallisto_index         // channel: [ meta, path(kallisto/index/) ]
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
