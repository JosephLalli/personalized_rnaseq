//
// Uncompress and prepare reference genome files
//
include { SEPARATE_AUTOSOMAL_CHROMS         } from '../../../modules/local/bcftools/separate_autosomal_chroms'
include { VCF2VCI                           } from '../../../modules/local/g2gtools/vcf2vci/main'
include { PATCH as PATCH_REF                } from '../../../modules/local/g2gtools/patch/main'
include { TRANSFORM                         } from '../../../modules/local/g2gtools/transform/main'
include { GTF2DB                            } from '../../../modules/local/g2gtools/gtf2db/main'
include { CONVERT                           } from '../../../modules/local/g2gtools/convert/main'
include { EXTRACT                           } from '../../../modules/local/g2gtools/extract/main'
include { MERGE_DIPLOID_HAPLOID_REFS        } from '../../../modules/local/g2gtools/merge_diploid_haploid_refs'

include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF              } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../../modules/nf-core/gunzip/main'

include { UNTAR as UNTAR_BBSPLIT_INDEX      } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_RSEM_INDEX         } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../../modules/nf-core/untar/main'

include { CUSTOM_GETCHROMSIZES              } from '../../../modules/nf-core/custom/getchromsizes/main'
include { STAR_GENOMEGENERATE               } from '../../../modules/nf-core/star/genomegenerate/main'
include { HISAT2_EXTRACTSPLICESITES         } from '../../../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD                      } from '../../../modules/nf-core/hisat2/build/main'
include { SALMON_INDEX                      } from '../../../modules/nf-core/salmon/index/main'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../../modules/nf-core/rsem/preparereference/main'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../../modules/nf-core/rsem/preparereference/main'

include { GTF2BED                           } from '../../../modules/local/gtf2bed'
include { EXTRACT_INTRONIC_REGIONS          } from '../../../modules/local/bedtools_extract_intronic_regions'

workflow PERSONALIZE_REFERENCES {
    take:
        samples      // meta information for each sample
        ch_ref_fasta // path: reference fasta
        ch_ref_fai   // path: reference fasta index
        ch_ref_gtf   // path: reference gtf
        ch_vcf
        ch_vcf_index
        ch_par_bed   // path: reference chrX par bed file

    main:
    ch_versions = Channel.empty()

    //
    // Separate haploid regions for g2gtools
    // Unfortuantely, g2gtools cannot convert genomes where a portion of the genome is haploid.
    // This causes problems when handling chrX, chrY, and chrM.
    // While chrM support is only partially implemented, chrX/Y are vital.
    // So we split these contigs out and separately run g2gtools on the diploid and hapl
    //

    def reference_contigs = ch_ref_fai.map{it->it[1]}.first().splitCsv(sep: '\t').collect { it0 -> it0[0] }
    def par_regions = ch_par_bed.map{it -> it[1]}.first().splitCsv(sep: '\t').collect { line -> "${line[0]}:${line[1]}-${line[2]}" }

    SEPARATE_AUTOSOMAL_CHROMS(ch_vcf, ch_vcf_index, ch_ref_fasta, ch_ref_fai, ch_ref_gtf, ch_par_bed, reference_contigs, par_regions)
    diploid_female_refs = SEPARATE_AUTOSOMAL_CHROMS.out.diploid_female_refs.collect()
    haploid_female_refs = SEPARATE_AUTOSOMAL_CHROMS.out.haploid_female_refs.collect()
    diploid_male_refs = SEPARATE_AUTOSOMAL_CHROMS.out.diploid_male_refs.collect()
    haploid_male_refs = SEPARATE_AUTOSOMAL_CHROMS.out.haploid_male_refs.collect()

    ch_versions = ch_versions.mix(SEPARATE_AUTOSOMAL_CHROMS.out.versions)

    //todo:
    // Sample 676's test fastq needs to be regenerated - different number of R1 and R2 reads
    // Does the below logic need to happen if we just artifically make haploid sequences diploid in the vcf?
        // That would allow for a more graceful way to handle situations where chrM or chrY is not used
    // Probably still want to drop sequences that are not valid for organism from reference before mapping
    // Maybe after the fact, if there are any haploid regions in the dataset, drop the second haplotype from all references (fasta, gtf, transcriptome) before making indexes

    samples.multiMap { meta ->
            haploid: [ meta + [id: meta.id + '-haploid',
                               sample: meta.id,
                               ploidy: 'haploid']]

            diploid: [ meta + [id: meta.id + '-diploid',
                               sample: meta.id,
                               ploidy: 'diploid']]
        }.set{samples}

    samples.diploid.branch{ meta ->
                        male: meta.sex[0]=='XY' | meta.sex[0]=='Male' | meta.sex[0]=='male'
                        female: meta.sex[0]=='XX' | meta.sex[0]=='Female' | meta.sex[0] == 'female'
                    }.set{meta_diploid}

    samples.haploid.branch{ meta ->
                        male: meta.sex[0]=='XY' | meta.sex[0]=='Male' | meta.sex[0]=='male'
                        female: meta.sex[0]=='XX' | meta.sex[0]=='Female' | meta.sex[0] == 'female'
                    }.set{meta_haploid}

    diploid_males = meta_diploid.male.combine(diploid_male_refs)
    haploid_males = meta_haploid.male.combine(haploid_male_refs)
    diploid_females = meta_diploid.female.combine(diploid_female_refs)
    haploid_females = meta_haploid.female.combine(haploid_female_refs)

    samples_with_refs = diploid_males.mix(diploid_females, haploid_males, haploid_females)

    VCF2VCI(samples_with_refs)
    ch_vci = VCF2VCI.out.vci


    PATCH_REF(ch_vci.join(samples_with_refs.map{[it[0], //meta
                                                it[1], //fasta
                                                it[2], //fai
                                                it[3]  //gzi
                                                ]}))
    patched_ref = PATCH_REF.out.reference_with_snps

    TRANSFORM(ch_vci.join(patched_ref))
    CONVERT(ch_vci.join(samples_with_refs.map{[it[0], //meta
                                              it[6]  //gtf
                                              ]}), 'gtf')

    ch_gtf = CONVERT.out.converted_file
    GTF2DB(ch_gtf)

    transformed_fasta = TRANSFORM.out.fasta // meta, diploid_ref

    EXTRACT(transformed_fasta.join(GTF2DB.out.db))
    ch_transcriptome = EXTRACT.out.transcripts

    ch_versions = ch_versions.mix( VCF2VCI.out.versions,
                                   PATCH_REF.out.versions,
                                   TRANSFORM.out.versions,
                                   CONVERT.out.versions,
                                   GTF2DB.out.versions,
                                   EXTRACT.out.versions)

    all_transcriptomes_refs_per_sample = ch_transcriptome
                                         .join(transformed_fasta)
                                         .join(ch_gtf)
                                         .join(ch_vci)
                                         .map { meta, transcriptome, ref, ref_fai, fai_gz, gtf, vci, vci_index ->
                                                def keys = meta.keySet().findAll { it != 'ploidy' }
                                                [ meta.subMap(keys) + [ id: meta.sample ],
                                                  transcriptome, ref, ref_fai, fai_gz, gtf, vci, vci_index] }
                                         .groupTuple(by: 0, size: 2, remainder:true, sort: {a, _b ->
                                                                        if (a.toString().contains('haploid')) {
                                                                            return -1
                                                                        } else {
                                                                            return 1
                                                                        }})

    MERGE_DIPLOID_HAPLOID_REFS( all_transcriptomes_refs_per_sample)

    ch_versions = ch_versions.mix(MERGE_DIPLOID_HAPLOID_REFS.out.versions)

    ch_fasta_fai_gzi    = MERGE_DIPLOID_HAPLOID_REFS.out.merged_refseq
    ch_fasta            = ch_fasta_fai_gzi.map { meta, fasta, _fai -> [meta, fasta] }
    ch_fai              = ch_fasta_fai_gzi.map { meta, _fasta, fai -> [meta, fai] }
    ch_transcript_fasta = MERGE_DIPLOID_HAPLOID_REFS.out.merged_transcriptome
    ch_gtf              = MERGE_DIPLOID_HAPLOID_REFS.out.merged_gtf
    ch_merged_vci       = MERGE_DIPLOID_HAPLOID_REFS.out.merged_vci


    // //
    // // Create gene BED annotation file from GTFs if required
    // //
    // ch_gene_bed = Channel.empty()
    // if (! params.skip_rseqc && params.rseqc_modules.size() > 0) {
    //     ch_gene_bed = GTF2BED ( ch_gtf ).bed
    //     ch_versions = ch_versions.mix(GTF2BED.out.versions)
    // }

    // //
    // // Create chromosome sizes file
    // //
    // GUNZIP_FASTA(ch_fasta.map{meta, fasta, _fai, _gzi -> [meta, fasta]})
    // ch_unzipped_fastas_no_fai = GUNZIP_FASTA.out.gunzip
    // CUSTOM_GETCHROMSIZES ( ch_unzipped_fastas_no_fai )
    // ch_fai         = CUSTOM_GETCHROMSIZES.out.fai
    // ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes
    // ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)
    // ch_fasta_fai_unzipped = ch_unzipped_fastas_no_fai.join(ch_fai)

    // //
    // // Uncompress BBSplit index or generate from scratch if required
    // //
    // ch_bbsplit_index = Channel.empty()
    // // BBSplit not implemented - I don't understand what's going on with indexer

    // //
    // // Uncompress STAR index or generate from scratch if required
    // //
    // ch_for_genome_generation = ch_fasta_fai_unzipped.join(ch_gtf)

    // ch_star_index = Channel.empty()
    // if (!params.skip_alignment && params.aligner == 'star_salmon') {
    //     ch_star_index = STAR_GENOMEGENERATE ( ch_for_genome_generation ).index
    //     ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    // }

    // //
    // // Uncompress RSEM index or generate from scratch if required
    // //
    // ch_rsem_index = Channel.empty()
    // if (!params.skip_alignment && params.aligner == 'star_rsem'){
    //         ch_rsem_index = RSEM_PREPAREREFERENCE_GENOME ( ch_fasta.join(ch_gtf) ).index
    //         ch_versions   = ch_versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
    // }

    // //
    // // Uncompress HISAT2 index or generate from scratch if required
    // //
    // ch_splicesites  = Channel.empty()
    // ch_hisat2_index = Channel.empty()
    // if (!params.skip_alignment && params.aligner == 'hisat2') {
    //     ch_splicesites  = HISAT2_EXTRACTSPLICESITES ( ch_gtf ).txt
    //     ch_hisat2_index = HISAT2_BUILD ( ch_fasta.join(ch_gtf).join(ch_splicesites) ).index
    //     ch_versions     = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
    //     ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
    // }


    //
    // Uncompress Salmon index or generate from scratch if required
    //
    // ch_salmon_index = Channel.empty()
    // if (params.pseudo_aligner == 'salmon') {
    //     SALMON_INDEX ( ch_fasta.map{[it[0], it[1]]}.join(ch_transcript_fasta.map{[it[0], it[1]]}) )
    //     ch_salmon_index = SALMON_INDEX.out.index
    //     ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
    // }


    // EXTRACT_INTRONIC_REGIONS(ch_gtf.join(ch_fai))

    // //outtuple: [ val(meta), [ reads, salmon_index, transcript.fasta, genome.gtf ] ]
    // // outtuple = ch_salmon_index.join(personalized_ref_fasta).map{meta, index, fasta-> [meta, index, fasta, []]}

    // ch_versions = ch_versions.mix(VCF2VCI.out.versions)

    emit:
        ch_fasta            = ch_fasta
        ch_fai              = ch_fai
        ch_gtf              = ch_gtf
        ch_gtf_unzipped     = ch_gtf
        ch_vci              = ch_merged_vci
        ch_transcript_fasta = ch_transcript_fasta
    // personal_transcriptome_tuple = outtuple

    versions                  = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
