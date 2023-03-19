//
// Uncompress and prepare reference genome files
//
include { SEPARATE_AUTOSOMAL_CHROMS         } from '../../modules/local/bcftools/separate_autosomal_chroms'
include { VCF2VCI                           } from '../../modules/local/g2gtools/vcf2vci/main'
include { PATCH as PATCH_REF                } from '../../modules/local/g2gtools/patch/main'
include { TRANSFORM                         } from '../../modules/local/g2gtools/transform/main'
include { GTF2DB                            } from '../../modules/local/g2gtools/gtf2db/main'
include { CONVERT                           } from '../../modules/local/g2gtools/convert/main'
include { EXTRACT                           } from '../../modules/local/g2gtools/extract/main'
include { MERGE_DIPLOID_HAPLOID_REFS        } from '../../modules/local/g2gtools/merge_diploid_haploid_refs'

include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GFF              } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../modules/nf-core/modules/gunzip/main'

include { UNTAR as UNTAR_BBSPLIT_INDEX      } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_RSEM_INDEX         } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../modules/nf-core/modules/untar/main'

include { CUSTOM_GETCHROMSIZES              } from '../../modules/nf-core/modules/custom/getchromsizes/main'
include { STAR_GENOMEGENERATE               } from '../../modules/nf-core/modules/star/genomegenerate/main'
include { HISAT2_EXTRACTSPLICESITES         } from '../../modules/nf-core/modules/hisat2/extractsplicesites/main'
include { HISAT2_BUILD                      } from '../../modules/nf-core/modules/hisat2/build/main'
include { SALMON_INDEX                      } from '../../modules/nf-core/modules/salmon/index/main'
include { RSEM_PREPAREREFERENCE as RSEM_PREPAREREFERENCE_GENOME } from '../../modules/nf-core/modules/rsem/preparereference/main'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../modules/nf-core/modules/rsem/preparereference/main'

include { GTF2BED                           } from '../../modules/local/gtf2bed'
include { EXTRACT_INTRONIC_REGIONS          } from '../../modules/local/bedtools_extract_intronic_regions'


workflow PERSONALIZE_REFERENCES {
    take:
        samples  // meta information for each sample
        ch_fasta // path: reference fasta
        ch_fai   // path: reference fasta index
        ch_gtf   // path: reference gtf

    main:
    ch_versions = Channel.empty()

    vcf = Channel.fromPath(params.vcf)
    vcf_index = Channel.fromPath(params.vcf_index)

    ref_fasta = ch_fasta.first().map{it[1]}
    ref_fai = ch_fai.first().map{it[1]}
    ref_gtf = ch_gtf.first().map{it[1]}

    SEPARATE_AUTOSOMAL_CHROMS(vcf, vcf_index, ref_fasta, ref_fai, ref_gtf)

    diploid_female_refs = SEPARATE_AUTOSOMAL_CHROMS.out.diploid_female_refs.collect()
    haploid_female_refs = SEPARATE_AUTOSOMAL_CHROMS.out.haploid_female_refs.collect()
    diploid_male_refs = SEPARATE_AUTOSOMAL_CHROMS.out.diploid_male_refs.collect()
    haploid_male_refs = SEPARATE_AUTOSOMAL_CHROMS.out.haploid_male_refs.collect()

    samples.multiMap { meta -> 
            haploid: [[id: meta.id + '-haploid',
                       sex: meta.sex,
                       sample: meta.id,
                       single_end: meta.single_end,
                       strandedness: meta.strandedness,
                       ploidy: 'haploid']]
            
            diploid: [[id: meta.id + '-diploid',
                       sex: meta.sex,
                       sample: meta.id,
                       single_end: meta.single_end,
                       strandedness: meta.strandedness,
                       ploidy: 'diploid']]
        }.set{samples}
    

    samples.diploid.branch{ meta ->
                        male: meta.sex[0]=='XY'
                        female: meta.sex[0]=='XX'
                    }.set{meta_diploid}
    
    samples.haploid.branch{ meta ->
                        male: meta.sex[0]=='XY' 
                        female: meta.sex[0]=='XX'
                    }.set{meta_haploid}

    diploid_males = meta_diploid.male.combine(diploid_male_refs)
    haploid_males = meta_haploid.male.combine(haploid_male_refs)
    diploid_females = meta_diploid.female.combine(diploid_female_refs)
    haploid_females = meta_haploid.female.combine(haploid_female_refs)

    samples_with_refs = diploid_males.mix(diploid_females, haploid_males, haploid_females)

    VCF2VCI(samples_with_refs)
    vci = VCF2VCI.out.vci
    

    PATCH_REF(vci.join(samples_with_refs.map{[it[0], //meta
                                                it[1], //fasta
                                                it[2], //fai
                                                it[3]  //gzi
                                                ]}))
    patched_ref = PATCH_REF.out.reference_with_snps

    TRANSFORM(vci.join(patched_ref))
    CONVERT(vci.join(samples_with_refs.map{[it[0], //meta
                                              it[6]  //gtf
                                              ]}), 'gtf')

    ch_gtf = CONVERT.out.converted_file
    GTF2DB(ch_gtf)

    transformed_fasta = TRANSFORM.out.fasta // meta, diploid_ref

    EXTRACT(transformed_fasta.join(GTF2DB.out.db))
    ch_transcriptome = EXTRACT.out.transcripts

    all_transcriptomes_refs_per_sample = ch_transcriptome
                                         .join(transformed_fasta)
                                         .join(ch_gtf)
                                         .join(vci)
                                         .map { meta, transcriptome, ref, ref_fai, fai_gz, gtf, vci, vci_index -> 
                                                [[ id: meta.sample,
                                                   sex: meta.sex, 
                                                   single_end: meta.single_end,
                                                   strandedness: meta.strandedness
                                                 ], transcriptome, ref, ref_fai, fai_gz, gtf, vci, vci_index] }
                                         .groupTuple(by: 0, size: 2, sort: { it[1] })
                                         .map { it.flatten() }

    MERGE_DIPLOID_HAPLOID_REFS( all_transcriptomes_refs_per_sample )
    
    ch_fasta = MERGE_DIPLOID_HAPLOID_REFS.out.merged_refseq
    ch_transcript_fasta = MERGE_DIPLOID_HAPLOID_REFS.out.merged_transcriptome
    ch_gtf = MERGE_DIPLOID_HAPLOID_REFS.out.merged_gtf
    ch_merged_vci = MERGE_DIPLOID_HAPLOID_REFS.out.merged_vci
    

    //
    // Create gene BED annotation file from GTFs if required
    //
    if (!params.skip_rseqc && params.rseqc_modules.size() > 0) {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    GUNZIP_FASTA(ch_fasta.map{meta, fasta, fai, gzi -> [meta, fasta]})
    ch_unzipped_fastas_no_fai = GUNZIP_FASTA.out.gunzip
    CUSTOM_GETCHROMSIZES ( ch_unzipped_fastas_no_fai )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)
    ch_fasta_fai_unzipped = ch_unzipped_fastas_no_fai.join(ch_fai)
    //
    // Uncompress BBSplit index or generate from scratch if required
    //
    ch_bbsplit_index = Channel.empty()
    // BBSplit not implemented - I don't understand what's going on with indexer

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_for_genome_generation = ch_fasta_fai_unzipped.join(ch_gtf)

    ch_star_index = Channel.empty()
    if (params.aligner == 'star_salmon') {
        ch_star_index = STAR_GENOMEGENERATE ( ch_for_genome_generation ).index
        ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    //
    // Uncompress RSEM index or generate from scratch if required
    //
    ch_rsem_index = Channel.empty()
    if (params.aligner == 'star_rsem'){
            ch_rsem_index = RSEM_PREPAREREFERENCE_GENOME ( ch_fasta.join(ch_gtf) ).index
            ch_versions   = ch_versions.mix(RSEM_PREPAREREFERENCE_GENOME.out.versions)
    }

    //
    // Uncompress HISAT2 index or generate from scratch if required
    //
    ch_splicesites  = Channel.empty()
    ch_hisat2_index = Channel.empty()
    if (params.aligner == 'hisat2') {
        ch_splicesites  = HISAT2_EXTRACTSPLICESITES ( ch_gtf ).txt
        ch_hisat2_index = HISAT2_BUILD ( ch_fasta.join(ch_gtf).join(ch_splicesites) ).index
        ch_versions     = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
    }
    

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index = Channel.empty()
    if (params.aligner == 'salmon') {
            SALMON_INDEX ( ch_fasta.join(ch_transcript_fasta) )
            ch_salmon_index = SALMON_INDEX.out.index
            ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
        }


    EXTRACT_INTRONIC_REGIONS(ch_gtf.join(ch_fai))
    
    //outtuple: [ val(meta), [ reads, salmon_index, transcript.fasta, genome.gtf ] ]
    // outtuple = ch_salmon_index.join(personalized_ref_fasta).map{meta, index, fasta-> [meta, index, fasta, []]}

    ch_versions = ch_versions.mix(VCF2VCI.out.versions)

    emit:
        ch_fasta            = ch_fasta
        ch_fasta_unzipped   = ch_unzipped_fastas_no_fai
        ch_fai              = ch_fai
        ch_gtf              = ch_gtf
        ch_gtf_unzipped     = ch_gtf
        vci                 = ch_merged_vci
        exon_bed            = EXTRACT_INTRONIC_REGIONS.out.exonic_regions
        intron_bed          = EXTRACT_INTRONIC_REGIONS.out.intronic_regions
        intergenic_bed      = EXTRACT_INTRONIC_REGIONS.out.intergenic_regions
        ch_gene_bed         = ch_gene_bed
        ch_transcript_fasta = ch_transcript_fasta
        ch_chrom_sizes      = ch_chrom_sizes
        ch_splicesites      = ch_splicesites
        ch_bbsplit_index    = ch_bbsplit_index
        ch_star_index       = ch_star_index
        ch_rsem_index       = ch_rsem_index
        ch_hisat2_index     = ch_hisat2_index
        ch_salmon_index     = ch_salmon_index
    // personal_transcriptome_tuple = outtuple

    versions                  = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
