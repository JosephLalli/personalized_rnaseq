prepare star genome
    - includes extracting snps from vcf and preparing all personalized diploid ref genomes
align STAR
index, sort, mark duplicates w/ biobambam2
annotate stringtie all samples
merge stringtie annotations
align STAR again
    - wasp, bootstrap, custom diploid ref genomes
index, sort, mark duplicates w/ biobambam2

count reads w/ salmon
ASEReadCounter to extract allele-specific counts

seperate leafcutter run on product of biobambam duplicates
collect splice events at end

Alignment QC subworkflow:
    gather QC reports from Star and Salmon
    run Bamdiff on ref gff and stringtie gff for transcriptome comparison
    produce bigwig files for figures
    run featureCount to determine relative alignment rate to different classes of genes
    run rSeqQC:
        junction_annotation (splicing statistics)
        read_distribution (reads/kb in different categories of sequence)



Use built-in feature of STAR to align reads to diploid reference genomes with individual’s mutations in place (should increase alignment accuracy, help resolve edge cases (eg pseudogenes), allow for allele-specific counts)
Run w/ wasp (accounts for differential alignment between alleles)
Run w/ bootstrap (data can be used with R/fishpond to inform uncertainties around DE, potentially can be used in xQTL analyses too)
Mark duplicate reads w/ biobambam2 (same as picard, faster)
Stringtie to identify previously unknown transcripts in dataset
Salmon to locally remap reads to updated transcriptome
Leafcutter to independently produce splicing event data
Tximport to extract transcript and gene counts
ASEReadCounter to extract allele-specific counts
Alignment QC
QC reports from Star and Salmon
Bamdiff for transcriptome comparison
Bigwig files for figures
FeatureCount to determine relative alignment rate to different classes of genes
rSeqQC
junction_annotation (splicing statistics)
read_distribution (reads/kb in different categories of sequence)

Run geneQTL, transcriptQTL, splicingQTL, ASE-aware transcriptQTL
Will be formatting these data so that they are all the same format; will allow for these analyses to be done w/ just different phenotypes into the same pipeline.
Use mash and/or mvSusie and/or coloc to integrate these analyses to identify effector SNPs. (A snp that alters splicing and gene expression is more likely to be an effector SNP than one that affects only one or the other.)
