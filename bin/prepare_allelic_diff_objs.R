library(tidyr)
library(interleave)
library(data.table)
#Create coldata csv of experiment
coldata = read.csv('/mnt/sas0/AD/lalli/nf_stage/rnaseq_JLL/trimmed_brainvar_fastqs_corrected_samples_with_sex.csv')
meta = read.csv('/mnt/sas0/AD/lalli/nf_stage/rnaseq_JLL/brainvar_sample_metadata.csv')
coldata = merge(coldata, meta[,c('Braincode','AgeDays')], by.x='sample',by.y ='Braincode', T, F)
root='/mnt/sas0/AD/lalli/reference_comparison_results/T2Tv2_leung_ncbi110_personalized_oldcalls/salmon'
coldata$files = file.path(root, coldata$sample, "quant.sf")
coldata = coldata[c('sample','sex','AgeDays', 'files')]
names(coldata) <- c('names','sex','age','files')
coldata = coldata[file.exists(coldata$files),]
coldata$sex <- factor(coldata$sex, levels=c("XX","XY"))

coldata$salmon_folder=file.path(root, coldata$names)

#Import annotations
###TODO/NOTE: Probably the way to go about this is to load each allelic count with own gtf, then cbind() the resulting SummarizedExperiments together.
refgtffile="/mnt/sas0/AD/lalli/nf_stage/genome_refs/T2T-CHM13_v2_ncbi110/NCBI110_Leung_T2T.gtf"
gtfPath="/mnt/sas0/AD/lalli/reference_comparison_results/T2Tv2_leung_ncbi110_personalized_oldcalls/personalized_transcriptomes/HSB105/HSB105.gtf.gz"
testgtf='/mnt/sas0/AD/lalli/reference_comparison_results/T2Tv2_leung_ncbi110_personalized_oldcalls/test_105.gtf'
library(Rgb)
gtf = read.gtf(testgtf)

se <- importAllelicCounts(coldata2, 'L','R', tx2gene=gene2tx[c('ASE_tx','ASE_gene')])

refgtf <- read.gtf(refgtffile)
gene2tx=generate_tx2gene_from_gtf_obj(refgtf)
gene2tx=make_gene2tx_diploid(gene2tx)
all_txs = gene2tx$ASE_tx
gene2tx[,c('ASE_tx', 'ASE_gene')]

make_gene2tx_diploid <- function(gene2tx){
  gene2tx_L=gene2tx
  gene2tx_R=gene2tx
  gene2tx_L$ASE_tx = paste0(gene2tx$transcript_id, '_L')
  gene2tx_R$ASE_tx = paste0(gene2tx$transcript_id, '_R')
  gene2tx_L$ASE_gene = paste0(gene2tx$gene_id, '_L')
  gene2tx_R$ASE_gene = paste0(gene2tx_R$gene_id, '_R')
  gene2tx_L$ASE_gene = paste0(gene2tx_L$gene_id, '_L')
  gene2tx = gdata::interleave(gene2tx_L, gene2tx_R)
  gene2tx = gene2tx[c('ASE_tx','ASE_gene','transcript_id','gene_id')]
  return (gene2tx)
}

fill_in_missing_transcripts <- function(salmon_folder, all_txs){
  edit_minfo(salmon_folder, all_txs)
  fix_quant(file.path(salmonpath, 'quant.sf'), all_txs)
  fix_bootstrap(salmon_folder, all_txs)
}
HSB105_indexdir='/mnt/sas0/AD/lalli/reference_comparison_results/T2Tv2_leung_ncbi110_personalized_oldcalls/salmon_indexes/HSB105_salmon_indexes'
HSB105_gtf='/mnt/sas0/AD/lalli/reference_comparison_results/T2Tv2_leung_ncbi110_personalized_oldcalls/personalized_transcriptomes/HSB105/HSB105_no_prefixes.gtf.gz'
HSB105_txome='/mnt/sas0/AD/lalli/reference_comparison_results/T2Tv2_leung_ncbi110_personalized_oldcalls/personalized_transcriptomes/HSB105/HSB105.transcriptome.fa.gz'

HSB121_indexdir='/mnt/sas0/AD/lalli/reference_comparison_results/T2Tv2_leung_ncbi110_personalized_oldcalls/salmon_indexes/HSB121_salmon_indexes'
HSB121_gtf='/mnt/sas0/AD/lalli/reference_comparison_results/T2Tv2_leung_ncbi110_personalized_oldcalls/personalized_transcriptomes/HSB121/HSB121_no_prefixes.gtf.gz'
HSB121_txome='/mnt/sas0/AD/lalli/reference_comparison_results/T2Tv2_leung_ncbi110_personalized_oldcalls/personalized_transcriptomes/HSB121/HSB121.transcriptome.fa.gz'
##TODO: if this becomes part of my workflow, tximeta uses the gtf name as the index name.
##So every iteration of the gtf will need to have a unique name.

create_linked_tx <- function(index_dir, txome_fasta, gtf){
  jsonFile <- file.path(index_dir, "tximeta.json")
  makeLinkedTxome(indexDir=index_dir,
                  source="RefSeq+Leung2021", organism="Homo sapiens",
                  release="110", genome="CHM13.v2",
                  fasta=txome_fasta, gtf=gtf,
                  jsonFile=jsonFile)
  return (jsonFile)
}

edit_minfo <- function(salmon_folder, all_txs){
  minfo <- jsonlite::fromJSON(file.path(salmon_folder, 'aux_info', "meta_info.json"))
  minfo$num_valid_targets_originally = minfo$num_valid_targets
  minfo$num_valid_targets = length(all_txs)
  minfo$num_targets = length(all_txs)
  write(jsonlite::toJSON(minfo), file.path(salmon_folder, 'aux_info', "meta_info.json"))
}

fix_quant <- function(file, all_txs){
  all_txs_df = as.data.frame(all_txs)
  colnames(all_txs_df) = c('Name')
  quant <- read.csv(file, sep='\t')
  if (sum(duplicated(quant$Name)) > 0){
    print ("Warning! Duplicate transcripts detected in quant.sf! Removing all but first entry.")
    quant = quant[!duplicated(quant$Name),]
  }
  quant <- left_join(all_txs_df, quant)
  quant[is.na(quant)] <- 0
  write.table(quant, file, sep='\t', row.names=FALSE, quote=FALSE)
  return (quant)
}

fix_bootstrap <- function(salmonpath, all_txs){
  boots <- read_bootstrap(salmonpath)
  all_txs_df = as.data.frame(all_txs)
  colnames(all_txs_df) = c('Name')
  boots <- left_join(all_txs_df, boots, on="Name")
  boots[is.na(boots)] <- 0
  boots <- boots[,-1]
  write_bootstrap(salmonpath, boots) #as.list(all_txs_df$Name)
  fwrite(as.list(all_txs), file.path(salmonpath, 'aux_info', 'bootstrap','names.tsv.gz'), sep='\t')
}

write_bootstrap <- function(salmonpath, atomic.boots){
  if (!is.vector(atomic.boots)){
    atomic.boots <- c(data.matrix(atomic.boots))
  }
  bootCon <- gzcon(file(file.path(salmonpath, 'aux_info','bootstrap', 'bootstraps.gz'), "wb"))
  writeBin(atomic.boots, bootCon)
  close(bootCon)
}

read_bootstrap <- function(salmonpath){
  minfo <- jsonlite::fromJSON(file.path(salmonpath, 'aux_info', "meta_info.json"))
  txnames <- colnames(fread(file.path(salmonpath, 'aux_info','bootstrap', 'names.tsv.gz')))

  bootCon <- gzcon(file(file.path(salmonpath, 'aux_info','bootstrap', 'bootstraps.gz'), "rb"))
  expected.n <- length(txnames) * minfo$num_bootstraps
  boots <- tryCatch({
    bootsIn <- readBin(bootCon, "double", n = expected.n)
    stopifnot(length(bootsIn) == expected.n)
    bootsIn
  }) 
  close(bootCon)
  dim(boots) <- c(length(txnames), minfo$num_bootstraps)
  boots = as.data.frame(boots)
  boots$Name <- txnames
  return (boots)
}

generate_tx2gene <- function(gtf_file){
  gtf = read.gtf(gtf_file)
  gene2tx=unique(gtf[c('transcript_id','gene_id')])
  gene2tx=gene2tx[complete.cases(gene2tx),]
  gene2tx[,2] = gsub('gene-','',gene2tx[,2])
  gene2tx[,1] = gsub('rna-','',gene2tx[,1])
  gene2tx = gene2tx[!duplicated(gene2tx[,1]),]
  return (gene2tx)
}

generate_tx2gene_from_gtf_obj <- function(gtf){
  gene2tx = gtf[gtf$feature=='transcript',]
  gene2tx=unique(gene2tx[c('transcript_id','gene_id')])
  gene2tx=gene2tx[complete.cases(gene2tx),]
  gene2tx[,2] = gsub('gene-','',gene2tx[,2])
  gene2tx[,1] = gsub('rna-','',gene2tx[,1])
  gene2tx = gene2tx[!duplicated(gene2tx[,1]) | ((duplicated(gene2tx[,1]) | duplicated(gene2tx[,1], fromLast=TRUE)) & !grepl('_', gene2tx$gene_id)),]
  return (gene2tx)
}