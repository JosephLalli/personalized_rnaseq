#!/usr/bin/env Rscript

library(SummarizedExperiment)
library(tximport)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Usage: salmon_tximport.r <coldata> <salmon_out>", call.=FALSE)
}

coldata = args[1]
path = args[2]
sample_name = args[3]

prefix = sample_name
tx2gene = "salmon_tx2gene.tsv"
info = file.info(tx2gene)
if (info$size == 0) {
    tx2gene = NULL
} else {
    rowdata = read.csv(tx2gene, sep="\t", header = FALSE)
    colnames(rowdata) = c("tx", "gene_id", "gene_name")
    tx2gene = rowdata[,1:2]
}

fns = list.files(path, pattern = "quant.sf", recursive = T, full.names = T)
names = basename(dirname(fns))
names(fns) = names

if (file.exists(coldata)) {
    coldata = read.csv(coldata, sep="\t")
    coldata = coldata[match(names, coldata[,1]),]
    coldata = cbind(files = fns, coldata)
} else {
    message("ColData not avaliable ", coldata)
    coldata = data.frame(files = fns, names = names)
}

txi = tximport(fns, type = "salmon", txOut = TRUE)
rownames(coldata) = coldata[["names"]]
extra = setdiff(rownames(txi[[1]]),  as.character(rowdata[["tx"]]))
if (length(extra) > 0) {
    rowdata = rbind(rowdata, data.frame(tx=extra, gene_id=extra, gene_name=extra))
}
rowdata = rowdata[match(rownames(txi[[1]]), as.character(rowdata[["tx"]])),]
rownames(rowdata) = rowdata[["tx"]]
se = SummarizedExperiment(assays = list(counts = txi[["counts"]], abundance = txi[["abundance"]], length = txi[["length"]]),
                        colData = DataFrame(coldata),
                        rowData = rowdata)
if (!is.null(tx2gene)) {
    gi = summarizeToGene(txi, tx2gene = tx2gene)
    gi.ls = summarizeToGene(txi, tx2gene = tx2gene,countsFromAbundance="lengthScaledTPM")
    gi.s = summarizeToGene(txi, tx2gene = tx2gene,countsFromAbundance="scaledTPM")
    growdata = unique(rowdata[,2:3])
    growdata = growdata[match(rownames(gi[[1]]), growdata[["gene_id"]]),]
    rownames(growdata) = growdata[["tx"]] # In my hands, resolves to NULL, but that's a good thing
    gse = SummarizedExperiment(assays = list(counts = gi[["counts"]], abundance = gi[["abundance"]], length = gi[["length"]]),
                                colData = DataFrame(coldata),
                                rowData = growdata)
    gse.ls = SummarizedExperiment(assays = list(counts = gi.ls[["counts"]], abundance = gi.ls[["abundance"]], length = gi.ls[["length"]]),
                                colData = DataFrame(coldata),
                                rowData = growdata)
    gse.s = SummarizedExperiment(assays = list(counts = gi.s[["counts"]], abundance = gi.s[["abundance"]], length = gi.s[["length"]]),
                                colData = DataFrame(coldata),
                                rowData = growdata)
}

txi.dtu = tximport(fns, type = "salmon", tx2gene=tx2gene, txOut = TRUE, countsFromAbundance="dtuScaledTPM")
se.dtu = SummarizedExperiment(assays = list(counts = txi.dtu[["counts"]], abundance = txi.dtu[["abundance"]], length = txi.dtu[["length"]]),
                        colData = DataFrame(coldata),
                        rowData = rowdata)

if ("infReps" %in% names(txi)){
    txi.ls.infReps = tximport(fns, type = "salmon", tx2gene=tx2gene, txOut = TRUE, countsFromAbundance="dtuScaledTPM", varReduce=TRUE, infRepStat=matrixStats::rowMedians)
    gi.ls.infReps = summarizeToGene(txi.ls.infReps, tx2gene = tx2gene, countsFromAbundance="lengthScaledTPM", varReduce=TRUE)
    # gi.ls.infReps = summarizeToGene(txi, tx2gene = tx2gene,countsFromAbundance="lengthScaledTPM", varReduce=TRUE)
    se.ls.infReps = SummarizedExperiment(assays = list(counts = txi.ls.infReps[["counts"]], abundance = txi.ls.infReps[["abundance"]], length = txi.ls.infReps[["length"]]),
                                colData = DataFrame(coldata),
                                rowData = rowdata)
    gse.ls.infReps = SummarizedExperiment(assays = list(counts = gi.ls.infReps[["counts"]], abundance = gi.ls.infReps[["abundance"]], length = gi.ls.infReps[["length"]]),
                                colData = DataFrame(coldata),
                                rowData = growdata)

}
# 
# txi.ls = tximport(fns, type = "salmon", tx2gene=tx2gene, txOut = TRUE, countsFromAbundance="lengthScaledTPM", varReduce=TRUE, infRepStat=matrixStats::rowMedians)


build_table = function(se.obj, slot) {
    cbind(rowData(se.obj)[,1:2], assays(se.obj)[[slot]])
}

if(exists("gse")){
    write.table(build_table(gse, "abundance"), paste(c(prefix, "gene_tpm.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse, "counts"), paste(c(prefix, "gene_counts.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse, "length"), paste(c(prefix, "gene_length.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse.ls, "abundance"), paste(c(prefix, "gene_tpm_length_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse.ls, "counts"), paste(c(prefix, "gene_counts_length_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse.ls, "length"), paste(c(prefix, "gene_length_length_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse.s, "abundance"), paste(c(prefix, "gene_tpm_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse.s, "counts"), paste(c(prefix, "gene_counts_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse.s, "length"), paste(c(prefix, "gene_length_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
}

write.table(build_table(se,"abundance"), paste(c(prefix, "transcript_tpm.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(build_table(se, "counts"), paste(c(prefix, "transcript_counts.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(build_table(se, "length"), paste(c(prefix, "transcript_length.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(build_table(se.dtu, "abundance"), paste(c(prefix, "transcript_tpm_dtuscaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(build_table(se.dtu, "counts"), paste(c(prefix, "transcript_counts_dtuscaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(build_table(se.dtu, "length"), paste(c(prefix, "transcript_length_dtuscaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)

if (exists("gse.ls.infReps")){
    write.table(build_table(gse.ls.infReps, "abundance"), paste(c(prefix, "gene_tpm_length_scaled_median_bootstraps.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse.ls.infReps, "counts"), paste(c(prefix, "gene_counts_length_scaled_median_bootstraps.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse.ls.infReps, "length"), paste(c(prefix, "gene_length_length_scaled_median_bootstraps.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(se.ls.infReps,"abundance"), paste(c(prefix, "transcript_tpm_length_scaled_median_bootstraps.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(se.ls.infReps, "counts"), paste(c(prefix, "transcript_counts_length_scaled_median_bootstraps.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(se.ls.infReps, "length"), paste(c(prefix, "transcript_length_length_scaled_median_bootstraps.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
}

saveRDS(se,     file = paste0("transcript_tpm.tximport_se.rds"))
saveRDS(se.dtu,   file = paste0("transcript_tpm_dtuscaled.tximport_se.rds"))
saveRDS(se.ls.infReps,   file = paste0("transcript_tpm_length_scaled.median_bootstrap.tximport_se.rds"))
saveRDS(gse,    file = paste0("gene_tpm.tximport_se.rds"))
saveRDS(gse.ls, file = paste0("gene_tpm_length_scaled.tximport_se.rds"))
saveRDS(gse.s,  file = paste0("gene_tpm_scaled.tximport_se.rds"))
saveRDS(gse.ls.infReps,  file = paste0("gene_tpm_length_scaled.median_bootstrap.tximport_se.rds"))

# Print sessioninfo to standard out
citation("tximeta")
sessionInfo()
