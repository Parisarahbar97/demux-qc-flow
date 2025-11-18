#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(Seurat)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: qc_intersect.R <seurat_rds> <demux_best> <out_tsv>")
}
seurat_rds <- args[1]
demux_best <- args[2]
out_tsv    <- args[3]

seu <- readRDS(seurat_rds)
dt  <- as.data.table(seu@meta.data, keep.rownames = TRUE)
setnames(dt, "rn", "barcode")
dt[, qc_keep := (doublet_class == "singlet" & nuclear_fraction >= 0.4 & percent.mt <= 20)]
qc_sing <- dt[qc_keep == TRUE, .(barcode)]

demux <- fread(demux_best)
demux_sng <- demux[DROPLET.TYPE == "SNG",
                   .(barcode = BARCODE,
                     Individual_ID = SNG.BEST.GUESS,
                     demux_NUM.SNPS = NUM.SNPS,
                     demux_RD.UNIQ = RD.UNIQ,
                     demux_SNG.POST = SNG.POSTERIOR)]

res <- merge(qc_sing, demux_sng, by = "barcode")
res[, Consensus_Tier := "QC∩DEMUX_SNG"]
fwrite(res, out_tsv, sep = "\t")
cat(sprintf("[INFO] Wrote %d consensus singlets to %s\n", nrow(res), out_tsv))
