#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(Seurat); library(ggplot2)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: report_counts.R <pool> <seurat_rds> <demux_best> <final_labels> <out_tsv> <out_png>")
}
pool <- args[1]; seurat_rds <- args[2]; demux_best <- args[3]; final_labels <- args[4]; out_tsv <- args[5]; out_png <- args[6]

# demuxlet counts
Demux <- fread(demux_best)
demux_total <- nrow(Demux)
demux_sng   <- Demux[DROPLET.TYPE == "SNG", .N]
demux_dbl   <- Demux[DROPLET.TYPE == "DBL", .N]

# Seurat QC counts
seu <- readRDS(seurat_rds)
dt  <- as.data.table(seu@meta.data)
qc_total <- nrow(dt)
dt[, qc_keep := (doublet_class == "singlet" & nuclear_fraction >= 0.4 & percent.mt <= 20)]
qc_keep_n <- dt[qc_keep == TRUE, .N]

# Intersection (already computed by qc_intersect)
Final <- fread(final_labels)
intersect_n <- nrow(Final)

# Percentages
pct <- function(a,b) ifelse(b>0, 100*a/b, 0)

summary_dt <- data.table(
  pool = pool,
  demux_total = demux_total,
  demux_sng = demux_sng,
  demux_dbl = demux_dbl,
  demux_sng_pct = pct(demux_sng, demux_total),
  demux_dbl_pct = pct(demux_dbl, demux_total),
  qc_total = qc_total,
  qc_keep = qc_keep_n,
  qc_keep_pct = pct(qc_keep_n, qc_total),
  intersect_n = intersect_n,
  intersect_vs_demux_pct = pct(intersect_n, demux_sng),
  intersect_vs_qc_pct = pct(intersect_n, qc_keep_n)
)

fwrite(summary_dt, out_tsv, sep="\t")

# Bar plot
plot_dt <- data.table(
  category = c("Demux singlets", "Demux doublets", "QC keep", "QC∩Demux singlets"),
  count    = c(demux_sng, demux_dbl, qc_keep_n, intersect_n)
)

g <- ggplot(plot_dt, aes(x=category, y=count, fill=category)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title=paste0("Counts for ", pool), x="", y="Count") +
  guides(fill="none")

ggsave(out_png, g, width=6, height=4, dpi=150)

cat(sprintf("[INFO] Wrote summary for %s: %s, %s\n", pool, out_tsv, out_png))
