suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

setwd("<SELF_NONSELF>/human/RNAseq/HCC1937")

# ---- Paths (edit if needed) ----
fc_file  <- "featureCount/HCC1937_featureCounts_peptide_table.tsv"
kal_file <- "kalisto/HCC1937_leadingIsoform_9mer_table.tsv"
out_prefix <- "HCC1937_unique9mers_per_gene_canonical_vs_leading"

# ---- Read ----
fc  <- fread(fc_file)
kal <- fread(kal_file)

stopifnot(all(c("gene","peptide_9mer") %in% names(fc)))
stopifnot(all(c("gene","peptide_9mer") %in% names(kal)))

# ---- Per-gene counts of UNIQUE peptides ----
fc_gene <- fc[, .(n9_fc = uniqueN(peptide_9mer)), by=gene]
# Processed 15350 groups out of 15350
kal_gene <- kal[, .(n9_kal = uniqueN(peptide_9mer)), by=gene]
# Processed 20076 groups out of 20076.

# merge on gene (keep genes present in either model)
m <- merge(fc_gene, kal_gene, by="gene", all=TRUE)
m[is.na(n9_fc),  n9_fc := 0L]
m[is.na(n9_kal), n9_kal := 0L]

# Derived metrics
m[, delta := n9_kal - n9_fc]
m[, ratio := (n9_kal + 1) / (n9_fc + 1)]  # +1 to avoid divide-by-zero blowups

# ---- Quick summary ----
cat("Genes in canonical (fc):", nrow(fc_gene), "\n") # 15350
cat("Genes in leading (kal):", nrow(kal_gene), "\n") # 20 076
cat("Genes with MORE unique 9-mers in leading vs canonical:", sum(m$delta > 0), "\n") # 7201
cat("Genes with SAME:", sum(m$delta == 0), "\n") # 9200
cat("Genes with FEWER:", sum(m$delta < 0), "\n") # 3675
cat("Median delta (kal - fc):", median(m$delta), "\n") # 0

fwrite(m, paste0(out_prefix, "_table.tsv"), sep="\t")

# ---- Plot 1: distribution of unique 9-mers per gene (two models) ----
long <- rbind(
  fc_gene[, .(gene, n9 = n9_fc, model="Canonical (featureCounts)")],
  kal_gene[, .(gene, n9 = n9_kal, model="Leading isoform (kallisto)")]
)

p1 <- ggplot(long, aes(x=n9, fill=model)) +
  geom_histogram(bins=80, position="identity", alpha=0.45) +
  scale_x_log10() +
  labs(
    title="HCC1937: unique 9-mers per gene",
    x="Number of unique 9-mers per gene (log10 scale)",
    y="Number of genes",
    fill=""
  ) +
  theme_bw()

ggsave(paste0(out_prefix, "_hist_logx.png"), p1, width=8, height=5, dpi=300)

# ---- Plot 2: scatter leading vs canonical per gene ----
p2 <- ggplot(m, aes(x=n9_fc, y=n9_kal)) +
  geom_point(alpha=0.25, size=0.8) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(
    title="HCC1937: unique 9-mers per gene (leading isoform vs canonical)",
    x="Canonical model: # unique 9-mers per gene (log10)",
    y="Leading isoform model: # unique 9-mers per gene (log10)"
  ) +
  theme_bw()

ggsave(paste0(out_prefix, "_scatter_loglog.png"), p2, width=6.5, height=6.5, dpi=300)

# ---- Plot 3 (optional): delta distribution ----
p3 <- ggplot(m, aes(x=delta)) +
  geom_histogram(bins=80) +
  labs(
    title="HCC1937: change in unique 9-mers per gene (leading - canonical)",
    x="Δ unique 9-mers per gene",
    y="Number of genes"
  ) +
  theme_bw()

ggsave(paste0(out_prefix, "_delta_hist.png"), p3, width=8, height=5, dpi=300)
