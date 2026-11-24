#!/usr/bin/env Rscript
# Compare 9-mer repertoires between JY and HCC1937:
#  - counts for shared / JY-only / HCC-only
#  - distributions of cumulative_est_counts for each category
#  - (optional) log10 distributions
#
# Input files (as provided):
#   HCC1937/kalisto/peptide_ninemers_ExprSum.tsv
#   JY/peptide_ninemers_ExprSum.tsv
#
# Expected columns:
#   ninemer, IDs, cumulative_est_counts
#
# Usage:
#   Rscript compare_peptide_repertoires.R
#
# If you want to change thresholds or output names, edit parameters below.

setwd("<SELF_NONSELF>/human/RNAseq/")

suppressWarnings({
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
  })
})

# ---------------------------
# Parameters
# ---------------------------
file_hcc <- "HCC1937/kalisto/peptide_9mers_ExprSum.tsv"
file_jy  <- "JY/peptide_9mers_ExprSum.tsv"

# Define whether a peptide is "present"
# If you want to require a minimum cumulative value, set > 0 to e.g. >= 1 or >= 10.
presence_threshold <- 0

# Pseudocount for log10 plots (avoid log10(0))
log_pseudocount <- 1e-6

# Output prefix
out_prefix <- "JY_vs_HCC1937_ninemer"

# ---------------------------
# Load data
# ---------------------------
read_peptide_tsv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  
  dt <- fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
  
  required <- c("ninemer", "cumulative_est_counts")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Missing required columns in ", path, ": ", paste(missing, collapse = ", "))
  }
  
  # Ensure correct types
  dt[, ninemer := as.character(ninemer)]
  dt[, cumulative_est_counts := as.numeric(cumulative_est_counts)]
  dt
}

jy  <- read_peptide_tsv(file_jy)
hcc <- read_peptide_tsv(file_hcc)

# If a 9-mer appears more than once (shouldn't, but safe), collapse by sum
jy  <- jy[,  .(cumulative_est_counts = sum(cumulative_est_counts, na.rm = TRUE)), by = ninemer]
hcc <- hcc[, .(cumulative_est_counts = sum(cumulative_est_counts, na.rm = TRUE)), by = ninemer]

# Presence sets
jy_present  <- jy[cumulative_est_counts > presence_threshold, ninemer]
hcc_present <- hcc[cumulative_est_counts > presence_threshold, ninemer]

set_jy  <- unique(jy_present)
set_hcc <- unique(hcc_present)

shared   <- intersect(set_jy, set_hcc)
jy_only  <- setdiff(set_jy, set_hcc)
hcc_only <- setdiff(set_hcc, set_jy)

cat("\n=== Repertoire counts (threshold > ", presence_threshold, ") ===\n", sep = "")
cat("Shared peptides:      ", length(shared), "\n", sep = "") # 8691361
cat("Unique to JY:         ", length(jy_only), "\n", sep = "") # 1009318 = 10.4 %
cat("Unique to HCC1937:    ", length(hcc_only), "\n", sep = "") # 773172 = 8.2 %

# ---------------------------
# Build a table of values by category
# ---------------------------
# Merge values (keep both columns; missing -> 0)
merged <- merge(jy, hcc, by = "ninemer", all = TRUE, suffixes = c("_JY", "_HCC"))
merged[is.na(cumulative_est_counts_JY),  cumulative_est_counts_JY := 0]
merged[is.na(cumulative_est_counts_HCC), cumulative_est_counts_HCC := 0]

merged[, category := fifelse(
  ninemer %in% shared, "Shared",
  fifelse(ninemer %in% jy_only, "JY-only",
          fifelse(ninemer %in% hcc_only, "HCC1937-only", "Other"))
)]

# Keep only the three categories of interest (should cover all "present" peptides)
vals <- merged[category %in% c("Shared", "JY-only", "HCC1937-only")]

# For distributions, use the cumulative value from the relevant sample:
# - Shared: you may want both; we will create two "shared" distributions (JY and HCC) and also a combined one.
vals_long <- rbind(
  vals[category == "Shared",      .(category = "Shared (JY value)",  value = cumulative_est_counts_JY)],
  vals[category == "Shared",      .(category = "Shared (HCC value)", value = cumulative_est_counts_HCC)],
  vals[category == "JY-only",     .(category = "JY-only",            value = cumulative_est_counts_JY)],
  vals[category == "HCC1937-only",.(category = "HCC1937-only",       value = cumulative_est_counts_HCC)]
)

# Drop non-positive values (shouldn't exist except for safety)
vals_long <- vals_long[value > presence_threshold]

# Add log10 value
vals_long[, log10_value := log10(value + log_pseudocount)]

# ---------------------------
# Plotting helpers
# ---------------------------
save_hist <- function(dt, xcol, outfile, title, xlab) {
  p <- ggplot(dt, aes(x = .data[[xcol]])) +
    geom_histogram(bins = 80) +
    facet_wrap(~ category, scales = "free_y") +
    labs(title = title, x = xlab, y = "Count of peptides") +
    theme_bw()
  
  ggsave(outfile, p, width = 10, height = 6, dpi = 300)
  cat("Wrote: ", outfile, "\n", sep = "")
}

# ---------------------------
# Venn counts plot (simple bar summary)
# ---------------------------
count_dt <- data.table(
  group = c("Shared", "JY-only", "HCC1937-only"),
  n = c(length(shared), length(jy_only), length(hcc_only))
)

p_counts <- ggplot(count_dt, aes(x = group, y = n)) +
  geom_col() +
  labs(
    title = "9-mer repertoire overlap: JY vs HCC1937",
    x = "",
    y = "Number of 9-mers"
  ) +
  theme_bw()

p_counts
# counts_png <- paste0(out_prefix, "_venn_counts.png")
# ggsave(counts_png, p_counts, width = 7, height = 5, dpi = 300)
# cat("Wrote: ", counts_png, "\n", sep = "")

# ---------------------------
# Distributions (raw scale)
# ---------------------------
hist_raw_png <- paste0(out_prefix, "_cumulative_counts_hist.png")
save_hist(
  vals_long,
  xcol = "value",
  outfile = hist_raw_png,
  title = "Distributions of cumulative_est_counts by repertoire category",
  xlab = "cumulative_est_counts"
)

# ---------------------------
# Distributions (log10 scale)
# ---------------------------
hist_log_png <- paste0(out_prefix, "_log10_cumulative_counts_hist.png")
save_hist(
  vals_long,
  xcol = "log10_value",
  outfile = hist_log_png,
  title = "Distributions of log10(cumulative_est_counts + pseudocount) by category",
  xlab = "log10(cumulative_est_counts + pseudocount)"
)

# ---------------------------
# Optional: write category membership tables
# ---------------------------
# fwrite(data.table(ninemer = shared),   paste0(out_prefix, "_shared_ninemers.tsv"),   sep = "\t")
# fwrite(data.table(ninemer = jy_only),  paste0(out_prefix, "_jy_only_ninemers.tsv"),  sep = "\t")
# fwrite(data.table(ninemer = hcc_only), paste0(out_prefix, "_hcc_only_ninemers.tsv"), sep = "\t")

################### HCC1937 ########################
setwd("HCC1937/")
fC <- read.csv("featureCount/fc_geneid_length_counts_ProportionWithMMs_RPK_TPM.tsv", sep="\t")
kal <- read.csv("kalisto/abundance.tsv", sep="\t")
fcFilt <- fC$Length_kb[which(fC$TPM >=1 )]
kalFilt <- kal$length[which(kal$tpm >=1)]/1000   
  
par(mfrow=c(2,1), mar=c(3.5,3.5,1.5,1), cex.main=1, cex.lab=0.9, cex.axis=0.9, mgp=c(1.5, 0.25,0))
hist(log10(fcFilt), main="14 167 genes (TPM>=1)", xlab="log10-scaled length (kb)", ylab="count", ylim=c(0,3000))
hist(log10(kalFilt), main="33 051 isoforms (TPM>=1)", xlab="log10-scaled length (kb)", ylab="count", ylim=c(0,8000))

summary(fcFilt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.041   1.882   3.549   4.793   6.159  91.667 
summary(kalFilt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.039   0.714   1.212   1.589   1.971  46.050 

# Fixed breaks shared across both plots
br <- seq(log10(0.03), log10(100), by = 0.05)  # 30 bp .. 100 kb in log10(kb)

par(mfrow=c(2,1), mar=c(3.5,3.8,1.5,1), mgp=c(1.8, 0.6, 0), cex.axis=0.9)

# Genes (featureCounts)
hist(log10(fcFilt), breaks=br,
     main="featureCounts: genes (TPM >= 1)",
     xlab="length (kb)", ylab="count",
     xaxt="n")
# Isoforms (kallisto)
hist(log10(kalFilt), breaks=br,
     main="kallisto: isoforms (TPM >= 1)",
     xlab="length (kb)", ylab="count",
     xaxt="n")

# Add human-readable tick marks in kb
ticks_kb <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)
axis(1, at=log10(ticks_kb), labels=ticks_kb)

################# Distribution of expressed isoforms per gene (TPM > 1) ##################
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- Inputs ----
abund <- "<SELF_NONSELF>/human/RNAseq/HCC1937/kalisto/abundance.tsv"  # adjust path
mapfile <- "<SELF_NONSELF>/human/References/GRCh38.p14/CDS_GRCh38p14_annotNames.csv"

tpm_thr <- 1

# ---- Read kallisto ----
kal <- fread(abund)
# expected columns: target_id, length, eff_length, est_counts, tpm
stopifnot(all(c("target_id","tpm") %in% names(kal)))

# Extract RefSeq CDS protein accession from target_id
# Example: lcl|NC_000001.11_cds_NP_001005484.2_1  -> NP_001005484.2
kal[, CDS := sub("^.*_cds_", "", target_id)]
kal[, CDS := sub("_[^_]+$", "", CDS)]  # drop last underscore suffix (e.g. _1)

# Keep expressed isoforms
kal_exp <- kal[tpm >= tpm_thr, .(CDS, tpm)]

# ---- Read mapping CDS -> gene ----
map <- fread(mapfile)
# columns: CDS, gene, xref
stopifnot(all(c("CDS","gene") %in% names(map)))

# Join to gene
kal_exp <- merge(kal_exp, map[, .(CDS, gene)], by="CDS", all.x=TRUE)

# Some CDS might not map (rare); keep track
cat("Expressed isoforms (TPM>=1):", nrow(kal_exp), "\n") # 39908 
cat("Missing gene mapping:", sum(is.na(kal_exp$gene)), "\n") # 8

kal_exp2 <- kal_exp[!is.na(gene)]

# ---- Isoforms per gene distribution ----
iso_per_gene <- kal_exp2[, .(
  n_isoforms = uniqueN(CDS),
  tpm_sum = sum(tpm)
), by = gene]

cat("Expressed genes (>=1 isoform with TPM>=1):", nrow(iso_per_gene), "\n") # 11559 
cat("Median isoforms per expressed gene:", median(iso_per_gene$n_isoforms), "\n") # 2
cat("Genes with >=2 isoforms:", sum(iso_per_gene$n_isoforms >= 2), "\n") # 6800
cat("Genes with >=5 isoforms:", sum(iso_per_gene$n_isoforms >= 5), "\n") # 1704

# Plot: histogram of isoforms per gene
p1 <- ggplot(iso_per_gene, aes(x = n_isoforms)) +
  geom_histogram(binwidth = 1, boundary = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    title = paste0("Isoforms per gene (kallisto TPM >= ", tpm_thr, ")"),
    x = "Number of expressed isoforms per gene",
    y = "Number of genes"
  ) +
  theme_bw()

ggsave("HCC1937_isoforms_per_gene_hist.png", p1, width=8, height=5, dpi=300)

# Optional: show heavy tail on log y-scale
p2 <- p1 + scale_y_log10() + labs(y="Number of genes (log10)")
ggsave("HCC1937_isoforms_per_gene_hist_logy.png", p2, width=8, height=5, dpi=300)

# Optional: fraction of genes with n_isoforms >= k (a CDF-like summary)
ks <- c(1,2,3,5,10)
frac <- data.table(k=ks, frac=sapply(ks, function(k) mean(iso_per_gene$n_isoforms >= k)))
print(frac)
# k       frac  
# 1 1.00000000
# 2 0.58828618
# 3 0.35954667
# 5 0.14741760
# 10 0.02889523
# fwrite(frac, "HCC1937_isoforms_per_gene_fraction_ge_k.tsv", sep="\t")

library(data.table)
library(ggplot2)

# iso_per_gene already computed earlier
genes_2iso <- iso_per_gene[n_isoforms == 2, gene]

kal_2iso <- kal_exp2[gene %in% genes_2iso]

# Compute per-gene proportions
kal_2iso[, total_tpm := sum(tpm), by = gene]
kal_2iso[, prop := tpm / total_tpm]

# Sort within gene by decreasing TPM
kal_2iso <- kal_2iso[order(gene, -prop)]
kal_2iso[, iso_rank := 1:.N, by = gene]

# Keep rank 1 and 2
kal_2iso <- kal_2iso[iso_rank %in% c(1,2)]

cat("Number of genes with exactly 2 isoforms:", length(unique(kal_2iso$gene)), "\n") # 2644

# ---- Plot stacked bar (percentage per gene) ----
p_stacked <- ggplot(kal_2iso, aes(x = gene, y = prop, fill = factor(iso_rank))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("salmon","cornflowerblue"),
                    labels=c("Major isoform","Minor isoform")) +
  labs(
    title="Genes with exactly 2 expressed isoforms (TPM ≥ 1)",
    x="Genes",
    y="Fraction of total gene TPM",
    fill=""
  ) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("HCC1937_two_isoform_genes_stacked.png",
       p_stacked, width=10, height=5, dpi=300)


###
major_iso <- kal_2iso[iso_rank == 1]

p_dist <- ggplot(major_iso, aes(x = prop)) +
  geom_histogram(bins=40, fill="cornflowerblue") +
  labs(
    title="Dominance of major isoform (genes with exactly 2 isoforms)",
    x="Fraction of total gene TPM contributed by major isoform",
    y="Number of genes"
  ) +
  theme_bw()

ggsave("HCC1937_two_isoform_major_fraction.png",
       p_dist, width=7, height=5, dpi=300)

mean(major_iso$prop > 0.9) # 0.1626324
mean(major_iso$prop < 0.7) # 0.5798033
mean(major_iso$prop < 0.5) # 0.05295008
mean(major_iso$prop > 0.5) # 0.7246596
