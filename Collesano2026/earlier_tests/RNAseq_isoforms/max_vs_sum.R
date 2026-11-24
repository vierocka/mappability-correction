setwd("<SELF_NONSELF>/human/RNAseq")

JYec <- read.table("JY/gene_sum_max_estcounts.tsv", sep="\t", header = TRUE)
HCT116ec <- read.table("HCT116/gene_sum_max_estcounts.tsv", sep="\t", header = TRUE)
HCC1937ec <- read.table("HCC1937/kalisto/gene_sum_max_estcounts.tsv", sep="\t", header = TRUE)

JYtpm <- read.table("JY/gene_sum_max_tpm.tsv", sep="\t", header = TRUE)
HCT116tpm <- read.table("HCT116/gene_sum_max_tpm.tsv", sep="\t", header = TRUE)
HCC1937tpm <- read.table("HCC1937/kalisto/gene_sum_max_tpm.tsv", sep="\t", header = TRUE)

# JY est counts
plot(log10(JYec$sum_est_counts+1), log10(JYec$max_est_counts+1), xlab="log10-scaled sum +1", ylab="mlog10-scaled max +1", main="JY, estimated counts: 20 076 genes")
# JY tpm
plot(log10(JYtpm$sum_tpm+1), log10(JYtpm$max_tpm+1), xlab="log10-scaled sum +1", ylab="mlog10-scaled max +1", main="JY, tpm: 20 076 genes")

# HCT116 est counts
plot(log10(HCT116ec$sum_est_counts+1), log10(HCT116ec$max_est_counts+1), xlab="log10-scaled sum +1", ylab="mlog10-scaled max +1", main="HCT116, estimated counts: 20 076 genes")
# HCT116 tpm
plot(log10(HCT116tpm$sum_tpm+1), log10(HCT116tpm$max_tpm+1), xlab="log10-scaled sum +1", ylab="mlog10-scaled max +1", main="HCT116, tpm: 20 076 genes")

# HCC1937 est counts
plot(log10(HCC1937ec$sum_est_counts+1), log10(HCC1937ec$max_est_counts+1), xlab="log10-scaled sum +1", ylab="mlog10-scaled max +1", main="HCC1937, estimated counts: 20 076 genes")
# HCC1937 tpm
plot(log10(HCC1937tpm$sum_tpm+1), log10(HCC1937tpm$max_tpm+1), xlab="log10-scaled sum +1", ylab="mlog10-scaled max +1", main="HCC1937, tpm: 20 076 genes")

### length of max
par(mfrow=c(3,1), mar=c(3,3,1,0.5), mgp=c(1.5,0.5,0), cex.main=1, cex.axis=1, cex.lab=1)
hist(log10(JYec$best_isoform_length), main="JY: length of the most abundant isoform (20 076 genes)", n=20, xlab="log10-scaled length (bp)")
hist(log10(HCC1937ec$best_isoform_length), main="HCC1937: length of the most abundant isoform (20 076 genes)", n=20, xlab="log10-scaled length (bp)")
hist(log10(HCT116ec$best_isoform_length), main="HCT116: length of the most abundant isoform (20 076 genes)", n=20, xlab="log10-scaled length (bp)")

diff <- JYec$best_isoform_length-HCC1937ec$best_isoform_length
myTab <- as.table(cbind(JYec$sum_est_counts, HCC1937ec$sum_est_counts, diff))
head(myTab)
colnames(myTab) <- c("JYecSum","HCC1937ecSum","diffL")

### 6 groups:
## both expressed, diff 0 - black - 11 115
## none expressed - white - 2185
## only JY expressed (HCC1937=0) - cornflowerblue - 1688
## only HCC1937 expressed - goldenrod - 1507
## both both expressed, JY longer - blue - 1911
## both expressed, HCC1937 longer - firebrick - 1670

jy  <- myTab[, "JYecSum"]
hcc <- myTab[, "HCC1937ecSum"]
d   <- myTab[, "diffL"]

colors <- ifelse(jy == 0 & hcc == 0, "white",
                 ifelse(jy > 0 & hcc == 0, "goldenrod",
                        ifelse(jy == 0 & hcc > 0, "cornflowerblue",
                               ifelse(d == 0, "black",
                                      ifelse(d > 0, "blue", "firebrick")))))

table(colors)

#  JY - HCC1937
plot(x=c(1:20076), y=log10(JYec$best_isoform_length)-log10(HCC1937ec$best_isoform_length), xlab="order", ylab="log10-scaled length (bp): JY - HCC1937", col=colors)


length(which(diff > 0))
length(which(diff == 0)) # 15319
length(which(diff < 0))

# isoform length: JY vs HCC1937
plot(log10(JYec$best_isoform_length), log10(HCC1937ec$best_isoform_length), main="JY vs HCC1937: prefered isoform length\n(20 076 genes)", pch=20, col="dimgray", xlab = "JY: log10-scaled length (bp)", ylab="HCC1937: log10-scaled length (bp)")


#######################################################
######## gene level vs isoform level ###############

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

###### HCC1937
# ---- paths ----
fc_file <- "<SELF_NONSELF>/human/RNAseq/HCC1937/featureCount/fc_geneid_length_counts_ProportionWithMMs_RPK_TPM.tsv"
k_gene  <- "<SELF_NONSELF>/human/RNAseq/HCC1937/kalisto/gene_sum_max_estcounts.tsv"

# ---- 1) featureCounts: keep Geneid, Length_kb, Counts ----
fc <- read_tsv(fc_file, col_types = cols()) %>%
  transmute(
    gene = Geneid,
    gene_length_kb = Length_kb,
    fc_counts = Counts
  )

# ---- 2) kallisto gene max: keep gene, max_est_counts, eff_length (col 5) ----
# If your 5th column already is the effective length of the best isoform,
# we read it as 'best_isoform_eff_length' by position.
kg_raw <- read_tsv(k_gene, col_types = cols())

kg <- read_tsv(k_gene, col_types = cols()) %>%
  transmute(
    gene = gene,
    kallisto_max_est_counts = max_est_counts,
    best_isoform_eff_length = best_isoform_length  # if column 5 is protein length
  )

# ---- 3) merge: keep only genes present in kallisto (coding set) ----
merged <- kg %>%
  inner_join(fc, by = "gene")

# ---- 4) write merged table ----
write_tsv(merged, "<SELF_NONSELF>/human/RNAseq/HCC1937/fc_vs_kallisto_max_merged.tsv")

# ---- 5) plot: featureCounts counts vs kallisto max est_counts (log-log) ----
p <- ggplot(merged, aes(x = fc_counts + 1, y = kallisto_max_est_counts + 1)) +
  geom_point(alpha = 0.35, size = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4) +
  labs(
    title = "FeatureCounts gene counts vs\n kallisto gene MAX est_counts (HCC1937)",
    x = "featureCounts Counts (+1), log10",
    y = "kallisto MAX est_counts (+1), log10"
  ) +
  theme_minimal(base_size = 13)

ggsave("fc_vs_kallisto_max_scatter.pdf", p, width = 6.2, height = 5.6)
ggsave("fc_vs_kallisto_max_scatter.png", p, width = 6.2, height = 5.6, dpi = 300)

print(p)
cor(merged$fc_counts, merged$kallisto_max_est_counts, method = "spearman", use = "complete.obs")
# 0.8646236 
# Gene-level featureCounts and isoform-aware kallisto MAX show strong concordance (Spearman ρ ≈ 0.86). Differences are largely confined to low-expression genes, consistent with different handling of multimapping reads.

plot(log10(merged$gene_length_kb), log10(merged$best_isoform_eff_length/1000), main="HCC1937 (20 076 genes)", xlab="log10-scaled gene length (kb)", ylab="log10-scaled the max isoform length (kb)")

# The effective coding length per gene is strongly constrained by isoform usage, 
# not by total genomic span.

###### JY
# ---- paths ----
fc_file <- "<SELF_NONSELF>/human/RNAseq/JY/SRR364065.counts.csv"
k_gene  <- "<SELF_NONSELF>/human/RNAseq/JY/gene_sum_max_estcounts.tsv"

# ---- 1) Read featureCounts output ----
# This file has a comment line starting with '#', then a tab-separated header.
fc <- read_tsv(
  fc_file,
  comment = "#",
  col_types = cols()
)

# The 7th column name is the BAM path; grab it programmatically
count_col <- names(fc)[7]

fc2 <- fc %>%
  transmute(
    gene = Geneid,
    gene_length_kb = Length / 1000,           # Length is in bp
    fc_counts = .data[[count_col]]
  )

# ---- 2) Read kallisto gene table (max isoform length is in aa) ----
kg <- read_tsv(k_gene, col_types = cols()) %>%
  transmute(
    gene = gene,
    best_isoform_length_aa = best_isoform_length
  ) %>%
  mutate(
    best_isoform_length_kbCDS = (best_isoform_length_aa * 3) / 1000  # optional conversion
  )

# ---- 3) Merge (restrict to genes present in kallisto coding set) ----
merged <- kg %>%
  inner_join(fc2, by = "gene")

write_tsv(merged, "<SELF_NONSELF>/human/RNAseq/JY/JY_gene_vs_bestIsoformLength_merged.tsv")

# ---- 4) Plot: gene length vs best isoform length ----
p <- ggplot(merged, aes(x = fc_counts + 1, y = kallisto_max_est_counts + 1)) +
  geom_point(alpha = 0.35, size = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4) +
  labs(
    title = "JY: FeatureCounts gene counts vs\n kallisto gene MAX est_counts",
    x = "featureCounts Counts (+1), log10",
    y = "kallisto MAX est_counts (+1), log10") +
  theme_minimal(base_size = 13)

ggsave("<SELF_NONSELF>/human/RNAseq/JY/JY_geneLength_vs_leadingIsoformLength.png", p, width = 7.2, height = 5.2, dpi = 300)
ggsave("<SELF_NONSELF>/human/RNAseq/JY/JY_geneLength_vs_leadingIsoformLength.pdf", p, width = 7.2, height = 5.2)

print(p)

cat("Spearman correlation (log lengths): ",
    cor(log10(merged$gene_length_kb + 1e-6),
        log10(merged$best_isoform_length_kbCDS + 1e-6),
        method = "spearman", use = "complete.obs"),
    "\n")

#  0.5852967 