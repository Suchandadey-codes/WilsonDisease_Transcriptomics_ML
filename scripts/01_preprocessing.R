############################################################
# 01_preprocessing.R
# Dataset: GSE197406 (Wilson's Disease)
# Purpose: Preprocessing & normalization ONLY
# Author: Suchanda Dey
############################################################

# -----------------------------
# 1. Load required libraries
# -----------------------------
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(dplyr)
})

set.seed(42)

# -----------------------------
# 2. Create output directories
# -----------------------------
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 3. Download GEO dataset
# -----------------------------
message("Downloading GSE197406 from GEO...")

gse <- getGEO(
  "GSE197406",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

gse <- gse[[1]]

# -----------------------------
# 4. Extract expression & metadata
# -----------------------------
expr_raw <- exprs(gse)
pheno <- pData(gse)

message("Expression matrix dimensions (raw):")
print(dim(expr_raw))

# Save raw data (optional but good practice)
write.csv(expr_raw, "data/raw/expr_matrix_raw.csv")
write.csv(pheno, "data/raw/pheno_data_raw.csv")

# -----------------------------
# 5. Log2 transform (if needed)
# -----------------------------
# Check if data looks already log-transformed
qx <- quantile(expr_raw, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)

if (qx[5] > 100) {
  message("Applying log2(x + 1) transformation...")
  expr_raw <- log2(expr_raw + 1)
} else {
  message("Data already appears log-transformed.")
}

# -----------------------------
# 6. Quantile normalization
# -----------------------------
message("Performing quantile normalization...")
expr_norm <- normalizeBetweenArrays(expr_raw, method = "quantile")

# -----------------------------
# 7. Remove control probes (AFFX)
# -----------------------------
message("Removing AFFX/control probes...")
expr_norm <- expr_norm[!grepl("^AFFX", rownames(expr_norm)), ]

# -----------------------------
# 8. Remove low-variance genes
# -----------------------------
message("Filtering low-variance genes...")

gene_variance <- apply(expr_norm, 1, var)
variance_cutoff <- quantile(gene_variance, 0.25)

expr_norm <- expr_norm[gene_variance > variance_cutoff, ]

message("Expression matrix dimensions (filtered):")
print(dim(expr_norm))

# -----------------------------
# 9. Extract disease labels
# -----------------------------
# IMPORTANT:
# Check pheno column names carefully if this fails
colnames(pheno)

# ---- MODIFY THIS LINE IF NEEDED ----
# Common names: disease_status, condition, group, phenotype

group_column <- "disease state:ch1"

# Inspect class distribution
print(table(pheno[[group_column]]))

labels <- ifelse(
  pheno[[group_column]] == "Human WD liver (Cirrhosis)",
  1, 0
)

names(labels) <- rownames(pheno)

# Final sanity check
print(table(labels))

# -----------------------------
# 10. Match samples (CRITICAL)
# -----------------------------
common_samples <- intersect(colnames(expr_norm), names(labels))

expr_norm <- expr_norm[, common_samples]
labels <- labels[common_samples]

stopifnot(
  all(colnames(expr_norm) == names(labels))
)

# -----------------------------
# 11. Save processed outputs
# -----------------------------
write.csv(
  expr_norm,
  "data/processed/expr_matrix_normalized.csv"
)

write.csv(
  data.frame(label = labels),
  "data/processed/sample_labels.csv"
)

# -----------------------------
# 12. Summary
# -----------------------------
message("Preprocessing completed successfully!")
message("Final samples: ", ncol(expr_norm))
message("Final genes: ", nrow(expr_norm))
message("Files saved in data/processed/")
