############################################################
# 02_limma_train_only.R
# Dataset: GSE197406 (Wilson's Disease)
# Purpose: Differential expression on TRAINING SET ONLY
# Author: Suchanda Dey
############################################################

# -----------------------------
# 1. Load libraries
# -----------------------------
suppressPackageStartupMessages({
  library(limma)
})

set.seed(42)

# -----------------------------
# 2. Load processed data
# -----------------------------
expr <- read.csv(
  "data/processed/expr_matrix_normalized.csv",
  row.names = 1,
  check.names = FALSE
)

labels_df <- read.csv(
  "data/processed/sample_labels.csv",
  row.names = 1
)

labels <- labels_df$label
names(labels) <- rownames(labels_df)

# sanity check
stopifnot(all(colnames(expr) == names(labels)))

message("Total samples: ", length(labels))
message("Class distribution:")
print(table(labels))

# -----------------------------
# 3. Train / Test split (base R, stratified)
# -----------------------------
message("Creating stratified train-test split (base R)...")

idx_control <- which(labels == 0)
idx_wilson  <- which(labels == 1)

train_control <- sample(
  idx_control,
  size = floor(0.8 * length(idx_control))
)

train_wilson <- sample(
  idx_wilson,
  size = floor(0.8 * length(idx_wilson))
)

train_idx <- c(train_control, train_wilson)
test_idx  <- setdiff(seq_along(labels), train_idx)

train_samples <- colnames(expr)[train_idx]
test_samples  <- colnames(expr)[test_idx]

message("Training set distribution:")
print(table(labels[train_samples]))

message("Test set distribution:")
print(table(labels[test_samples]))

# save sample IDs
write.csv(
  train_samples,
  "data/processed/train_samples.csv",
  row.names = FALSE
)

write.csv(
  test_samples,
  "data/processed/test_samples.csv",
  row.names = FALSE
)

# -----------------------------
# 4. Prepare training matrix
# -----------------------------
expr_train <- expr[, train_samples]
labels_train <- labels[train_samples]

group <- factor(
  labels_train,
  levels = c(0, 1),
  labels = c("Control", "Wilson")
)

design <- model.matrix(~ group)

message("Design matrix:")
print(design)

# -----------------------------
# 5. limma differential expression
# -----------------------------
message("Running limma on TRAINING SET ONLY...")

fit <- lmFit(expr_train, design)
fit <- eBayes(fit)

deg_all <- topTable(
  fit,
  coef = "groupWilson",
  adjust.method = "BH",
  number = Inf,
  sort.by = "P"
)

# -----------------------------
# 6. DEG filtering (FEATURE FILTERING ONLY)
# -----------------------------
deg_filtered <- deg_all[
  deg_all$adj.P.Val < 0.05 &
    abs(deg_all$logFC) > 1,
]

message("DEGs after filtering: ", nrow(deg_filtered))

# -----------------------------
# 7. Save DEG results
# -----------------------------
dir.create("results", showWarnings = FALSE)

write.csv(
  deg_all,
  "results/DEGs_training_only_all.csv"
)

write.csv(
  deg_filtered,
  "results/DEGs_training_only_filtered.csv"
)

# -----------------------------
# 8. Summary
# -----------------------------
message("limma analysis completed successfully!")
message("IMPORTANT:")
message("- limma was run ONLY on training samples")
message("- DEGs are for FEATURE FILTERING, not final biomarkers")
message("- Test set has NOT been touched")

##################################################
# 9. visualization
library(ggplot2)
library(dplyr)

# -----------------------------
# 9.1. Load limma results (training only)
# -----------------------------
deg <- read.csv(
  "results/DEGs_training_only_all.csv",
  stringsAsFactors = FALSE
)

# -----------------------------
# 9.2. Prepare data for volcano plot
# -----------------------------
deg <- deg %>%
  mutate(
    negLogP = -log10(adj.P.Val),
    regulation = case_when(
      adj.P.Val < 0.05 & logFC > 1  ~ "Up in Wilson disease",
      adj.P.Val < 0.05 & logFC < -1 ~ "Up in Control",
      TRUE                          ~ "Not significant"
    )
  )

table(deg$regulation)

# -----------------------------
# 9.3. Volcano plot
# -----------------------------
volcano <- ggplot(deg, aes(x = logFC, y = negLogP)) +
  geom_point(aes(color = regulation), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c(
      "Up in Wilson disease" = "red",
      "Up in Control" = "blue",
      "Not significant" = "grey70"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano plot of differentially expressed genes (training set)",
    x = "log2 Fold Change (Wilson vs Control)",
    y = "-log10 adjusted P-value"
  ) +
  theme_minimal()

print(volcano)

# -----------------------------
# 9.4. Save plot
# -----------------------------
ggsave(
  "results/volcano_plot_training_only.png",
  volcano,
  width = 7,
  height = 6,
  dpi = 300
)
#------------------------------------
# Heatmap---------------------------
#--------------------------------------
#Step 1: Load required libraries
library(pheatmap)
library(dplyr)
#Step 2: Load expression data and biomarkers
# expression matrix
expr <- read.csv(
  "data/processed/expr_matrix_normalized.csv",
  row.names = 1,
  check.names = FALSE
)

# final annotated biomarkers
bio <- read.csv(
  "results/final_candidate_biomarkers_annotated.csv",
  stringsAsFactors = FALSE
)

# use gene SYMBOLS (remove NA)
bio_genes <- bio %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(SYMBOL)

genes_to_plot <- bio_genes$SYMBOL
#Step 3: Map SYMBOL ??? expression rows
# add SYMBOL annotation to expression
expr$PROBEID <- rownames(expr)

annotation_map <- bio %>%
  select(gene, SYMBOL)

expr_annot <- expr %>%
  inner_join(annotation_map, by = c("PROBEID" = "gene"))

# collapse probes ??? gene-level
expr_gene <- expr_annot %>%
  select(-PROBEID) %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), mean))

expr_gene <- as.data.frame(expr_gene)
rownames(expr_gene) <- expr_gene$SYMBOL
expr_gene$SYMBOL <- NULL
#Step 4: Prepare sample annotation
labels <- read.csv(
  "data/processed/sample_labels.csv",
  row.names = 1
)

annotation_col <- data.frame(
  Condition = factor(
    labels$label,
    levels = c(0, 1),
    labels = c("Control", "Wilson disease")
  )
)

rownames(annotation_col) <- rownames(labels)
#Step 5: Scale expression
expr_scaled <- t(scale(t(expr_gene)))
#Step 6: Generate heatmap
pheatmap(
  expr_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Heatmap of ML-derived candidate biomarkers"
)
#Step 7: Save heatmap
pheatmap(
  expr_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Heatmap of ML-derived candidate biomarkers",
  filename = "results/biomarker_heatmap.png",
  width = 8,
  height = 6
)

