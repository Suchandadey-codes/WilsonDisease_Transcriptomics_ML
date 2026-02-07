# Probe to gene annotation
library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)

# -----------------------------
# 1. Load ML-derived biomarkers
# -----------------------------
biomarkers <- read.csv(
  "results/final_candidate_biomarkers.csv",
  stringsAsFactors = FALSE
)

# sanity check
stopifnot("gene" %in% colnames(biomarkers))

# -----------------------------
# 2. Annotate probes (ALLOW duplicates)
# -----------------------------
annotation_table <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = unique(biomarkers$gene),
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "PROBEID"
)

# DO NOT set rownames here
head(annotation_table)

# -----------------------------
# 3. Merge annotation with biomarkers
# -----------------------------
biomarkers_annotated <- biomarkers %>%
  left_join(annotation_table, by = c("gene" = "PROBEID"))

# -----------------------------
# 4. Save annotated biomarkers
# -----------------------------
write.csv(
  biomarkers_annotated,
  "results/final_candidate_biomarkers_annotated.csv",
  row.names = FALSE
)

print(biomarkers_annotated)
