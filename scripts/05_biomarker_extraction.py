import pandas as pd
import numpy as np
import os

# -----------------------------
# 1. Define paths
# -----------------------------
RESULTS_DIR = "/mnt/d/BiotecNika/New_analyses/R/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

INPUT_FILE = os.path.join(
    RESULTS_DIR,
    "elasticnet_selected_genes_training.csv"
)

OUTPUT_FILE = os.path.join(
    RESULTS_DIR,
    "final_candidate_biomarkers.csv"
)

# -----------------------------
# 2. Load ElasticNet-selected genes
# -----------------------------
biomarkers = pd.read_csv(
    INPUT_FILE,
    index_col=0
)

# rename column safely
biomarkers.columns = ["coefficient"]

print("Total ML-selected genes:", biomarkers.shape[0])
print(biomarkers.head())

# -----------------------------
# 3. Rank by importance
# -----------------------------
biomarkers["abs_coefficient"] = biomarkers["coefficient"].abs()

biomarkers_ranked = biomarkers.sort_values(
    by="abs_coefficient",
    ascending=False
)

# -----------------------------
# 4. Select top biomarkers
# -----------------------------
TOP_K = 15
top_biomarkers = biomarkers_ranked.head(TOP_K)

# make gene names explicit (VERY IMPORTANT)
top_biomarkers = top_biomarkers.reset_index()
top_biomarkers.columns = ["gene", "coefficient", "abs_coefficient"]

print("\nTop candidate biomarkers:")
print(top_biomarkers)

# -----------------------------
# 5. Save final biomarker table
# -----------------------------
top_biomarkers.to_csv(
    OUTPUT_FILE,
    index=False
)

print("\nSaved final biomarkers to:")
print(OUTPUT_FILE)
