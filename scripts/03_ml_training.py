############################################################
# 03_ml_training.py
# Dataset: GSE197406 (Wilson's Disease)
# Purpose: ML training & cross-validation (TRAINING SET ONLY)
# Author: Suchanda Dey
############################################################

import numpy as np
import pandas as pd

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

# -----------------------------
# 1. Load data
# -----------------------------
expr = pd.read_csv(
    "/mnt/d/BiotecNika/New_analyses/R/data/processed/expr_matrix_normalized.csv",
    index_col=0
)

labels = pd.read_csv(
    "/mnt/d/BiotecNika/New_analyses/R/data/processed/sample_labels.csv",
    index_col=0
)["label"]

train_samples = pd.read_csv(
    "/mnt/d/BiotecNika/New_analyses/R/data/processed/train_samples.csv"
).iloc[:, 0].values

# -----------------------------
# 2. Load DEGs (training-only)
# -----------------------------
deg = pd.read_csv(
    "/mnt/d/BiotecNika/New_analyses/R/results/DEGs_training_only_filtered.csv",
    index_col=0
)

print(f"Total DEGs from limma (filtered): {deg.shape[0]}")

# ---- choose number of genes (important) ----
TOP_N_GENES = min(200, deg.shape[0])
top_genes = deg.index[:TOP_N_GENES]

print(f"Using top {len(top_genes)} genes for ML")

# -----------------------------
# 3. Build TRAINING feature matrix
# -----------------------------
X_train = expr.loc[top_genes, train_samples].T
y_train = labels.loc[train_samples].values

print("Training feature matrix shape:", X_train.shape)
print("Training labels distribution:", np.bincount(y_train))

# -----------------------------
# 4. Define ML pipeline
# -----------------------------
pipeline = Pipeline([
    ("scaler", StandardScaler()),
    ("clf", LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        l1_ratio=0.5,        # balance between L1 & L2
        max_iter=5000,
        random_state=RANDOM_STATE
    ))
])

# -----------------------------
# 5. Cross-validation (TRAIN ONLY)
# -----------------------------
cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=RANDOM_STATE
)

cv_scores = cross_val_score(
    pipeline,
    X_train,
    y_train,
    cv=cv,
    scoring="roc_auc"
)

print("\nCross-validation ROC-AUC scores:")
print(cv_scores)
print("Mean CV ROC-AUC:", cv_scores.mean())
print("Std CV ROC-AUC:", cv_scores.std())

# -----------------------------
# 6. Fit model on FULL training set
# -----------------------------
pipeline.fit(X_train, y_train)

# -----------------------------
# 7. Extract selected genes (ElasticNet)
# -----------------------------
model = pipeline.named_steps["clf"]
coef = pd.Series(
    model.coef_[0],
    index=top_genes
)

selected_genes = coef[coef != 0].sort_values(
    key=np.abs,
    ascending=False
)

print(f"\nSelected genes (non-zero coefficients): {selected_genes.shape[0]}")

# -----------------------------
# 8. Save outputs
# -----------------------------
selected_genes.to_csv(
    "/mnt/d/BiotecNika/New_analyses/R/results/elasticnet_selected_genes_training.csv",
    header=["coefficient"]
)

pd.DataFrame({
    "cv_auc": cv_scores
}).to_csv(
    "/mnt/d/BiotecNika/New_analyses/R/results/cv_auc_scores_training.csv",
    index=False
)

print("\nML training completed successfully!")
print("Saved:")
print("- /mnt/d/BiotecNika/New_analyses/R/results/elasticnet_selected_genes_training.csv")
print("- /mnt/d/BiotecNika/New_analyses/R/results/cv_auc_scores_training.csv")
