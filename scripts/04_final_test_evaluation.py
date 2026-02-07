############################################################
# 04_final_test_evaluation.py
# Dataset: GSE197406 (Wilson's Disease)
# Purpose: FINAL evaluation on HELD-OUT TEST SET
# Author: Suchanda Dey
############################################################

import numpy as np
import pandas as pd

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score

print(">>> 04_final_test_evaluation.py started <<<")


RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

# -----------------------------
# 1. Load processed data
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

test_samples = pd.read_csv(
    "/mnt/d/BiotecNika/New_analyses/R/data/processed/test_samples.csv"
).iloc[:, 0].values

# -----------------------------
# 2. Load DEGs (training-only)
# -----------------------------
deg = pd.read_csv(
    "/mnt/d/BiotecNika/New_analyses/R/results/DEGs_training_only_filtered.csv",
    index_col=0
)

TOP_N_GENES = min(200, deg.shape[0])
top_genes = deg.index[:TOP_N_GENES]

# -----------------------------
# 3. Build feature matrices
# -----------------------------
X_train = expr.loc[top_genes, train_samples].T
y_train = labels.loc[train_samples].values

X_test = expr.loc[top_genes, test_samples].T
y_test = labels.loc[test_samples].values

print("Train shape:", X_train.shape)
print("Test shape:", X_test.shape)
print("Test label distribution:", np.bincount(y_test))

# -----------------------------
# 4. Define model (SAME AS TRAINING)
# -----------------------------
pipeline = Pipeline([
    ("scaler", StandardScaler()),
    ("clf", LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        l1_ratio=0.5,
        max_iter=5000,
        random_state=RANDOM_STATE
    ))
])

# -----------------------------
# 5. Train on FULL training set
# -----------------------------
pipeline.fit(X_train, y_train)

# -----------------------------
# 6. Final test evaluation
# -----------------------------
y_prob = pipeline.predict_proba(X_test)[:, 1]
y_pred = pipeline.predict(X_test)

test_auc = roc_auc_score(y_test, y_prob)
test_acc = accuracy_score(y_test, y_pred)

print("\nFINAL TEST RESULTS")
print("------------------")
print("Test ROC-AUC:", test_auc)
print("Test Accuracy:", test_acc)

# -----------------------------
# 7. Save predictions
# -----------------------------
import os

RESULTS_DIR = "/mnt/d/BiotecNika/New_analyses/R/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

results_df = pd.DataFrame({
    "sample": test_samples,
    "true_label": y_test,
    "predicted_prob": y_prob,
    "predicted_label": y_pred
})

output_path = os.path.join(
    RESULTS_DIR,
    "test_set_predictions.csv"
)

results_df.to_csv(output_path, index=False)

print("\nSaved:")
print(output_path)

print(">>> Script finished <<<")
