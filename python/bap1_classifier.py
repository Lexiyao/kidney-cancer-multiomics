"""Module 4 non-circular classifier: predict BAP1 mutation status from expression.

Labels are BAP1 somatic-mutation status on the n=417 mutation subset; features
are log-transformed normalised expression. Because the label is not derived
from the expression features, this task is non-circular (spec §6b, §3.5).

Scope of the held-out figure: the caller passes the top-variable-gene subset,
and that variance filter is computed over ALL samples, including the rows this
module later scores. The filter never sees the label, so there is no leak and
the non-circularity claim stands -- but `heldout_auroc` is transductive in the
feature-selection step and is therefore not fully out-of-sample. The
standardiser, by contrast, lives inside the pipeline and is refitted on the
training rows of every fold.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

RANDOM_STATE = 20160128
N_SPLITS = 5
HELDOUT_FRACTION = 0.30


def build_classifier() -> Pipeline:
    """Standardise features then fit L2 logistic regression (balanced classes)."""
    return Pipeline([
        ("scale", StandardScaler()),
        ("clf", LogisticRegression(penalty="l2", C=1.0, max_iter=1000,
                                   class_weight="balanced",
                                   random_state=RANDOM_STATE)),
    ])


def train_bap1_classifier(expr: pd.DataFrame, labels, random_state: int = RANDOM_STATE) -> dict:
    """Return CV and held-out AUROC for BAP1-from-expression prediction."""
    y = np.asarray(labels)
    if expr.shape[0] != y.shape[0]:
        raise ValueError(
            f"expr has {expr.shape[0]} samples but labels has {y.shape[0]}")
    y = y.astype(int)
    if set(np.unique(y).tolist()) - {0, 1}:
        raise ValueError("labels must be binary 0/1 BAP1 mutation status")
    minority = min(int(y.sum()), int((1 - y).sum()))
    if minority < N_SPLITS:
        raise ValueError(
            f"too few samples in a class ({minority}) for {N_SPLITS}-fold CV")

    x = expr.to_numpy(dtype=float)

    # StandardScaler lives INSIDE the pipeline, so cross_val_predict refits it
    # on each training fold: the scoring rows never contribute to the centring
    # or scaling. Standardising `x` up front would be the classic silent leak.
    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True,
                         random_state=random_state)
    oof = cross_val_predict(build_classifier(), x, y, cv=cv,
                            method="predict_proba")[:, 1]
    cv_auroc = float(roc_auc_score(y, oof))

    x_tr, x_te, y_tr, y_te = train_test_split(
        x, y, test_size=HELDOUT_FRACTION, stratify=y,
        random_state=random_state)
    model = build_classifier().fit(x_tr, y_tr)
    heldout = model.predict_proba(x_te)[:, 1]
    heldout_auroc = float(roc_auc_score(y_te, heldout))

    return {
        "cv_auroc": cv_auroc,
        "heldout_auroc": heldout_auroc,
        "n_samples": len(y),
        "n_bap1_mutant": int(y.sum()),
        "n_features": int(x.shape[1]),
    }
