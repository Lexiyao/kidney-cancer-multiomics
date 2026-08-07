"""Module 4: the non-circular BAP1-from-expression classifier.

The label (BAP1 somatic-mutation status) is EXTERNAL to the features
(log-transformed normalised expression), so this task is not the tautology
spec section 6b forbids. The tests below check the reported numbers, the input
guards, and — most importantly — that a run on pure noise does NOT come back
with a good AUROC, which is what a label leak would look like.
"""
import numpy as np
import pytest

from python.bap1_classifier import build_classifier, train_bap1_classifier


def test_build_classifier_returns_pipeline_with_scaler_and_lr():
    pipe = build_classifier()
    steps = dict(pipe.named_steps)
    assert "scale" in steps and "clf" in steps


def test_train_bap1_classifier_reports_cv_and_heldout_auroc(bap1_fixture):
    expr, labels = bap1_fixture
    out = train_bap1_classifier(expr, labels)
    assert set(out) == {"cv_auroc", "heldout_auroc", "n_samples",
                        "n_bap1_mutant", "n_features"}
    assert 0.0 <= out["cv_auroc"] <= 1.0
    assert 0.0 <= out["heldout_auroc"] <= 1.0
    # signal is planted in the fixture -> better than chance
    assert out["cv_auroc"] > 0.7
    assert out["n_samples"] == expr.shape[0]
    assert out["n_bap1_mutant"] == int(np.asarray(labels).sum())


def test_train_bap1_classifier_rejects_shape_mismatch(bap1_fixture):
    expr, labels = bap1_fixture
    with pytest.raises(ValueError, match="samples"):
        train_bap1_classifier(expr.iloc[:-3], labels)


def test_train_bap1_classifier_rejects_non_binary_labels(bap1_fixture):
    expr, labels = bap1_fixture
    bad = np.asarray(labels).copy()
    bad[0] = 2
    with pytest.raises(ValueError, match="binary"):
        train_bap1_classifier(expr, bad)


# --- added guards (not in the plan) -----------------------------------------


def test_cv_auroc_is_near_chance_when_features_carry_no_signal(bap1_null_fixture):
    """The leak detector.

    Features are pure noise and the labels are unrelated to them, so the only
    honest answer is an out-of-fold AUROC near 0.5. If the standardiser were
    fitted outside the folds, or the labels reached the scoring rows in any
    other way, this figure would climb towards 1.0. A generous ceiling is used
    so the test flags leakage, not ordinary sampling noise.
    """
    expr, labels = bap1_null_fixture
    out = train_bap1_classifier(expr, labels)
    assert out["cv_auroc"] < 0.75


def test_reported_counts_match_the_input_matrix(bap1_fixture):
    expr, labels = bap1_fixture
    out = train_bap1_classifier(expr, labels)
    assert out["n_features"] == expr.shape[1]
    assert out["n_samples"] == len(labels)


def test_same_seed_reproduces_identical_metrics(bap1_fixture):
    expr, labels = bap1_fixture
    first = train_bap1_classifier(expr, labels, random_state=7)
    second = train_bap1_classifier(expr, labels, random_state=7)
    assert first == second


def test_rejects_a_minority_class_too_small_for_cv(bap1_fixture):
    """The `minority < N_SPLITS` guard: 5-fold stratified CV cannot be built
    from 3 mutants, and sklearn would only warn rather than stop."""
    expr, labels = bap1_fixture
    scarce = np.asarray(labels).copy()
    scarce[:] = 0
    scarce[:3] = 1
    with pytest.raises(ValueError, match="too few"):
        train_bap1_classifier(expr, scarce)


def test_accepts_labels_as_a_plain_python_list(bap1_fixture):
    """The reticulate calling convention (Task 4.9).

    `reticulate::r_to_py()` turns an R integer vector into a Python *list*, not
    a numpy array, so `train_bap1_classifier` has to accept one. This is the
    only part of the `bap1_auroc` wiring that can be executed outside the
    container, so it is asserted here rather than assumed.
    """
    expr, labels = bap1_fixture
    as_list = [int(v) for v in labels]
    out = train_bap1_classifier(expr, as_list)
    assert out["n_bap1_mutant"] == sum(as_list)
    assert out["n_samples"] == len(as_list)
    assert out == train_bap1_classifier(expr, labels)


def test_does_not_mutate_its_inputs(bap1_fixture):
    expr, labels = bap1_fixture
    expr_before = expr.copy(deep=True)
    labels_before = np.asarray(labels).copy()
    train_bap1_classifier(expr, labels)
    assert expr.equals(expr_before)
    assert np.array_equal(np.asarray(labels), labels_before)
