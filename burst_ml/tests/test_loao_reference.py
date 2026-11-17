"""Scoring checked against independent scipy computations (written by the orchestrator)."""

import numpy as np
import pandas as pd
import pytest
from scipy.special import expit, logsumexp
from scipy.stats import betabinom

from bml.loao import calibration, compare, predict
from bml.models import MODELS


def _design(dose, k, n):
    N = len(dose)
    return {"mouse": np.full(N, 9), "target": np.arange(1, N + 1), "animal": np.zeros(N, int),
            "n_animals": 1, "k": np.asarray(k), "n": np.asarray(n),
            "dose": np.asarray(dose, float), "h2": np.ones((N, 2)), "t": np.zeros((N, 2)),
            "mask": np.ones((N, 2), bool), "Z": np.zeros((N, 2, 3)), "h_ref": 1.0}


def _ll(p, phi, k, n):
    return betabinom.logpmf(k, n, p * phi, (1 - p) * phi)


def test_predict_pooled_matches_scipy():
    d = _design([0.5, 1.0, 2.0], [30, 200, 900], [1000, 1000, 2000])
    s = {"beta0": np.array([-1.0, -1.5, -0.5]), "beta1": np.array([1.0, 1.4, 0.8]),
         "log_phi": np.log(np.array([8.0, 12.0, 20.0]))}
    out = predict(MODELS["m2"], s, d, n_a=5, seed=0)
    ll = np.array([[_ll(expit(s["beta0"][i] + s["beta1"][i] * np.log(d["dose"][j])),
                        np.exp(s["log_phi"][i]), d["k"][j], d["n"][j]) for j in range(3)]
                   for i in range(3)])                                    # (S, N)
    assert np.isclose(out["joint"], logsumexp(ll.sum(1)) - np.log(3), atol=1e-6)
    assert np.allclose(out["pointwise"], logsumexp(ll, axis=0) - np.log(3), atol=1e-6)
    assert out["frac"].shape[1] == 3
    assert ((out["frac"] >= 0) & (out["frac"] <= 1)).all()


def test_predict_hier_with_zero_sigma_matches_fixed_effect():
    """sigma_a = 0 makes every a_new 0, so m3a must equal the pooled formula."""
    d = _design([0.5, 2.0], [50, 700], [1000, 1500])
    s = {"beta0": np.array([-1.2, -0.8]), "beta1": np.array([1.1, 1.3]),
         "log_phi": np.log(np.array([10.0, 15.0])), "sigma_a": np.array([0.0, 0.0]),
         "z": np.zeros((2, 6))}
    out = predict(MODELS["m3a"], s, d, n_a=7, seed=1)
    ll = np.array([[_ll(expit(s["beta0"][i] + s["beta1"][i] * np.log(d["dose"][j])),
                        np.exp(s["log_phi"][i]), d["k"][j], d["n"][j]) for j in range(2)]
                   for i in range(2)])
    assert np.isclose(out["joint"], logsumexp(ll.sum(1)) - np.log(2), atol=1e-6)


def test_predict_hier_spreads_with_sigma():
    """A large sigma_a must widen the predictive distribution of the fraction."""
    d = _design([1.0], [200], [1000])
    base = {"beta0": np.full(50, -1.4), "beta1": np.full(50, 1.0),
            "log_phi": np.full(50, np.log(50.0)), "z": np.zeros((50, 6))}
    narrow = predict(MODELS["m3a"], {**base, "sigma_a": np.full(50, 0.01)}, d, n_a=40, seed=2)
    wide = predict(MODELS["m3a"], {**base, "sigma_a": np.full(50, 2.0)}, d, n_a=40, seed=2)
    assert wide["frac"].std() > 2 * narrow["frac"].std()


def test_compare_decision_rule():
    folds = pd.DataFrame({
        "model": ["m3a"] * 6 + ["m4"] * 6 + ["m5a"] * 6,
        "mouse": list(range(1, 7)) * 3,
        "elpd": [-10.0] * 6 + [-9.0, -9.1, -8.9, -9.0, -9.2, -8.8] + [-12, -8, -11, -9, -10, -10.0],
    })
    c = compare(folds, base="m3a").set_index("model")
    diff4 = np.array([1.0, 0.9, 1.1, 1.0, 0.8, 1.2])
    assert np.isclose(c.loc["m4", "delta"], diff4.sum())
    assert np.isclose(c.loc["m4", "se"], np.sqrt(6) * diff4.std(ddof=1))
    assert bool(c.loc["m4", "beats"]) is True
    assert bool(c.loc["m5a", "beats"]) is False          # delta 0, noisy
    assert np.isnan(c.loc["m3a", "z"])
    assert c["elpd"].is_monotonic_decreasing


def test_calibration_coverage():
    t = pd.DataFrame({"model": ["m"] * 4, "mouse": [1, 1, 2, 2], "target": [1, 2, 1, 2],
                      "frac": [0.1, 0.5, 0.2, 0.9], "pred_median": [0.1, 0.4, 0.3, 0.8],
                      "q05": [0.0, 0.3, 0.25, 0.1], "q25": [0.05, 0.45, 0.28, 0.6],
                      "q75": [0.2, 0.6, 0.35, 0.85], "q95": [0.3, 0.7, 0.5, 0.95],
                      "pit": [0.5, 0.6, 0.1, 0.9], "lpd": [-1.0, -2.0, -3.0, -4.0]})
    c = calibration(t).set_index("model")
    assert np.isclose(c.loc["m", "cover50"], 0.5)        # 0.1 and 0.5 inside; 0.2 and 0.9 not
    assert np.isclose(c.loc["m", "cover90"], 0.75)       # 0.2 is below q05 = 0.25
    assert np.isclose(c.loc["m", "mean_lpd"], -2.5)
