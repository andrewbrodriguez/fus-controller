"""Model formulas checked against independent numpy/scipy (written by the orchestrator)."""

import numpy as np
from scipy.special import expit
from scipy.stats import betabinom

from bml.models import MODELS, dose_m4, dose_m5a, dose_m5b, loglik


def _design(rng, N=4, T=6):
    h2 = rng.uniform(20, 200, (N, T))
    mask = np.ones((N, T), bool)
    mask[0, 4:] = False
    h2[~mask] = 0.0
    Z = rng.normal(size=(N, T, 3))
    Z[~mask] = 0.0
    t = np.tile(np.arange(T, dtype=float), (N, 1))
    dose = (h2 * mask).sum(1) / 1e4
    return {"animal": np.array([0, 0, 1, 1]), "n_animals": 2, "k": np.array([10, 200, 50, 900]),
            "n": np.array([1000, 1000, 800, 2000]), "dose": dose, "h2": h2, "t": t,
            "mask": mask, "Z": Z, "h_ref": 77.5}


def test_loglik_matches_scipy():
    d = _design(np.random.default_rng(0))
    th = {"beta0": -1.0, "beta1": 1.2, "log_phi": np.log(12.0)}
    a = np.array([0.3, 0.3, -0.2, -0.2])
    p = expit(-1.0 + a + 1.2 * np.log(d["dose"]))
    ref = betabinom.logpmf(d["k"], d["n"], p * 12, (1 - p) * 12)
    got = np.asarray(loglik(MODELS["m3a"], th, d, a, d["k"]))
    assert np.allclose(got, ref, atol=1e-8)


def test_m3b_reduces_to_m3a():
    """p0 = 0, pmax = 1, h = beta1, log_ec50 = -beta0/beta1: the two curves coincide
    (with the animal effect's sign and scale flipped: a_m3b = -a_m3a / beta1)."""
    d = _design(np.random.default_rng(1))
    b0, b1 = -1.1, 1.3
    a = np.array([0.4, 0.4, -0.5, -0.5])
    p3a = MODELS["m3a"].prob({"beta0": b0, "beta1": b1}, d, a)
    th3b = {"p0": 0.0, "pmax": 1.0, "log_ec50": -b0 / b1, "log_h": np.log(b1)}
    p3b = MODELS["m3b"].prob(th3b, d, -a / b1)
    assert np.allclose(np.asarray(p3a), np.asarray(p3b), atol=1e-10)


def test_burst_doses_by_hand():
    d = _design(np.random.default_rng(2))
    m = d["mask"]
    g = 0.7
    w = np.exp(g * (d["t"] / 100 - 0.5))
    assert np.allclose(np.asarray(dose_m4({"gamma": g}, d)), (w * d["h2"] * m).sum(1) / 1e4)
    rho = 1.4
    x = np.where(m, d["h2"] / 77.5, 1.0)
    ref5a = 77.5 / 1e4 * np.where(m, x ** rho, 0).sum(1)
    assert np.allclose(np.asarray(dose_m5a({"log_rho": np.log(rho)}, d)), ref5a)
    bb = np.array([0.1, -0.2, 0.05])
    ref5b = 77.5 / 1e4 * np.where(m, x ** rho * np.exp(d["Z"] @ bb), 0).sum(1)
    got5b = dose_m5b({"log_rho": np.log(rho), "beta_bands": bb}, d)
    assert np.allclose(np.asarray(got5b), ref5b)


def test_identity_burst_models_equal_dose():
    d = _design(np.random.default_rng(3))
    assert np.allclose(np.asarray(dose_m4({"gamma": 0.0}, d)), d["dose"])
    assert np.allclose(np.asarray(dose_m5a({"log_rho": 0.0}, d)), d["dose"])
    assert np.allclose(np.asarray(dose_m5b({"log_rho": 0.0, "beta_bands": np.zeros(3)}, d)),
                       d["dose"])
