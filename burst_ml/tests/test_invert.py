"""Dose inversion against hand-computed values (orchestrator)."""

import numpy as np
from scipy.special import expit, logit

from bml.invert import dose_at_plateau_fraction, dose_for_level


def test_m3a_no_animal_spread_is_exact():
    s = {"beta0": np.array([-1.0]), "beta1": np.array([2.0]), "sigma_a": np.array([0.0])}
    r = dose_for_level("m3a", s, level=0.2, n_a=10)
    want = np.exp((logit(0.2) + 1.0) / 2.0)
    assert np.isclose(r["median"], want) and np.isclose(r["q05"], want)
    assert r["p_unreachable"] == 0
    assert np.isclose(expit(-1.0 + 2.0 * np.log(r["median"])), 0.2)


def test_m3a_negative_slope_is_unreachable():
    s = {"beta0": np.array([-1.0, -1.0]), "beta1": np.array([2.0, -1.0]),
         "sigma_a": np.array([0.0, 0.0])}
    assert np.isclose(dose_for_level("m3a", s, n_a=4)["p_unreachable"], 0.5)


def test_m3b_exact_unreachable_and_floor():
    base = {"p0": np.array([0.01]), "log_ec50": np.array([0.0]), "log_h": np.log([1.5]),
            "sigma_a": np.array([0.0])}
    r = dose_for_level("m3b", {**base, "pmax": np.array([0.6])}, level=0.2, n_a=5)
    H = (0.2 - 0.01) / (0.6 - 0.01)
    assert np.isclose(r["median"], np.exp(logit(H) / 1.5))
    assert dose_for_level("m3b", {**base, "pmax": np.array([0.15])}, level=0.2)["p_unreachable"] == 1
    assert dose_for_level("m3b", {**base, "pmax": np.array([0.6])}, level=0.005)["median"] == 0.0


def test_plateau_fraction():
    s = {"log_ec50": np.array([np.log(1.2)]), "log_h": np.log([2.0]), "sigma_a": np.array([0.0])}
    r = dose_at_plateau_fraction(s, frac=0.9, n_a=3)
    assert np.isclose(r["median"], 1.2 * 9 ** 0.5)
