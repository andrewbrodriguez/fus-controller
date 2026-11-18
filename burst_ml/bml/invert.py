"""Dose needed for a delivery level in a new animal (README "Dose for a delivery level").

Both functions take full-data ``samples`` (dict of arrays with leading draw axis S) and draw
``n_a`` fresh animal effects per draw, ``a_new ~ Normal(0, sigma_a[s])``, giving S·n_a values
of the required dose. Report the median, 5% and 95% quantiles over the *reachable* values,
and the share that is unreachable.

M3a: ``logit p = beta0 + a + beta1·log D``  →  ``log D = (logit(level) − beta0 − a) / beta1``.
     Unreachable when ``beta1 <= 0``.
M3b: ``p = p0 + (pmax − p0)·expit(h·(log D − log_ec50 − a))``, ``h = exp(log_h)``.
     With ``H = (level − p0)/(pmax − p0)``: ``log D = log_ec50 + a + logit(H)/h``.
     Unreachable when ``pmax <= level`` (H >= 1). When ``level <= p0`` the dose is 0.
     Dose at 90% of the plateau: ``H = 0.9`` → ``log D = log_ec50 + a + log(9)/h``.
"""

from __future__ import annotations

import numpy as np
from scipy.special import logit


def _a_new(samples: dict, n_a: int, seed: int) -> np.ndarray:
    sigma = np.asarray(samples["sigma_a"], float)
    return np.random.default_rng(seed).standard_normal((len(sigma), n_a)) * sigma[:, None]


def _summary(dose: np.ndarray, reachable: np.ndarray) -> dict:
    """Median and 5-95% of the reachable doses, and the unreachable share."""
    vals = dose[reachable]
    q = np.quantile(vals, [0.5, 0.05, 0.95]) if vals.size else [np.nan] * 3
    return {"median": float(q[0]), "q05": float(q[1]), "q95": float(q[2]),
            "p_unreachable": float(1 - reachable.mean()), "n": int(dose.size)}


def dose_for_level(model: str, samples: dict, level: float = 0.2, n_a: int = 200,
                   seed: int = 0) -> dict:
    """``{"model", "level", "median", "q05", "q95", "p_unreachable", "n"}`` for model
    ``"m3a"`` or ``"m3b"`` (module doc)."""
    a = _a_new(samples, n_a, seed)                                    # (S, n_a)
    if model == "m3a":
        b0 = np.asarray(samples["beta0"], float)[:, None]
        b1 = np.asarray(samples["beta1"], float)[:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            dose = np.exp((logit(level) - b0 - a) / b1)
        reachable = np.broadcast_to(b1 > 0, a.shape)
    elif model == "m3b":
        p0 = np.asarray(samples["p0"], float)[:, None]
        pmax = np.asarray(samples["pmax"], float)[:, None]
        h = np.exp(np.asarray(samples["log_h"], float))[:, None]
        H = np.broadcast_to((level - p0) / (pmax - p0), a.shape)
        reachable = H < 1                                             # pmax > level
        with np.errstate(divide="ignore", invalid="ignore"):
            dose = np.exp(np.asarray(samples["log_ec50"], float)[:, None] + a
                          + logit(np.clip(H, 1e-12, 1 - 1e-12)) / h)
        dose = np.where(H <= 0, 0.0, dose)                            # level at or below p0
    else:
        raise ValueError(f"unknown model {model!r}")
    return {"model": model, "level": level, **_summary(dose, reachable)}


def dose_at_plateau_fraction(samples: dict, frac: float = 0.9, n_a: int = 200,
                             seed: int = 0) -> dict:
    """M3b only: dose reaching ``frac`` of the way from p0 to pmax in a new animal.
    ``{"model": "m3b", "plateau_frac", "median", "q05", "q95", "n"}``."""
    a = _a_new(samples, n_a, seed)
    h = np.exp(np.asarray(samples["log_h"], float))[:, None]
    dose = np.exp(np.asarray(samples["log_ec50"], float)[:, None] + a + logit(frac) / h)
    out = _summary(dose, np.ones_like(dose, bool))
    out.pop("p_unreachable")
    return {"model": "m3b", "plateau_frac": frac, **out}
