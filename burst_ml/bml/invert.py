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


def dose_for_level(model: str, samples: dict, level: float = 0.2, n_a: int = 200,
                   seed: int = 0) -> dict:
    """``{"model", "level", "median", "q05", "q95", "p_unreachable", "n"}`` for model
    ``"m3a"`` or ``"m3b"`` (module doc)."""
    raise NotImplementedError


def dose_at_plateau_fraction(samples: dict, frac: float = 0.9, n_a: int = 200,
                             seed: int = 0) -> dict:
    """M3b only: dose reaching ``frac`` of the way from p0 to pmax in a new animal.
    ``{"model": "m3b", "plateau_frac", "median", "q05", "q95", "n"}``."""
    raise NotImplementedError
