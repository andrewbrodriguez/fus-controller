"""Fitting and leave-one-animal-out scoring (burst_ml/README.md "Fitting and scoring").

Conventions:

* ``samples``: dict ``name -> array`` with a leading draw axis S (chains flattened), as from
  ``MCMC.get_samples()``. It holds every site of the model, including ``log_phi`` and, for
  hierarchical models, ``sigma_a`` and ``z``.
* A draw ``th`` is ``{name: samples[name][s]}`` for one s; vectorise over s with ``jax.vmap``.
* Held-out animals get fresh effects: for each draw s, ``n_a`` values
  ``a_new ~ Normal(0, sigma_a[s])``. Non-hierarchical models use ``a_new = 0`` (n_a = 1).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .data import subset
from .models import MODELS, Spec, loglik, make_model


def fit(spec: Spec, d: dict, seed: int = 0, chains: int = 4, warmup: int = 1000,
        draws: int = 1000, target_accept: float = 0.95) -> tuple[dict, dict]:
    """NUTS fit of ``make_model(spec)`` to ``d`` (observed ``k = d["k"]``).

    Returns ``(samples, diag)``. ``diag`` = ``{"max_rhat": float, "divergences": int,
    "draws": int}``; ``max_rhat`` over every scalar of every site, using
    ``numpyro.diagnostics.split_gelman_rubin`` on ``get_samples(group_by_chain=True)``.
    Use ``chain_method="parallel"`` when ``jax.local_device_count() >= chains``, else
    ``"sequential"``; ``progress_bar=False``."""
    raise NotImplementedError


def predict(spec: Spec, samples: dict, d_test: dict, n_a: int = 200, seed: int = 0,
            max_draws: int = 1000) -> dict:
    """Held-out prediction for one animal (or any set of rows sharing one new animal).

    Thin ``samples`` evenly to at most ``max_draws`` draws (S). For each draw and each of
    ``n_a`` animal effects, compute ``ll = loglik(spec, th, d_test, a_new, d_test["k"])``,
    shape (S, n_a, N_test). Return:

    * ``"joint"``: ``logsumexp(ll.sum(-1)) − log(S·n_a)`` -- the animal's targets jointly.
    * ``"pointwise"``: (N_test,) ``logsumexp(ll, axis=(0, 1)) − log(S·n_a)``.
    * ``"frac"``: (S·n_a, N_test) predicted GFP+ fractions, one BetaBinomial draw of k per
      (draw, a_new) divided by ``d_test["n"]``.
    """
    raise NotImplementedError


def loao(spec: Spec, d: dict, seed: int = 0, **fit_kw) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Six folds: hold out each mouse, ``fit`` on the rest, ``predict`` the held-out mouse.

    Returns ``(folds, targets)``:

    * ``folds``: one row per held-out mouse: model, mouse, elpd (the ``joint`` value),
      max_rhat, divergences, draws.
    * ``targets``: one row per held-out target: model, mouse, target, dose, k, n, frac
      (observed k/n), lpd (``pointwise``), pred_median, q05, q25, q75, q95 (of ``frac`` draws),
      pit (share of ``frac`` draws below the observed fraction, plus half the ties).
    """
    raise NotImplementedError


def compare(folds: pd.DataFrame, base: str = "m3a") -> pd.DataFrame:
    """Columns model, elpd, delta, se, z, beats. Per model: elpd (sum over folds), delta vs ``base`` (sum of per-mouse differences),
    se = sqrt(n_folds)·std(per-mouse differences, ddof=1), z = delta/se (NaN for the base),
    beats = delta > 0 and z > 2 (the README's decision rule). Sorted by elpd, descending."""
    raise NotImplementedError


def calibration(targets: pd.DataFrame) -> pd.DataFrame:
    """One row per model with columns model, cover50, cover90, mean_pit, mean_lpd, spearman:
    the share of observed fractions inside the 50% (q25-q75) and 90% (q05-q95) intervals
    (inclusive), the mean pit, the mean lpd, and the mean over mice of the within-mouse
    Spearman ρ between pred_median and frac."""
    raise NotImplementedError
