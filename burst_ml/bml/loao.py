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
from scipy.special import logsumexp

from .data import subset
from .models import MODELS, Spec, loglik, make_model


def fit(spec: Spec, d: dict, seed: int = 0, chains: int = 4, warmup: int = 1000,
        draws: int = 1000, target_accept: float = 0.95) -> tuple[dict, dict]:
    """NUTS fit of ``make_model(spec)`` to ``d`` (observed ``k = d["k"]``).

    Returns ``(samples, diag)``. ``diag`` = ``{"max_rhat": float, "divergences": int,
    "draws": int}``; ``max_rhat`` over every scalar of every site (split R-hat on
    ``get_samples(group_by_chain=True)``). Chains run in parallel when there are enough
    host devices (``numpyro.set_host_device_count``), else sequentially."""
    import jax
    from numpyro.diagnostics import split_gelman_rubin
    from numpyro.infer import MCMC, NUTS

    method = "parallel" if jax.local_device_count() >= chains else "sequential"
    mcmc = MCMC(NUTS(make_model(spec), target_accept_prob=target_accept), num_warmup=warmup,
                num_samples=draws, num_chains=chains, chain_method=method, progress_bar=False)
    mcmc.run(jax.random.PRNGKey(seed), d, k=d["k"])

    rhats = []
    for v in mcmc.get_samples(group_by_chain=True).values():
        v = np.asarray(v)
        v = v.reshape(v.shape[0], v.shape[1], -1)            # (chains, draws, scalars)
        rhats.extend(np.atleast_1d(np.asarray(split_gelman_rubin(v))).ravel())
    diag = {"max_rhat": float(np.nanmax(rhats)) if rhats else float("nan"),
            "divergences": int(np.asarray(mcmc.get_extra_fields()["diverging"]).sum()),
            "draws": int(chains * draws)}
    return {k: np.asarray(v) for k, v in mcmc.get_samples().items()}, diag


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
    import jax
    import jax.numpy as jnp
    from numpyro.distributions import BetaBinomial

    S_full = len(next(iter(samples.values())))
    idx = np.unique(np.linspace(0, S_full - 1, min(max_draws, S_full)).round().astype(int))
    th = {k: jnp.asarray(np.asarray(v)[idx]) for k, v in samples.items() if k != "z"}
    S = len(idx)
    rng = np.random.default_rng(seed)
    if spec.hier:
        a = rng.standard_normal((S, n_a)) * np.asarray(th["sigma_a"])[:, None]
    else:
        n_a, a = 1, np.zeros((S, 1))
    a = jnp.asarray(a)
    dt = {k: (jnp.asarray(v) if isinstance(v, np.ndarray) else v) for k, v in d_test.items()}
    k_obs = dt["k"]

    def per_draw(f):
        return jax.vmap(lambda th_s, a_s: jax.vmap(lambda aa: f(th_s, aa))(a_s))(th, a)

    ll = np.asarray(per_draw(lambda th_s, aa: loglik(spec, th_s, dt, aa, k_obs)))  # (S, n_a, N)
    # broadcast: M1's p doesn't depend on the targets and comes back as a scalar
    p = per_draw(lambda th_s, aa: jnp.broadcast_to(
        jnp.clip(spec.prob(th_s, dt, aa), 1e-6, 1 - 1e-6), k_obs.shape))
    phi = jnp.exp(th["log_phi"])[:, None, None]
    k_draw = BetaBinomial(p * phi, (1 - p) * phi, total_count=dt["n"]).sample(
        jax.random.PRNGKey(seed))                                             # (S, n_a, N)
    frac = np.asarray(k_draw).reshape(S * n_a, -1) / np.asarray(d_test["n"])

    return {"joint": float(logsumexp(ll.sum(-1)) - np.log(S * n_a)),
            "pointwise": logsumexp(ll, axis=(0, 1)) - np.log(S * n_a),
            "frac": frac}


def loao(spec: Spec, d: dict, seed: int = 0, n_a: int = 200, **fit_kw) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Six folds: hold out each mouse, ``fit`` on the rest, ``predict`` the held-out mouse.

    Returns ``(folds, targets)``:

    * ``folds``: one row per held-out mouse: model, mouse, elpd (the ``joint`` value),
      max_rhat, divergences, draws.
    * ``targets``: one row per held-out target: model, mouse, target, dose, k, n, frac
      (observed k/n), lpd (``pointwise``), pred_median, q05, q25, q75, q95 (of ``frac`` draws),
      pit (share of ``frac`` draws below the observed fraction, plus half the ties).
    """
    mice = sorted(set(d["mouse"]))
    folds_rows = []
    targets_rows = []

    for mouse in mice:
        # Hold out this mouse
        holdout_mask = d["mouse"] == mouse
        d_train = subset(d, ~holdout_mask)
        d_holdout = subset(d, holdout_mask)

        # Fit on training data
        samples, diag_out = fit(spec, d_train, seed=seed, **fit_kw)

        # Predict held-out mouse
        pred = predict(spec, samples, d_holdout, n_a=n_a, seed=seed + int(mouse))

        folds_rows.append({
            "model": spec.name,
            "mouse": int(mouse),
            "elpd": pred["joint"],
            "max_rhat": diag_out["max_rhat"],
            "divergences": diag_out["divergences"],
            "draws": diag_out["draws"],
        })

        # Build target rows
        frac_draws = pred["frac"]  # (S * n_a, N_holdout)
        n_holdout = len(d_holdout["k"])

        for i in range(n_holdout):
            frac_i = frac_draws[:, i]
            frac_obs = d_holdout["k"][i] / d_holdout["n"][i]

            # pit: share of frac draws below observed fraction, plus half the ties
            below = (frac_i < frac_obs).sum()
            ties = (frac_i == frac_obs).sum()
            pit = (below + 0.5 * ties) / len(frac_i)

            targets_rows.append({
                "model": spec.name,
                "mouse": int(mouse),
                "target": int(d_holdout["target"][i]),
                "dose": float(d_holdout["dose"][i]),
                "k": int(d_holdout["k"][i]),
                "n": int(d_holdout["n"][i]),
                "frac": frac_obs,
                "lpd": pred["pointwise"][i],
                "pred_median": float(np.median(frac_i)),
                "q05": float(np.percentile(frac_i, 5)),
                "q25": float(np.percentile(frac_i, 25)),
                "q75": float(np.percentile(frac_i, 75)),
                "q95": float(np.percentile(frac_i, 95)),
                "pit": float(pit),
            })

    folds = pd.DataFrame(folds_rows)
    targets = pd.DataFrame(targets_rows)
    return folds, targets


def compare(folds: pd.DataFrame, base: str = "m3a") -> pd.DataFrame:
    """Columns model, elpd, delta, se, z, beats. Per model: elpd (sum over folds), delta vs ``base`` (sum of per-mouse differences),
    se = sqrt(n_folds)·std(per-mouse differences, ddof=1), z = delta/se (NaN for the base),
    beats = delta > 0 and z > 2 (the README's decision rule). Sorted by elpd, descending."""
    n_folds = len(folds["mouse"].unique())

    results = []
    for model in folds["model"].unique():
        model_folds = folds[folds["model"] == model].set_index("mouse")
        elpd = float(model_folds["elpd"].sum())

        if model == base:
            results.append({
                "model": model,
                "elpd": elpd,
                "delta": 0.0,
                "se": 0.0,
                "z": np.nan,
                "beats": False,
            })
        else:
            base_folds = folds[folds["model"] == base].set_index("mouse")
            diffs = model_folds["elpd"] - base_folds["elpd"]
            delta = float(diffs.sum())
            se = float(np.sqrt(n_folds) * diffs.std(ddof=1))
            z = delta / se if se > 0 else np.nan
            beats = bool(delta > 0 and z > 2)

            results.append({
                "model": model,
                "elpd": elpd,
                "delta": delta,
                "se": se,
                "z": z,
                "beats": beats,
            })

    out = pd.DataFrame(results).sort_values("elpd", ascending=False)
    return out.reset_index(drop=True)


def calibration(targets: pd.DataFrame) -> pd.DataFrame:
    """One row per model with columns model, cover50, cover90, mean_pit, mean_lpd, spearman:
    the share of observed fractions inside the 50% (q25-q75) and 90% (q05-q95) intervals
    (inclusive), the mean pit, the mean lpd, and the mean over mice of the within-mouse
    Spearman ρ between pred_median and frac."""
    from scipy.stats import spearmanr

    results = []
    for model in targets["model"].unique():
        mt = targets[targets["model"] == model]

        # cover50: share where frac is inside [q25, q75]
        cover50 = float(((mt["frac"] >= mt["q25"]) & (mt["frac"] <= mt["q75"])).mean())

        # cover90: share where frac is inside [q05, q95]
        cover90 = float(((mt["frac"] >= mt["q05"]) & (mt["frac"] <= mt["q95"])).mean())

        # mean_pit
        mean_pit = float(mt["pit"].mean())

        # mean_lpd
        mean_lpd = float(mt["lpd"].mean())

        # spearman: mean over mice of within-mouse Spearman rho between pred_median and frac
        rho_list = []
        for mouse in mt["mouse"].unique():
            mouse_data = mt[mt["mouse"] == mouse]
            if len(mouse_data) >= 2:
                rho, _ = spearmanr(mouse_data["pred_median"], mouse_data["frac"])
                rho_list.append(rho)
        spearman = float(np.mean(rho_list)) if rho_list else np.nan

        results.append({
            "model": model,
            "cover50": cover50,
            "cover90": cover90,
            "mean_pit": mean_pit,
            "mean_lpd": mean_lpd,
            "spearman": spearman,
        })

    return pd.DataFrame(results)