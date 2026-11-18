"""Aggregate figures from ``burst_ml/results/*.csv`` and ``cache/samples_*.npz``.

From ``burst_ml/``: ``.venv/bin/python -m bml.figures [--results DIR] [--samples DIR] [--out DIR]``
(defaults ``results/``, ``cache/``, ``results/``; for development use
``--results cache/quick --samples cache/quick --out cache/quick``). Writes PNGs (dpi 150).
CSV paths below are relative to ``--results`` and ``samples_*.npz`` to ``--samples``. Every figure pools all six mice. Use matplotlib's Agg backend,
colour mice with ``tab10`` in mouse order (1-6), label axes with units, and put the source
CSV's name in each figure's lower-right corner in small grey text.

1. ``comparison.png`` (from ``comparison.csv`` rows with base m3a, and ``loao_folds.csv``):
   one row per model, x = Δelpd vs m3a, a dot at delta with a ±2·se bar, and the six
   per-mouse differences as small grey dots. Vertical line at 0.
2. ``loao_predictions.png`` (from ``loao_targets.csv``): two panels, m3a and the model with
   the highest elpd other than m3a. x = observed GFP+ fraction, y = pred_median with q05-q95
   bars, coloured by mouse, y = x line, both axes 0-1.
3. ``dose_response.png`` (from ``cache/samples_m3a.npz``, ``cache/samples_m3b.npz`` and the
   labels via ``bml.data``): x = dose (log axis, 0.2-4), y = GFP+ fraction. For each of m3a and
   m3b, the new-animal predictive median and 50% / 90% bands of p (draw a_new per posterior
   draw as in ``bml.invert``); observed targets as points coloured by mouse; horizontal dotted
   line at 0.2.
4. ``burst_params.png``: posterior histograms of ``gamma`` (m4), ``exp(log_rho)`` (m5a) and the
   three ``beta_bands`` (m5b) from the full-data samples, each with its prior density as a
   line and a vertical line at the value that reduces the model to m3a (0, 1, 0).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

_MOUSE_ORDER = list(range(1, 7))


def _source_label(csv_name: str) -> str:
    """Return the label to place in the lower-right corner."""
    return f"Source: {csv_name}"


def _load_csv(results: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(results / name)


def _load_npz(samples: Path, name: str) -> dict:
    return dict(np.load(str(samples / name)))


# ---------------------------------------------------------------------------
# 1. comparison.png
# ---------------------------------------------------------------------------

def _fig_comparison(results: Path, out: Path) -> None:
    csv_name = "comparison.csv, loao_folds.csv"
    comp = _load_csv(results, "comparison.csv")
    folds = _load_csv(results, "loao_folds.csv")

    # Rows with base == "m3a" (first block)
    base = comp[comp["base"] == "m3a"].copy()

    # Per-mouse deltas from the same base block
    # The per-mouse rows have the same model/mouse combos but from the second block
    # where base == the model being compared against. We need the per-mouse deltas
    # from the block where base == the model being compared.
    # Actually, looking at the data: rows 7-14 have base="m1", rows 0-6 have base="m3a".
    # The per-mouse dots should come from the comparison block where base matches the
    # reference model. For the m3a-based comparison (rows 0-6), we need per-mouse deltas.
    # Those are in the rows where base == the model name, e.g. for comparing m4 vs m3a,
    # the per-mouse deltas are in rows where model=="m4" and base=="m4" (second block).
    # But that's complicated. Let's look at the structure more carefully.
    #
    # Actually the second block has base == model name for each row.
    # Row 7: model=m4, base=m1 -> this is m4 vs m1 (per-mouse)
    # Row 8: model=m3a, base=m1 -> m3a vs m1 (per-mouse)
    # etc.
    #
    # For the "base m3a" comparison, the per-mouse differences are:
    # In the second block, find rows where model == m3a's model name (none, m3a is the ref)
    # or where base == model (per-mouse) and the model appears in the first block.
    #
    # Simpler: the per-mouse dots are from rows where base == model (i.e. beats is True
    # or the row is a per-mouse comparison). Let's use rows where base == model name
    # and the model is in the first block.
    w = folds.pivot(index="mouse", columns="model", values="elpd")
    order = [m for m in ["m1", "m2", "m3a", "m3b", "m4", "m5a", "m5b"] if m in set(base["model"])]
    base = base.set_index("model").loc[order].reset_index()

    fig, ax = plt.subplots(figsize=(6, 5))

    palette = dict(zip(_MOUSE_ORDER, plt.get_cmap("tab10").colors))

    # Six per-mouse differences (held-out mouse elpd minus m3a's) as small grey dots
    for m in order:
        diffs = (w[m] - w["m3a"]).to_numpy()
        ax.plot(diffs, [m] * len(diffs), "o", color="grey", markersize=3, alpha=0.6, zorder=2)

    # Plot aggregate dot + error bar for each model in first block
    for _, row in base.iterrows():
        model = row["model"]
        delta = row["delta"]
        se = row["se"]
        if np.isnan(se):
            ax.plot(delta, model, "o", color="C0", markersize=6, zorder=3)
        else:
            ax.errorbar(
                delta, model, xerr=2 * se, fmt="o", color="C0",
                markersize=6, capsize=3, zorder=3,
            )

    ax.axvline(0, color="black", linewidth=0.8, linestyle="--", zorder=1)
    ax.set_xlabel("Δ elpd vs m3a (log predictive density, held-out mice)")
    ax.set_ylabel("Model")
    ax.set_title("Model comparison (base: m3a)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(
        0.99, 0.03, _source_label(csv_name),
        transform=ax.transAxes, fontsize=7, ha="right", va="bottom",
        color="grey",
    )
    fig.tight_layout()
    fig.savefig(out / "comparison.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 2. loao_predictions.png
# ---------------------------------------------------------------------------

def _fig_loao_predictions(results: Path, out: Path) -> None:
    csv_name = "loao_targets.csv"
    targets = _load_csv(results, csv_name)
    folds = _load_csv(results, "loao_folds.csv")

    # Find the best model (highest elpd) other than m3a
    # elpd is per-mouse in folds; we need the overall best.
    # Use the overall elpd from the folds table.
    # Actually, looking at loao_folds.csv, each row is one model-mouse combo.
    # We want the model with the highest (least negative) elpd overall.
    # But elpd is per-mouse. Let's sum elpd by model.
    elpd_by_model = folds.groupby("model")["elpd"].sum()
    others = elpd_by_model.drop("m3a", errors="ignore")
    if len(others) == 0:
        best_model = "m3a"
    else:
        best_model = others.idxmax()

    # Also make a combined figure with two panels
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    for ax, model in zip(axes, ["m3a", best_model]):
        mdata = targets[targets["model"] == model].copy()
        palette = dict(zip(_MOUSE_ORDER, plt.get_cmap("tab10").colors))

        for mouse in _MOUSE_ORDER:
            msub = mdata[mdata["mouse"] == mouse]
            color = palette.get(mouse, "C0")
            ax.errorbar(
                msub["frac"], msub["pred_median"],
                yerr=[
                    msub["pred_median"] - msub["q05"],
                    msub["q95"] - msub["pred_median"],
                ],
                fmt="o", color=color, markersize=4, capsize=2,
                label=f"Mouse {mouse}", alpha=0.8,
            )

        ax.plot([0, 1], [0, 1], "k--", linewidth=0.8, alpha=0.5)
        ax.set_xlabel("Observed GFP+ fraction")
        ax.set_ylabel("Predicted median GFP+ fraction")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title(f"{model}")
        ax.legend(fontsize=7, loc="upper left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(
            0.99, 0.03, _source_label(csv_name),
            transform=ax.transAxes, fontsize=7, ha="right", va="bottom",
            color="grey",
        )

    fig.suptitle("LOAO Predictions", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / "loao_predictions.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 3. dose_response.png
# ---------------------------------------------------------------------------

def _fig_dose_response(samples: Path, out: Path) -> None:
    csv_name = "samples_m3a.npz, samples_m3b.npz"

    # Load samples
    s3a = _load_npz(samples, "samples_m3a.npz")
    s3b = _load_npz(samples, "samples_m3b.npz")

    from bml.data import load_design, load_labels

    design = load_design()
    labels = load_labels()

    fig, ax = plt.subplots(figsize=(7, 5))

    palette = dict(zip(_MOUSE_ORDER, plt.get_cmap("tab10").colors))

    # --- m3a predictive curves ---
    # For each posterior draw: draw a_new ~ Normal(0, sigma_a), then
    # p = expit(beta0 + a_new + beta1 * log(dose))
    # Compute median, 50%, 90% bands over draws for each dose value.
    doses = np.geomspace(0.2, 4.0, 200)   # sorted grid; the observed doses are drawn as points
    beta0 = s3a["beta0"]  # (S,)
    beta1 = s3a["beta1"]  # (S,)
    sigma_a = s3a["sigma_a"]  # (S,)

    rng = np.random.default_rng(42)
    n_draws = len(beta0)
    n_targets = len(doses)

    # Draw a_new for each draw and each target dose
    a_new = rng.standard_normal((n_draws, n_targets)) * sigma_a[:, None]  # (S, N)

    from jax.scipy.special import expit as _expit
    import jax.numpy as jnp

    dose_arr = jnp.asarray(doses)
    log_dose = jnp.log(dose_arr)

    # Compute p for each draw: (S, N)
    p_m3a = _expit(beta0[:, None] + a_new + beta1[:, None] * log_dose[None, :])

    # Also compute population-level (a=0) for reference
    p_m3a_pop = _expit(beta0[:, None] + beta1[:, None] * log_dose[None, :])

    # Median, 50%, 90% bands
    p_med = np.quantile(p_m3a, 0.5, axis=0)
    p_50lo = np.quantile(p_m3a, 0.25, axis=0)
    p_50hi = np.quantile(p_m3a, 0.75, axis=0)
    p_90lo = np.quantile(p_m3a, 0.05, axis=0)
    p_90hi = np.quantile(p_m3a, 0.95, axis=0)

    ax.fill_between(
        doses, p_90lo, p_90hi, alpha=0.15, color="C0", label="m3a 90% band (new animal)"
    )
    ax.fill_between(
        doses, p_50lo, p_50hi, alpha=0.3, color="C0", label="m3a 50% band (new animal)"
    )
    ax.plot(doses, p_med, "C0-", linewidth=1.5, label="m3a median (new animal)")

    # --- m3b predictive curves ---
    log_ec50 = s3b["log_ec50"]
    log_h = s3b["log_h"]
    p0 = s3b["p0"]
    pmax = s3b["pmax"]
    sigma_a_3b = s3b["sigma_a"]

    a_new_3b = rng.standard_normal((n_draws, n_targets)) * sigma_a_3b[:, None]

    h = np.exp(log_h)[:, None]  # (S, 1)
    log_ec50_arr = log_ec50[:, None]  # (S, 1)

    # p = p0 + (pmax - p0) * expit(h * (log_dose - log_ec50 - a_new))
    p_m3b = p0[:, None] + (pmax[:, None] - p0[:, None]) * _expit(
        h * (log_dose[None, :] - log_ec50_arr - a_new_3b)
    )

    p_med_3b = np.quantile(p_m3b, 0.5, axis=0)
    p_50lo_3b = np.quantile(p_m3b, 0.25, axis=0)
    p_50hi_3b = np.quantile(p_m3b, 0.75, axis=0)
    p_90lo_3b = np.quantile(p_m3b, 0.05, axis=0)
    p_90hi_3b = np.quantile(p_m3b, 0.95, axis=0)

    ax.fill_between(
        doses, p_90lo_3b, p_90hi_3b, alpha=0.15, color="C1",
        label="m3b 90% band (new animal)",
    )
    ax.fill_between(
        doses, p_50lo_3b, p_50hi_3b, alpha=0.3, color="C1",
        label="m3b 50% band (new animal)",
    )
    ax.plot(
        doses, p_med_3b, "C1-", linewidth=1.5, label="m3b median (new animal)"
    )

    # --- observed targets as points, coloured by mouse ---
    labels = labels.copy()
    labels["frac"] = labels["k"] / labels["n"]
    for _, row in labels.iterrows():
        mouse = int(row["mouse"])
        color = palette.get(mouse, "C0")
        ax.plot(
            row["dose"], row["frac"],
            "o", color=color, markersize=4, alpha=0.7, zorder=5,
        )

    ax.axhline(0.2, color="grey", linewidth=0.8, linestyle=":", zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(0.2, 4)
    ax.set_xlabel("Dose: cumulative 2nd harmonic (lab units, log scale)")
    ax.set_ylabel("GFP+ fraction of NeuN neurons")
    ax.set_title("Dose–response: model predictions vs observed")
    ax.legend(fontsize=7, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(
        0.99, 0.03, _source_label(csv_name),
        transform=ax.transAxes, fontsize=7, ha="right", va="bottom",
        color="grey",
    )
    fig.tight_layout()
    fig.savefig(out / "dose_response.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 4. burst_params.png
# ---------------------------------------------------------------------------

def _fig_burst_params(samples: Path, out: Path) -> None:
    csv_name = "samples_m4/m5a/m5b.npz"

    s4 = _load_npz(samples, "samples_m4.npz")
    s5a = _load_npz(samples, "samples_m5a.npz")
    s5b = _load_npz(samples, "samples_m5b.npz")

    from bml.models import MODELS

    fig, axes = plt.subplots(1, 5, figsize=(20, 4))

    palette = dict(zip(_MOUSE_ORDER, plt.get_cmap("tab10").colors))

    # --- 1. gamma (m4) ---
    gamma = s4["gamma"]
    ax = axes[0]
    prior_gamma = MODELS["m4"].priors["gamma"]
    # Plot histogram
    ax.hist(gamma, bins=30, density=True, alpha=0.7, color="C0", edgecolor="white")
    # Plot prior density
    from scipy.stats import norm
    x = np.linspace(gamma.min() - 0.5, gamma.max() + 0.5, 200)
    ax.plot(x, norm.pdf(x, 0, 1), "C0-", linewidth=1.5, label="Prior")
    ax.axvline(0, color="black", linewidth=1.2, linestyle="--", label="Reduces to m3a")
    ax.set_xlabel("γ, time weighting per 100 bursts (0 = flat)")
    ax.set_ylabel("Density")
    ax.set_title("m4: γ")
    ax.legend(fontsize=7)

    # --- 2. exp(log_rho) (m5a) ---
    log_rho = s5a["log_rho"]
    rho = np.exp(log_rho)
    ax = axes[1]
    prior_rho = MODELS["m5a"].priors["log_rho"]
    ax.hist(rho, bins=30, density=True, alpha=0.7, color="C1", edgecolor="white")
    # Prior on log_rho is Normal(0, 0.5), so rho is log-normal
    x = np.linspace(0, rho.max(), 200)
    from scipy.stats import lognorm
    # lognorm with s=0.5, loc=0, scale=exp(0)=1
    ax.plot(x, lognorm.pdf(x, s=0.5, scale=np.exp(0)), "C1-", linewidth=1.5, label="Prior")
    ax.axvline(1.0, color="black", linewidth=1.2, linestyle="--", label="Reduces to m3a")
    ax.set_xlabel("ρ, per-burst power (1 = plain dose)")
    ax.set_ylabel("Density")
    ax.set_title("m5a: ρ")
    ax.legend(fontsize=7)

    # --- 3-5. beta_bands (m5b) ---
    beta_bands = s5b["beta_bands"]  # (S, 3)
    band_labels = ["β broadband", "β subharmonic", "β ultraharmonic"]
    for i in range(3):
        ax = axes[2 + i]
        bb = beta_bands[:, i]
        prior_bb = MODELS["m5b"].priors["beta_bands"]
        ax.hist(bb, bins=30, density=True, alpha=0.7, color="C2", edgecolor="white")
        x = np.linspace(bb.min() - 0.2, bb.max() + 0.2, 200)
        ax.plot(x, norm.pdf(x, 0, 0.25), "C2-", linewidth=1.5, label="Prior")
        ax.axvline(0, color="black", linewidth=1.2, linestyle="--", label="Reduces to m3a")
        ax.set_xlabel(band_labels[i])
        ax.set_ylabel("Density")
        ax.set_title(f"m5b: {band_labels[i]}")
        ax.legend(fontsize=7)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(
            0.99, 0.03, _source_label(csv_name),
            transform=ax.transAxes, fontsize=7, ha="right", va="bottom",
            color="grey",
        )

    fig.suptitle("Burst-model parameters: posteriors vs priors", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(out / "burst_params.png", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Generate aggregate figures")
    parser.add_argument(
        "--results", type=Path, default=None,
        help="Directory with CSV outputs (default: results/)",
    )
    parser.add_argument(
        "--samples", type=Path, default=None,
        help="Directory with samples_*.npz (default: cache/)",
    )
    parser.add_argument(
        "--out", type=Path, default=None,
        help="Output directory (default: results/)",
    )
    args = parser.parse_args(argv)

    from bml import RESULTS, CACHE

    results = args.results or RESULTS
    samples = args.samples or CACHE
    out = args.out or RESULTS

    out.mkdir(parents=True, exist_ok=True)

    _fig_comparison(results, out)
    _fig_loao_predictions(results, out)
    _fig_dose_response(samples, out)
    _fig_burst_params(samples, out)

    print(f"Wrote 4 figures to {out}")


if __name__ == "__main__":
    main()