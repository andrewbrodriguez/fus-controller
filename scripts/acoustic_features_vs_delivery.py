"""Does the full acoustic record predict delivery better than the 2nd-harmonic AUC alone?

    .venv/bin/python scripts/acoustic_features_vs_delivery.py

Features: the fixed list in ``fus.features.FEATURES`` (raw waveforms + controller log),
computed for every target of mice 1-6 and cached in ``target_features.csv``. Only mice
1-2 have delivery, so the comparison runs on their **11 sonicated targets** (Mouse 1 T3,
the no-FUS control, has no recording and so no features).

Delivery: cumulative over slices, relative to the animal's mean target, from
``results/histology/slices/cumulative_dose_delivery.csv`` (run
``scripts/slices_dose_delivery.py`` first). Pipeline A (GFP+ area) is primary; B
(GFP+ cells) is reported alongside.

With 11 points, any multi-feature model fits better in-sample, so every model is scored
on **leave-one-target-out** predictions: Q² = 1 - PRESS / SS. Models:

* null -- the mean of the other targets
* 2f AUC -- least squares on ``cum_2nd_harmonic`` (the current analysis)
* ridge, all features -- standardised, penalty chosen by an inner leave-one-out
* PLS, all features -- one component
* ridge / PLS on log features -- the same, with every positive continuous feature
  log-transformed. **Added after the first run**, in which ridge on raw features
  failed (Q² -20) on one left-out target with an extreme ``h2_per_f0``; energies and
  ratios span orders of magnitude, so log is the usual scale. Reported as post hoc.
* 2f AUC + one other feature -- exploratory; 13 of them, so expect one to look good

A permutation test (shuffle delivery, rerun everything) asks whether ridge's Q² gain
over 2f AUC exceeds what 11 random labels give. Leave-one-mouse-out (fit one animal,
predict the other) is reported as a harder check.

Writes to ``results/acoustic_features/``: target_features.csv, model_comparison.csv,
feature_screen.csv, loo_predictions.csv, model_comparison.png, README.md.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from fus.features import FEATURES, target_features

ROOT = Path(__file__).resolve().parents[1]
ACOUSTIC = ROOT / "data/acoustic/20260611"
DELIVERY = ROOT / "results/histology/slices/cumulative_dose_delivery.csv"
OUT = ROOT / "results/acoustic_features"
CACHE = OUT / "target_features.csv"
N_PERM = 2000
ALPHAS = np.logspace(-3, 3, 25)

# Recordings per target; Mouse 1 was the operator deviation (data/README.md).
RECORDINGS = {1: {1: ["Target1", "Target2", "Target3"], 2: ["Target2_Repeat"],
                  4: ["Target4"], 5: ["Target5"], 6: ["Target6"]}}
for _m in range(2, 7):
    RECORDINGS[_m] = {t: [f"Target{t}"] for t in range(1, 7)}

OUTCOMES = {"gfp_area_mm2_rel": ("A", "GFP+ area"), "gfp_cells_rel": ("B", "GFP+ cells")}
FEATS = list(FEATURES)

SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
MOUSE_COLOR = {1: "#2a78d6", 2: "#eb6834"}
MODEL_COLOR = "#2a78d6"


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


def features() -> pd.DataFrame:
    if CACHE.exists():
        return pd.read_csv(CACHE)
    rows = []
    for mouse, targets in RECORDINGS.items():
        for target, files in targets.items():
            paths = [ACOUSTIC / f"Mouse_Cntr_{mouse:02d}_{f}.mat" for f in files]
            if not all(p.exists() for p in paths):
                print(f"  mouse {mouse} target {target}: recording missing, skipped")
                continue
            print(f"  mouse {mouse} target {target} ({', '.join(files)})")
            rows.append({"mouse": mouse, "target": target, **target_features(paths)})
    d = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    d.to_csv(CACHE, index=False)
    return d


def labelled(f: pd.DataFrame) -> pd.DataFrame:
    deliv = pd.read_csv(DELIVERY)[["mouse", "target", *OUTCOMES]]
    return f.merge(deliv, on=["mouse", "target"]).sort_values(["mouse", "target"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------


def _null():
    class Mean:
        def fit(self, X, y):
            self.m = float(np.mean(y))
            return self

        def predict(self, X):
            return np.full(len(X), self.m)
    return Mean()


def _ridge():
    return make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS))


def _pls():
    return make_pipeline(StandardScaler(), PLSRegression(n_components=1))


#: Kept on their own scale under the log variant: a 0/1 flag and a bounded setpoint.
NOT_LOGGED = {"interlock", "harmonic_goal"}

# name -> (model factory, feature columns, log-transform?)
MODELS = {
    "null (mean)": (_null, ["cum_2nd_harmonic"], False),
    "2f AUC": (LinearRegression, ["cum_2nd_harmonic"], False),
    "ridge, all features": (_ridge, FEATS, False),
    "PLS, all features": (_pls, FEATS, False),
    "ridge, all features (log)": (_ridge, FEATS, True),
    "PLS, all features (log)": (_pls, FEATS, True),
}
RIDGES = ["ridge, all features", "ridge, all features (log)"]


def design(d: pd.DataFrame, cols: list[str], log: bool) -> np.ndarray:
    X = d[cols].astype(float).copy()
    if log:
        for c in cols:
            if c not in NOT_LOGGED:
                X[c] = np.log10(X[c])
    return X.to_numpy()


def loo(make, X: np.ndarray, y: np.ndarray, groups: np.ndarray | None = None) -> np.ndarray:
    """Out-of-sample predictions: leave one target out, or one group (mouse) out."""
    keys = np.arange(len(y)) if groups is None else groups
    pred = np.empty(len(y))
    for k in np.unique(keys):
        test = keys == k
        model = make().fit(X[~test], y[~test])
        pred[test] = np.ravel(model.predict(X[test]))
    return pred


def q2(y: np.ndarray, pred: np.ndarray) -> float:
    return float(1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2))


def compare(d: pd.DataFrame, outcome: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    y = d[outcome].to_numpy()
    rows, preds = [], {}
    for name, (make, cols, log) in MODELS.items():
        X = design(d, cols, log)
        p = loo(make, X, y)
        preds[name] = p
        lomo = loo(make, X, y, groups=d.mouse.to_numpy())
        rows.append({"outcome": outcome, "model": name, "n_features": len(cols), "q2_loo": q2(y, p),
                     "rmse_loo": float(np.sqrt(np.mean((y - p) ** 2))),
                     # the null's LOO predictions anti-correlate with y by construction
                     "rho_loo": stats.spearmanr(y, p).statistic if name != "null (mean)" else np.nan,
                     "q2_leave_mouse_out": q2(y, lomo)})
    res = pd.DataFrame(rows)

    # Permutation: shuffle delivery, rerun. p = how often random labels do at least as well.
    rng = np.random.default_rng(0)
    X2 = design(d, ["cum_2nd_harmonic"], False)
    Xr = {m: design(d, MODELS[m][1], MODELS[m][2]) for m in RIDGES}
    q = res.set_index("model").q2_loo
    null = {m: [] for m in ["2f AUC", *RIDGES]}
    for _ in range(N_PERM):
        yp = rng.permutation(y)
        null["2f AUC"].append(q2(yp, loo(LinearRegression, X2, yp)))
        for m in RIDGES:
            null[m].append(q2(yp, loo(_ridge, Xr[m], yp)))
    null = {k: np.array(v) for k, v in null.items()}
    res["perm_p_q2"] = res.model.map({m: float(np.mean(v >= q[m])) for m, v in null.items()})
    res["perm_p_gain_over_2f"] = res.model.map(
        {m: float(np.mean(null[m] - null["2f AUC"] >= q[m] - q["2f AUC"])) for m in RIDGES})

    p = pd.DataFrame({"mouse": d.mouse, "target": d.target, "observed": y,
                      **{f"pred_{k}": v for k, v in preds.items()}})
    p.insert(0, "outcome", outcome)
    return res, p


def screen(d: pd.DataFrame, f: pd.DataFrame) -> pd.DataFrame:
    """Per feature: redundancy with 2f AUC (36 targets) and LOO gain when added to it (11)."""
    rows = []
    for col in FEATS:
        row = {"feature": col, "what": FEATURES[col],
               "rho_with_2f_auc_36": stats.spearmanr(f[col], f.cum_2nd_harmonic).statistic
               if f[col].nunique() > 1 else np.nan}
        for outcome, (pipe, _) in OUTCOMES.items():
            y = d[outcome].to_numpy()
            row[f"rho_delivery_{pipe}"] = (stats.spearmanr(d[col], y).statistic
                                           if d[col].nunique() > 1 else np.nan)
            if col != "cum_2nd_harmonic" and d[col].nunique() > 1:
                X = d[["cum_2nd_harmonic", col]].to_numpy(float)
                base = q2(y, loo(LinearRegression, d[["cum_2nd_harmonic"]].to_numpy(float), y))
                row[f"q2_gain_{pipe}"] = q2(y, loo(lambda: LinearRegression(), X, y)) - base
            else:
                row[f"q2_gain_{pipe}"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------


def _clean(ax, grid_axis="both"):
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis=grid_axis, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(length=0)


def figure(res: pd.DataFrame, preds: pd.DataFrame):
    plt.rcParams.update({"font.family": ["Helvetica Neue", "Arial", "DejaVu Sans"], "font.size": 9.5,
                         "text.color": INK, "axes.labelcolor": INK_2, "xtick.color": MUTED,
                         "ytick.color": MUTED, "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE})
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.4), gridspec_kw={"width_ratios": [1.15, 1, 1]})
    fig.subplots_adjust(left=0.17, right=0.98, top=0.85, bottom=0.08, hspace=0.45, wspace=0.28)
    for r, (outcome, (pipe, name)) in enumerate(OUTCOMES.items()):
        rs = res[res.outcome == outcome].set_index("model")
        ax = axes[r, 0]
        _clean(ax, "x")
        names = list(MODELS)[::-1]
        vals = [rs.loc[n, "q2_loo"] for n in names]
        floor = -1.0  # bars below this are clipped; the label keeps the true value
        ax.barh(names, [max(v, floor) for v in vals], color=MODEL_COLOR, height=0.55, zorder=2)
        ax.axvline(0, color=AXIS, linewidth=1)
        for i, v in enumerate(vals):
            x = max(v, floor)
            ax.text(x + (0.03 if v >= 0 else -0.03), i, f"{v:+.2f}" + (" (clipped)" if v < floor else ""),
                    va="center", ha="left" if v >= 0 else "right", fontsize=8.5, color=INK_2)
        ax.set_xlim(floor - 0.75, 1.0)
        ax.set_title(f"Pipeline {pipe} ({name}): leave-one-out Q²", loc="left", fontsize=10,
                     weight="bold", color=INK)
        ax.tick_params(axis="y", labelcolor=INK)
        pg = rs.loc["ridge, all features (log)", "perm_p_gain_over_2f"]
        ax.text(1.0, -0.2, f"ridge (log) gain over 2f AUC: permutation p = {pg:.2f}", transform=ax.transAxes,
                ha="right", fontsize=8, color=MUTED)
        p = preds[preds.outcome == outcome]
        for c, model in enumerate(["2f AUC", "ridge, all features (log)"], start=1):
            ax = axes[r, c]
            _clean(ax)
            hi = max(p.observed.max(), p[f"pred_{model}"].max()) * 1.1
            lo = min(0, p[f"pred_{model}"].min()) - 0.05
            ax.plot([lo, hi], [lo, hi], color=AXIS, linestyle=":", linewidth=1.2, zorder=1)
            for mouse in (1, 2):
                m = p[p.mouse == mouse]
                ax.scatter(m[f"pred_{model}"], m.observed, s=170, color=MOUSE_COLOR[mouse],
                           edgecolors=SURFACE, linewidths=2, zorder=3,
                           label=f"Mouse {mouse}" if (r, c) == (0, 1) else None)
                for t, x, yv in zip(m.target, m[f"pred_{model}"], m.observed):
                    ax.annotate(str(t), (x, yv), ha="center", va="center_baseline", fontsize=8,
                                weight="bold", color="#ffffff", zorder=4)
            ax.set_xlim(lo, hi)
            ax.set_ylim(lo, hi)
            ax.set_title(f"{model}: Q² {rs.loc[model, 'q2_loo']:+.2f}", loc="left", fontsize=10,
                         weight="bold", color=INK)
            ax.set_xlabel("predicted (left out)", color=INK_2)
            if c == 1:
                ax.set_ylabel("observed, relative to animal's mean", color=INK_2)
            if (r, c) == (0, 1):
                ax.legend(loc="upper left", frameon=False, fontsize=8.5)
    fig.text(0.012, 0.985, "Full acoustic features vs 2nd-harmonic AUC alone, predicting delivery",
             ha="left", va="top", fontsize=14, weight="bold", color=INK)
    fig.text(0.012, 0.95, "11 sonicated targets, mice 1–2. Each prediction is made with that target "
             "left out. Q² = 1 − PRESS/SS: 0 is no better than the mean, negative is worse.",
             ha="left", va="top", fontsize=9.5, color=INK_2)
    fig.text(0.012, 0.008, "Delivery: summed over slices, ÷ the animal's mean target. Number = target. "
             "Dotted: perfect prediction.", ha="left", va="bottom", fontsize=8, color=MUTED)
    return fig


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------


def write_readme(res: pd.DataFrame, scr: pd.DataFrame, f: pd.DataFrame, d: pd.DataFrame) -> None:
    def model_table(outcome):
        rs = res[res.outcome == outcome]
        lines = ["| Model | Features | Q² (LOO) | RMSE (LOO) | ρ (LOO) | Q² leave-mouse-out | perm. p |",
                 "|---|---|---|---|---|---|---|"]
        for r in rs.itertuples():
            pp = (f"{r.perm_p_gain_over_2f:.2f} (gain over 2f)" if np.isfinite(r.perm_p_gain_over_2f)
                  else f"{r.perm_p_q2:.2f}" if np.isfinite(r.perm_p_q2) else "")
            lines.append(f"| {r.model} | {r.n_features} | {r.q2_loo:+.2f} | {r.rmse_loo:.2f} | "
                         f"{r.rho_loo:+.2f} | {r.q2_leave_mouse_out:+.2f} | {pp} |")
        return "\n".join(lines)

    scr_lines = ["| Feature | What | ρ with 2f AUC (36 targets) | ρ with delivery A / B (11) | "
                 "Q² gain added to 2f AUC, A / B |", "|---|---|---|---|---|"]
    for r in scr.itertuples():
        g = (f"{r.q2_gain_A:+.2f} / {r.q2_gain_B:+.2f}" if np.isfinite(r.q2_gain_A) else "—")
        scr_lines.append(f"| `{r.feature}` | {r.what} | {r.rho_with_2f_auc_36:+.2f} | "
                         f"{r.rho_delivery_A:+.2f} / {r.rho_delivery_B:+.2f} | {g} |")

    text = f"""# Full acoustic features vs 2nd-harmonic AUC

Generated by `scripts/acoustic_features_vs_delivery.py`; features from `src/fus/features.py`.

**Question:** does the whole acoustic record predict delivery better than the cumulative
2nd-harmonic AUC alone?

**Data:** {len(f)} targets of mice 1–6 have features; {len(d)} of them (mice 1–2, sonicated
targets) have delivery. Mouse 1 T3, the no-FUS control, has no recording and is left out, so
these numbers can't lean on it. Delivery is summed over slices and divided by the animal's mean
target (`results/histology/slices/`). Features are raw, as the dose is in the pooled dose figure.

**Method:** every model is scored on leave-one-target-out predictions (Q² = 1 − PRESS/SS; 0 = no
better than predicting the mean). Ridge picks its penalty by an inner leave-one-out inside each
fold. The permutation test reruns all of it on {N_PERM} shuffles of delivery.

## Result

Pipeline A (GFP+ area), primary:

{model_table("gfp_area_mm2_rel")}

Pipeline B (GFP+ cells):

{model_table("gfp_cells_rel")}

Figure: **`model_comparison.png`**, with Q² per model and predicted-vs-observed for 2f AUC and ridge.

## Feature screen (exploratory)

ρ with 2f AUC uses all 36 targets of mice 1–6, so it says which features restate the dose and
which carry something else, without touching delivery. The last two columns use the 11 labelled
targets. With 14 candidates, one or two will look useful by chance; treat a gain here as a
hypothesis for the next animals, not a finding.

{chr(10).join(scr_lines)}

## Read with care

- **n = 11 from two animals.** Leave-one-out Q² on 11 points is noisy: a single target can
  move it by 0.2. The permutation p treats targets as exchangeable, which they aren't quite.
- **Leave-mouse-out** fits on one animal (5–6 targets) and predicts the other. It is the
  honest test of generalisation and is expected to be poor at this size.
- **Mouse 1 T1** combines three recordings at one position; its burst-level features run over
  the concatenated bursts. Its cumulative features are the largest in the set, so it has
  leverage in every fit.
- **Feature list fixed before the comparison** (`fus.features.FEATURES`). Adding features after
  seeing this table would need the same out-of-sample treatment.
- **Broadband** is measured between the harmonic lines from the raw waveforms, not the lab's
  2.2 kHz window at 1.7 MHz. Band energies are folds over the zero-drive baseline bursts (the
  receiver noise floor), as in the lab's normalisation; they don't correct for the skull.
"""
    (OUT / "README.md").write_text(text)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    f = features()
    d = labelled(f)
    print(f"{len(f)} targets with features, {len(d)} with delivery")
    res, preds = [], []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for outcome in OUTCOMES:
            r, p = compare(d, outcome)
            res.append(r)
            preds.append(p)
        res, preds = pd.concat(res, ignore_index=True), pd.concat(preds, ignore_index=True)
        scr = screen(d, f)
    res.round(4).to_csv(OUT / "model_comparison.csv", index=False)
    preds.round(4).to_csv(OUT / "loo_predictions.csv", index=False)
    scr.round(4).to_csv(OUT / "feature_screen.csv", index=False)
    fig = figure(res, preds)
    fig.savefig(OUT / "model_comparison.png", dpi=170)
    plt.close(fig)
    write_readme(res, scr, f, d)
    print(res.round(3).to_string(index=False))
    print(scr.drop(columns="what").round(2).to_string(index=False))


if __name__ == "__main__":
    main()
