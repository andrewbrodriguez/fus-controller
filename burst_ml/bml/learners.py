"""Exploratory: standard supervised learners, sonication record -> delivery, per held-out mouse.

From ``burst_ml/``::

    .venv/bin/python -m bml.learners      # needs cache/spectra_bands.npz (``bml.encoder``)

**Not pre-registered.** Every learner and setting below was fixed before the first run, and
every result is reported, not only the best. Many learners are tried, so the best of them is
optimistic by construction; the question is whether *any* beats the dose line clearly.

**Outcome.** logit of ΣGFP+ ÷ ΣNeuN for the 35 sonicated targets of mice 1-6 (the burst_ml
labels), clipped to [0.001, 0.999] first.

**Feature sets** (per target; Mouse 1 T1 joins its three recordings, as in ``bml.data``):

* ``dose``: log dose
* ``engineered``: log dose, log n bursts, drive voltage (mean while on, max, log sum), log
  ``h2`` (mean, SD, max, first-third mean, last-third mean), mean log ``wb``, ``bb``, ``sub``,
  ``ultra``, ``f0``, max log ``bb``, and the 3 + 3 encoder latents (``bml.encoder``) -- 23
* ``spectrum``: the mean over the target's bursts of each of the encoder's 161 z-scored inputs
  (80 band sums, 80 band peaks, log voltage), plus log n bursts -- 162

**Learners** (scikit-learn, inputs z-scored inside each training fold, fixed seeds):
``linear`` (least squares), ``ridge`` (RidgeCV, alphas 1e-3..1e3), ``lasso`` (LassoCV, 5-fold),
``forest`` (500 trees, min leaf 3), ``boosting`` (200 stumps-of-depth-2, rate 0.05, subsample
0.8), ``mlp`` (32-16 tanh, L2 1e-2, lbfgs), ``gp`` (constant × RBF + white noise), ``knn``
(k = 5). ``linear`` is run on ``dose`` only (with 23 or 162 inputs and 30 rows it is
undetermined or nearly so); every other learner on all three sets.

**Scoring, hold out each mouse (6 folds).** Train on five mice, predict the sixth's targets.
``rmse``: root mean square error of the logit over all 35 held-out targets. ``r2``: 1 − SSE ÷ SSE
of predicting the training mean (no acoustics). ``spearman``: mean over held-out mice of the
within-mouse rank correlation, predicted vs observed. ``mice_better``: in how many held-out
mice the squared error is lower than ``linear`` on ``dose``.

Added after a first run crashed at the MLP (an argument error; the dose-only rows had printed):
``rmse_vs_own_dose``, each learner's RMSE minus the same learner's on ``dose`` alone, since a
flexible learner can beat the dose *line* just by bending the dose curve.

Writes ``results/exploratory_learners/``: ``scores.csv``, ``predictions.csv``, ``README.md``.
"""

from __future__ import annotations

import argparse
import warnings

import numpy as np
import pandas as pd
from scipy.special import logit
from scipy.stats import spearmanr

from . import CACHE, RESULTS
from .data import load_bursts, load_labels, recordings_for
from .encoder import encoder_inputs, load_spectra

OUT = RESULTS / "exploratory_learners"
ENC_FEATURES = RESULTS / "exploratory_encoder/target_features.csv"


def _log(x):
    return np.log(np.clip(np.asarray(x, float), 1e-12, None))


def engineered(labels: pd.DataFrame) -> pd.DataFrame:
    b = load_bursts()
    enc = pd.read_csv(ENC_FEATURES).set_index(["mouse", "target"])
    rows = []
    for m, t, dose in zip(labels.mouse, labels.target, labels.dose):
        g = pd.concat([b[(b.mouse == m) & (b.recording == r)] for r in recordings_for(m, t)])
        h = _log(g.h2)
        third = max(len(h) // 3, 1)
        on = g.volt > 0
        row = {"log_dose": np.log(dose), "log_n_bursts": np.log(len(g)),
               "volt_mean": g.volt[on].mean(), "volt_max": g.volt.max(),
               "log_volt_sum": np.log(g.volt.sum()),
               "log_h2_mean": h.mean(), "log_h2_sd": h.std(), "log_h2_max": h.max(),
               "log_h2_first": h[:third].mean(), "log_h2_last": h[-third:].mean(),
               **{f"log_{c}_mean": _log(g[c]).mean() for c in ("wb", "bb", "sub", "ultra", "f0")},
               "log_bb_max": _log(g.bb).max()}
        row.update(enc.loc[(m, t), [c for c in enc.columns if "_z" in c]].to_dict())
        rows.append(row)
    return pd.DataFrame(rows)


def spectrum(labels: pd.DataFrame) -> pd.DataFrame:
    s = load_spectra()
    X = encoder_inputs(s)
    key = pd.Series(list(zip(s["mouse"].tolist(), s["recording"].tolist())))
    rows = []
    for m, t in zip(labels.mouse, labels.target):
        sel = key.isin([(m, r) for r in recordings_for(m, t)]).to_numpy()
        rows.append(np.append(X[sel].mean(0), np.log(sel.sum())))
    cols = [f"sum_{i}" for i in range(80)] + [f"peak_{i}" for i in range(80)] + ["log_volt",
                                                                                 "log_n_bursts"]
    return pd.DataFrame(rows, columns=cols)


def learners() -> dict:
    from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
    from sklearn.linear_model import LassoCV, LinearRegression, RidgeCV
    from sklearn.neighbors import KNeighborsRegressor
    from sklearn.neural_network import MLPRegressor
    return {
        "linear": lambda: LinearRegression(),
        "ridge": lambda: RidgeCV(alphas=np.logspace(-3, 3, 13)),
        "lasso": lambda: LassoCV(cv=5, random_state=0, max_iter=50_000),
        "forest": lambda: RandomForestRegressor(n_estimators=500, min_samples_leaf=3, random_state=0),
        "boosting": lambda: GradientBoostingRegressor(n_estimators=200, max_depth=2,
                                                      learning_rate=0.05, subsample=0.8,
                                                      random_state=0),
        "mlp": lambda: MLPRegressor(hidden_layer_sizes=(32, 16), activation="tanh", alpha=1e-2, solver="lbfgs",
                                    max_iter=5000, random_state=0),
        "gp": lambda: GaussianProcessRegressor(ConstantKernel() * RBF() + WhiteKernel(),
                                               normalize_y=True, random_state=0),
        "knn": lambda: KNeighborsRegressor(n_neighbors=5),
    }


def loao(X: np.ndarray, y: np.ndarray, mouse: np.ndarray, make) -> np.ndarray:
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    pred = np.empty_like(y)
    for m in np.unique(mouse):
        test = mouse == m
        model = make_pipeline(StandardScaler(), make())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model.fit(X[~test], y[~test])
        pred[test] = model.predict(X[test])
    return pred


def baseline_mean(y, mouse):
    pred = np.empty_like(y)
    for m in np.unique(mouse):
        pred[mouse == m] = y[mouse != m].mean()
    return pred


def score(name, feats, y, pred, base_sq, ref_sq, mouse) -> dict:
    sq = (pred - y) ** 2
    with warnings.catch_warnings():   # a learner may predict one value for a whole mouse
        warnings.simplefilter("ignore")
        rho = [spearmanr(pred[mouse == m], y[mouse == m]).statistic for m in np.unique(mouse)]
    per_m = pd.Series(sq).groupby(mouse).sum()
    ref_m = pd.Series(ref_sq).groupby(mouse).sum()
    return {"learner": name, "features": feats, "rmse": float(np.sqrt(sq.mean())),
            "r2": float(1 - sq.sum() / base_sq.sum()),
            "spearman": float(np.nanmean(rho)) if np.isfinite(rho).any() else np.nan,
            "mice_better": int((per_m < ref_m - 1e-12).sum())}


def _md(df: pd.DataFrame) -> str:
    df = df.copy()
    for c in df.columns:
        if df[c].dtype.kind == "f":
            df[c] = df[c].map(lambda v: "—" if pd.isna(v) else
                              f"{v:+.2f}" if c in ("r2", "rmse_vs_own_dose") else f"{v:.2f}")
    head = "| " + " | ".join(df.columns) + " |"
    rule = "|" + "|".join("---" for _ in df.columns) + "|"
    return "\n".join([head, rule] + ["| " + " | ".join(map(str, r)) + " |" for r in df.values])


def main(argv=None) -> None:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    labels = load_labels().reset_index(drop=True)
    y = logit(np.clip(labels.k / labels.n, 1e-3, 1 - 1e-3)).to_numpy(float)
    mouse = labels.mouse.to_numpy()
    eng = engineered(labels)
    sets = {"dose": eng[["log_dose"]], "engineered": eng, "spectrum": spectrum(labels)}

    base = baseline_mean(y, mouse)
    base_sq = (base - y) ** 2
    ref = loao(sets["dose"].to_numpy(float), y, mouse, learners()["linear"])
    ref_sq = (ref - y) ** 2

    rows = [score("training mean (no acoustics)", "none", y, base, base_sq, ref_sq, mouse)]
    preds = pd.DataFrame({"mouse": mouse, "target": labels.target, "y_logit": y,
                          "none/mean": base})
    for fname, F in sets.items():
        for lname, make in learners().items():
            if lname == "linear" and fname != "dose":
                continue
            p = loao(F.to_numpy(float), y, mouse, make)
            preds[f"{fname}/{lname}"] = p
            rows.append(score(lname, fname, y, p, base_sq, ref_sq, mouse))
            print(f"  {fname:10s} {lname:9s} rmse {rows[-1]['rmse']:.2f}  r2 {rows[-1]['r2']:+.2f}  "
                  f"spearman {rows[-1]['spearman']:+.2f}  better in {rows[-1]['mice_better']}/6",
                  flush=True)
    scores = pd.DataFrame(rows)
    own = scores[scores.features == "dose"].set_index("learner").rmse
    scores["rmse_vs_own_dose"] = [r.rmse - own[r.learner] if r.learner in own.index
                                  and r.features != "dose" else np.nan
                                  for r in scores.itertuples()]
    scores.round(4).to_csv(OUT / "scores.csv", index=False)
    preds.round(4).to_csv(OUT / "predictions.csv", index=False)

    ref_row = scores[(scores.learner == "linear") & (scores.features == "dose")].iloc[0]
    more = scores[~scores.features.isin(["dose", "none"])]
    beat = more[more.rmse < ref_row.rmse]
    worse_own = int((more.rmse_vs_own_dose > 0).sum())
    text = f"""# Exploratory: standard learners, sonication record -> delivery

Generated by `bml/learners.py` (design in its docstring). **Not pre-registered**; every
learner and setting was fixed before the first run, and all of them are reported. Labels:
`results/histology/slices/cumulative_dose_delivery.csv` as of this run (35 targets, 6 mice).

Outcome: logit GFP+ fraction. Each mouse held out in turn; the reference is a straight line in
log dose (`linear`, `dose`): RMSE {ref_row.rmse:.2f}, R² {ref_row.r2:+.2f} over predicting the
training mean, within-mouse Spearman {ref_row.spearman:+.2f}.

Learners given more than the dose that beat the reference's RMSE: **{len(beat)} of
{len(more)}**{(': ' + ', '.join(beat.learner + ' (' + beat.features + ')')) if len(beat) else ''}.
Learners that do worse with more than the dose than with the dose alone
(`rmse_vs_own_dose` > 0, a column added after a crashed first run): **{worse_own} of {len(more)}**.

{_md(scores.sort_values(['features', 'rmse']))}

`mice_better`: held-out mice (of 6) where the learner's squared error is below the reference's.
"""
    (OUT / "README.md").write_text(text)
    print(scores.sort_values("rmse").round(2).to_string(index=False))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
