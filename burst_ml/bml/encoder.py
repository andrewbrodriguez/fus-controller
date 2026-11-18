"""Exploratory: one encoder over every burst spectrum, then delivery from its latent.

From ``burst_ml/``::

    .venv/bin/python -m bml.encoder            # cache spectra (~5 min first run), fit, score
    .venv/bin/python -m bml.encoder --quick    # smoke test, writes to cache/quick/encoder/
    .venv/bin/python -m bml.encoder --report   # README.md from the CSVs, no refit

**Not pre-registered.** The design below was written before any model saw delivery and is run
once as written; nothing is tuned against delivery. It is still a post hoc idea tested after
the burst models' verdict, so a win here would need a fresh pre-registered test.

**Step 1, label-free (all 24 mice, no delivery data).** Every sonicating burst's magnitude
spectrum (``fus.extract``, 131,072 bins to 2.5 MHz) is cut into 80 equal bands over
0.2-2.5 MHz. Each band gives two numbers, both as log fold over the mean of that recording's
pre-bubble baseline bursts: the band's **sum** (broadband energy) and its **peak** (a narrow line
such as a harmonic). With log drive voltage that is 161 inputs per burst, z-scored over all
bursts. Two encoders, each with a **3-number latent**, are fitted to all ~18,000 bursts:

* ``ae``: autoencoder, 161 -> 64 -> 3 -> 64 -> 161, tanh, Adam (scikit-learn
  ``MLPRegressor``, ``random_state=0``, early stopping on a 10% split)
* ``pca``: the linear version, first 3 principal components

**Step 2, per target.** A target's features are the mean of each latent number over its bursts
plus log(number of bursts), z-scored over all 145 recordings. Mouse 1 T1 joins its three
recordings, as in ``bml.data``.

**Step 3, delivery.** The burst_ml likelihood, animal effect, held-out-animal scoring and
decision rule (README "Fitting and scoring"), on the 35 targets of mice 1-6:

* ``m3a``: dose only (re-fitted on the current labels, so the comparison is like for like)
* ``ae_only``, ``pca_only``: ``logit p = beta0 + a_j + w . u`` with no dose, ``w ~ N(0, 1)``
* ``m3a_ae``, ``m3a_pca``: ``m3a + w . u``, ``w ~ N(0, 0.25)`` (as the band weights in M5b)

An encoder model **beats the dose** only if ``delta elpd vs m3a > 0`` and ``delta / SE > 2``.
Four models are tested against one baseline, so even that bar is loose.

Writes ``results/exploratory_encoder/``: ``comparison.csv``, ``calibration.csv``,
``loao_folds.csv``, ``loao_targets.csv``, ``target_features.csv``, ``encoder_checks.csv`` and
a generated ``README.md``. The spectra cache is ``cache/spectra_bands.npz`` (git-ignored).
"""

from __future__ import annotations

import numpyro

numpyro.set_host_device_count(4)   # before any jax array exists

import argparse  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import numpyro.distributions as dist  # noqa: E402

from . import CACHE, RESULTS  # noqa: E402
from .data import load_labels, recordings_for  # noqa: E402
from .loao import calibration, compare, fit, predict  # noqa: E402
from .models import MODELS, Spec, _intercept_slope  # noqa: E402

SPECTRA = CACHE / "spectra_bands.npz"
OUT = RESULTS / "exploratory_encoder"
EDGES = np.linspace(0.2e6, 2.5e6, 81)
LATENT = 3


# --- step 1: spectra -> bands, and the encoders -------------------------------------------

def _bands(member) -> dict:
    import fus.archive as A
    from fus.extract import extract
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with A.local_path(member) as path:
            r = extract(path).recording
    idx = np.searchsorted(r.freq_axis, EDGES)
    S = r.spectra
    sums = np.stack([S[:, a:b].sum(1) for a, b in zip(idx[:-1], idx[1:])], 1)
    peaks = np.stack([S[:, a:b].max(1) for a, b in zip(idx[:-1], idx[1:])], 1)
    nb = r.n_baseline
    f_sum = np.log(np.clip(sums[nb:], 1e-30, None) / sums[:nb].mean(0))
    f_peak = np.log(np.clip(peaks[nb:], 1e-30, None) / peaks[:nb].mean(0))
    volt = (r.drive_voltage[nb:, 0] if r.drive_voltage is not None
            else np.full(len(S) - nb, np.nan))
    return {"x": np.hstack([f_sum, f_peak]), "volt": volt}


def build_spectra(path=SPECTRA) -> None:
    """One row per sonicating burst of every recording in the archive (all 24 mice)."""
    import fus.archive as A
    members = A.recordings()
    xs, volts, mice, recs, ts = [], [], [], [], []
    t0 = time.time()
    for i, m in enumerate(members, 1):
        b = _bands(m)
        n = len(b["x"])
        xs.append(b["x"])
        volts.append(b["volt"])
        mice.append(np.full(n, m.mouse))
        recs.append(np.full(n, m.recording, dtype=object))
        ts.append(np.arange(n))
        print(f"  [{i}/{len(members)}] mouse {m.mouse} {m.recording}: {n} bursts "
              f"({time.time() - t0:.0f} s)", flush=True)
    np.savez_compressed(path, x=np.vstack(xs), volt=np.concatenate(volts),
                        mouse=np.concatenate(mice), recording=np.concatenate(recs).astype(str),
                        t=np.concatenate(ts), edges=EDGES)


def load_spectra(path=SPECTRA) -> dict:
    z = np.load(path, allow_pickle=False)
    return {k: z[k] for k in z.files}


def encoder_inputs(s: dict) -> np.ndarray:
    v = np.log(np.clip(s["volt"], 1e-4, None))
    X = np.column_stack([s["x"], v])
    X = np.where(np.isfinite(X), X, np.nan)
    X = np.where(np.isnan(X), np.nanmean(X, 0), X)
    return (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1.0)


def fit_encoders(X: np.ndarray, quick: bool = False) -> tuple[dict, list[dict]]:
    """Latent (n_bursts, LATENT) per encoder, plus label-free checks."""
    from sklearn.decomposition import PCA
    from sklearn.neural_network import MLPRegressor

    checks = []
    pca = PCA(LATENT, random_state=0).fit(X)
    z_pca = pca.transform(X)
    checks.append({"encoder": "pca", "check": "variance explained (all inputs)",
                   "value": float(pca.explained_variance_ratio_.sum())})

    ae = MLPRegressor(hidden_layer_sizes=(64, LATENT, 64), activation="tanh", solver="adam",
                      max_iter=40 if quick else 400, early_stopping=True,
                      validation_fraction=0.1, n_iter_no_change=20, random_state=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ae.fit(X, X)
    h = X
    for W, b in zip(ae.coefs_[:2], ae.intercepts_[:2]):
        h = np.tanh(h @ W + b)
    recon = ae.predict(X)
    checks.append({"encoder": "ae", "check": "variance explained (all inputs)",
                   "value": float(1 - ((X - recon) ** 2).sum() / ((X - X.mean(0)) ** 2).sum())})
    checks.append({"encoder": "ae", "check": "best validation R2 (sklearn)",
                   "value": float(ae.best_validation_score_)})
    checks.append({"encoder": "ae", "check": "epochs", "value": float(ae.n_iter_)})
    return {"ae": h, "pca": z_pca}, checks


# --- step 2: per-target features ----------------------------------------------------------

def recording_features(s: dict, Z: dict) -> pd.DataFrame:
    """Per recording: mean of each latent number and n bursts (all 145 recordings)."""
    key = pd.DataFrame({"mouse": s["mouse"], "recording": s["recording"]})
    cols = {}
    for enc, z in Z.items():
        for j in range(LATENT):
            cols[f"{enc}_z{j + 1}"] = z[:, j]
    df = pd.concat([key, pd.DataFrame(cols)], axis=1)
    g = df.groupby(["mouse", "recording"], sort=False)
    out = g.sum()
    out["n_bursts"] = g.size()
    return out.reset_index()   # sums; divided by n after recordings are joined per target


def target_rows(rec: pd.DataFrame, mouse: int, target: int) -> pd.Series:
    r = rec[(rec.mouse == mouse) & rec.recording.isin(recordings_for(mouse, target))]
    s = r.drop(columns=["mouse", "recording"]).sum()
    n = s["n_bursts"]
    out = s.drop("n_bursts") / n
    out["log_n_bursts"] = np.log(n)
    return out


def standardise(rec: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Mean and SD of each target feature over all 145 recordings (label-free constants)."""
    per = rec.drop(columns=["mouse", "recording"]).div(rec.n_bursts, axis=0)
    per["log_n_bursts"] = np.log(rec.n_bursts)
    per = per.drop(columns="n_bursts")
    return per.mean(), per.std(ddof=0)


# --- step 3: delivery ---------------------------------------------------------------------

def prob_only(th, d, a):
    import jax.numpy as jnp
    from jax.scipy.special import expit
    return expit(th["beta0"] + a + jnp.asarray(d["U"]) @ th["w"])


def prob_plus(th, d, a):
    import jax.numpy as jnp
    from jax.scipy.special import expit
    return expit(th["beta0"] + a + th["beta1"] * jnp.log(d["dose"]) + jnp.asarray(d["U"]) @ th["w"])


def specs(n_u: int) -> dict[str, Spec]:
    b0 = _intercept_slope()["beta0"]
    return {
        "m3a": MODELS["m3a"],
        "ae_only": Spec("ae_only", prob_only, {"beta0": b0, "w": dist.Normal(0, 1).expand([n_u])}),
        "pca_only": Spec("pca_only", prob_only, {"beta0": b0, "w": dist.Normal(0, 1).expand([n_u])}),
        "m3a_ae": Spec("m3a_ae", prob_plus, {**_intercept_slope(),
                                             "w": dist.Normal(0, 0.25).expand([n_u])}),
        "m3a_pca": Spec("m3a_pca", prob_plus, {**_intercept_slope(),
                                               "w": dist.Normal(0, 0.25).expand([n_u])}),
    }


def design(labels: pd.DataFrame, feats: pd.DataFrame, enc: str | None) -> dict:
    cols = [f"{enc}_z{j + 1}" for j in range(LATENT)] + ["log_n_bursts"] if enc else []
    mice = sorted(labels.mouse.unique())
    d = {"mouse": labels.mouse.to_numpy(np.int64), "target": labels.target.to_numpy(np.int64),
         "animal": labels.mouse.map({m: i for i, m in enumerate(mice)}).to_numpy(np.int64),
         "n_animals": len(mice), "k": labels.k.to_numpy(np.int64),
         "n": labels.n.to_numpy(np.int64), "dose": labels.dose.to_numpy(float),
         "U": (feats.loc[list(zip(labels.mouse, labels.target)), cols].to_numpy(float)
               if enc else np.zeros((len(labels), 0)))}
    return d


def subset(d: dict, rows: np.ndarray) -> dict:
    out = {k: v[rows] for k, v in d.items() if isinstance(v, np.ndarray)}
    first = {a: i for i, a in enumerate(dict.fromkeys(out["animal"].tolist()))}
    out["animal"] = np.array([first[a] for a in out["animal"]], dtype=np.int64)
    out["n_animals"] = len(first)
    return out


def loao(spec: Spec, d: dict, fit_kw: dict, n_a: int, seed: int = 0):
    folds, targets = [], []
    for mouse in sorted(set(d["mouse"])):
        hold = d["mouse"] == mouse
        samples, diag = fit(spec, subset(d, ~hold), seed=seed, **fit_kw)
        dt = subset(d, hold)
        pred = predict(spec, samples, dt, n_a=n_a, seed=seed + int(mouse))
        folds.append({"model": spec.name, "mouse": int(mouse), "elpd": pred["joint"], **diag})
        for i in range(len(dt["k"])):
            f = pred["frac"][:, i]
            obs = dt["k"][i] / dt["n"][i]
            targets.append({"model": spec.name, "mouse": int(mouse), "target": int(dt["target"][i]),
                            "dose": float(dt["dose"][i]), "k": int(dt["k"][i]), "n": int(dt["n"][i]),
                            "frac": obs, "lpd": float(pred["pointwise"][i]),
                            "pred_median": float(np.median(f)),
                            "q05": float(np.percentile(f, 5)), "q25": float(np.percentile(f, 25)),
                            "q75": float(np.percentile(f, 75)), "q95": float(np.percentile(f, 95)),
                            "pit": float((f < obs).mean() + 0.5 * (f == obs).mean())})
    return pd.DataFrame(folds), pd.DataFrame(targets)


# --- report -------------------------------------------------------------------------------

def _md(df: pd.DataFrame, nd: int = 2) -> str:
    df = df.copy()
    for c in df.columns:
        if df[c].dtype.kind == "f":
            df[c] = df[c].map(lambda v: "" if pd.isna(v) else f"{v:.{nd}f}")
    head = "| " + " | ".join(df.columns) + " |"
    rule = "|" + "|".join("---" for _ in df.columns) + "|"
    return "\n".join([head, rule] + ["| " + " | ".join(map(str, r)) + " |" for r in df.values])


def write_readme(out, comp, cal, folds, checks, n_bursts, n_rec) -> None:
    per = folds.pivot(index="mouse", columns="model", values="elpd")
    per = per.sub(per["m3a"], axis=0).drop(columns="m3a").reset_index()
    diag = (folds.groupby("model").agg(max_rhat=("max_rhat", "max"),
                                       divergences=("divergences", "sum"),
                                       draws_per_fold=("draws", "first"))
            .reset_index())
    flagged = diag[(diag.max_rhat > 1.01) | (diag.divergences > 0.01 * diag.draws_per_fold)]
    flag = ("Every fit passes the burst_ml bar (R-hat <= 1.01, at most 1% divergent draws)."
            if flagged.empty else
            "Flagged by the burst_ml bar (R-hat > 1.01 or more than 1% divergent draws): "
            + ", ".join(flagged.model) + ". Folds: see `loao_folds.csv`.")
    text = f"""# Exploratory: one encoder over every burst spectrum

Generated by `bml/encoder.py` (design in its docstring). **Not pre-registered**; run once as
written. Labels: `results/histology/slices/cumulative_dose_delivery.csv` as of this run.

Encoders were fitted to {n_bursts:,} sonicating bursts from {n_rec} recordings in all 24 mice,
with no delivery data. Delivery models are scored on each held-out mouse (6 folds, 35 targets).

## Against dose alone (m3a)

A model beats the dose only if delta > 0 and delta / SE > 2.

{_md(comp[comp.base == "m3a"].drop(columns="base"))}

## Against no acoustics at all (m1)

{_md(comp[comp.base == "m1"].drop(columns="base"))}

## Per held-out mouse, elpd minus m3a's

{_md(per)}

## Calibration (held-out)

{_md(cal)}

## Encoder checks (label-free)

{_md(checks, 3)}

## Fit diagnostics (6 held-out folds per model)

{_md(diag, 3)}

{flag}
"""
    (out / "README.md").write_text(text)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--rebuild", action="store_true", help="re-read every recording")
    ap.add_argument("--report", action="store_true", help="rewrite README.md from the CSVs only")
    args = ap.parse_args(argv)
    out = CACHE / "quick/encoder" if args.quick else OUT
    out.mkdir(parents=True, exist_ok=True)
    if args.report:
        s = load_spectra()
        write_readme(out, pd.read_csv(out / "comparison.csv"), pd.read_csv(out / "calibration.csv"),
                     pd.read_csv(out / "loao_folds.csv"), pd.read_csv(out / "encoder_checks.csv"),
                     len(s["mouse"]), len(set(zip(s["mouse"].tolist(), s["recording"].tolist()))))
        print(f"wrote {out / 'README.md'}")
        return
    fit_kw = (dict(chains=2, warmup=300, draws=300) if args.quick
              else dict(chains=4, warmup=1000, draws=1000))
    n_a = 50 if args.quick else 200

    if args.rebuild or not SPECTRA.exists():
        print("reading spectra from the archive ...", flush=True)
        build_spectra()
    s = load_spectra()
    X = encoder_inputs(s)
    print(f"{len(X):,} bursts x {X.shape[1]} inputs; fitting encoders ...", flush=True)
    Z, checks = fit_encoders(X, quick=args.quick)

    rec = recording_features(s, Z)
    mu, sd = standardise(rec)
    labels = load_labels()
    keys = sorted(set(zip(labels.mouse, labels.target)))
    # all targets of mice 7-24 too, for the label-free "does the latent hold the dose" check
    feats = pd.DataFrame({k: target_rows(rec, *k) for k in keys}).T
    feats = (feats - mu) / sd
    feats.index = pd.MultiIndex.from_tuples(feats.index, names=["mouse", "target"])

    # label-free checks: how much of log dose and log mean voltage each latent carries
    per_rec = rec.drop(columns=["mouse", "recording"]).div(rec.n_bursts, axis=0)
    per_rec["log_n_bursts"] = np.log(rec.n_bursts)
    rd = pd.read_csv(CACHE / "recordings.csv").set_index(["mouse", "recording"])
    log_dose = np.log(rd.loc[list(zip(rec.mouse, rec.recording)), "dose"].to_numpy())
    for enc in Z:
        A = np.column_stack([np.ones(len(rec))]
                            + [per_rec[f"{enc}_z{j + 1}"] for j in range(LATENT)]
                            + [per_rec.log_n_bursts])
        resid = log_dose - A @ np.linalg.lstsq(A, log_dose, rcond=None)[0]
        checks.append({"encoder": enc, "check": "R2 of log dose from target features (145 recordings)",
                       "value": float(1 - resid.var() / log_dose.var())})
    checks = pd.DataFrame(checks)
    print(checks.round(3).to_string(index=False), flush=True)

    fe = feats.reset_index()
    fe.round(4).to_csv(out / "target_features.csv", index=False)
    checks.round(4).to_csv(out / "encoder_checks.csv", index=False)

    folds, targets = [], []
    for name, spec in specs(LATENT + 1).items():
        enc = None if name == "m3a" else name.split("_")[-1] if name.startswith("m3a_") \
            else name.split("_")[0]
        d = design(labels, feats, enc)
        t0 = time.time()
        f, t = loao(spec, d, fit_kw, n_a)
        print(f"  {name}: elpd {f.elpd.sum():.2f}, max R-hat {f.max_rhat.max():.3f}, "
              f"divergences {int(f.divergences.sum())} ({time.time() - t0:.0f} s)", flush=True)
        folds.append(f)
        targets.append(t)
    m1_f, m1_t = loao(MODELS["m1"], design(labels, feats, None), fit_kw, n_a)
    folds = pd.concat(folds + [m1_f], ignore_index=True)
    targets = pd.concat(targets + [m1_t], ignore_index=True)
    comp = pd.concat([compare(folds, "m3a").assign(base="m3a"),
                      compare(folds, "m1").assign(base="m1")], ignore_index=True)
    cal = calibration(targets)
    folds.round(4).to_csv(out / "loao_folds.csv", index=False)
    targets.round(5).to_csv(out / "loao_targets.csv", index=False)
    comp.round(4).to_csv(out / "comparison.csv", index=False)
    cal.round(4).to_csv(out / "calibration.csv", index=False)
    write_readme(out, comp, cal, folds, checks, len(X), rec.shape[0])
    print(comp.round(2).to_string(index=False))
    print(cal.round(2).to_string(index=False))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
