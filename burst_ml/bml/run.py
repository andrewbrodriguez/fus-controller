"""Run the pre-registered analysis and write CSVs to ``burst_ml/results/``.

From ``burst_ml/``::

    .venv/bin/python -m bml.run                      # every model, LOAO + full-data fits
    .venv/bin/python -m bml.run --models m3a m4 --quick

Steps, each skipped for a model already in its output (``--force`` redoes it):

1. **LOAO** for each model: ``loao.loao`` → ``loao_folds.csv`` and ``loao_targets.csv``
   (saved after every model, so a killed run resumes).
2. **Comparison**: ``loao.compare`` against m3a and against m1 → ``comparison.csv`` (column
   ``base``); ``loao.calibration`` → ``calibration.csv``.
3. **Full-data fits** of every model on all 35 targets → ``posterior_summary.csv`` (model,
   param, index, mean, sd, q05, q50, q95, plus ``diag/max_rhat`` and ``diag/divergences``
   rows; ``z`` skipped) and ``cache/samples_<model>.npz`` for the figures.
4. **Dose inversion** from the full-data m3a and m3b samples → ``dose_for_delivery.csv``:
   the dose for 20% GFP+ (both) and for 90% of the plateau (m3b), in a new animal.

``--exclude 1:1 --out DIR`` drops those mouse:target rows (exploratory runs only; samples go
to ``DIR`` too). ``--quick`` (2 chains, 300 + 300 draws, n_a 50) is for smoke tests: it writes to
``cache/quick/`` unless ``--out`` is given, so quick numbers never land in ``results/``.
Defaults are the README's: 4 chains, 1000 + 1000, n_a 200, seed 0.
"""

from __future__ import annotations

import numpyro

numpyro.set_host_device_count(4)   # must run before anything creates a jax array

import argparse  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from . import CACHE, RESULTS, invert  # noqa: E402
from .data import load_design, subset  # noqa: E402
from .loao import calibration, compare, fit, loao  # noqa: E402
from .models import MODELS  # noqa: E402

ORDER = ["m1", "m2", "m3a", "m3b", "m4", "m5a", "m5b"]


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def _append(path: Path, new: pd.DataFrame, model: str) -> None:
    old = pd.read_csv(path) if path.exists() else pd.DataFrame()
    if len(old):
        old = old[old["model"] != model]
    pd.concat([old, new], ignore_index=True).to_csv(path, index=False)


def _has(path: Path, model: str) -> bool:
    return path.exists() and model in set(pd.read_csv(path)["model"])


def step_loao(models, d, out, fit_kw, n_a, force):
    for m in models:
        if _has(out / "loao_folds.csv", m) and not force:
            log(f"  {m}: LOAO done, skipping")
            continue
        log(f"  {m}: LOAO ...")
        folds, targets = loao(MODELS[m], d, seed=0, n_a=n_a, **fit_kw)
        _append(out / "loao_folds.csv", folds, m)
        _append(out / "loao_targets.csv", targets, m)
        log(f"  {m}: elpd {folds.elpd.sum():.2f}, max R-hat {folds.max_rhat.max():.3f}, "
            f"divergences {folds.divergences.sum()}")


def step_compare(out):
    folds = pd.read_csv(out / "loao_folds.csv")
    blocks = [compare(folds, base=b).assign(base=b) for b in ("m3a", "m1")
              if b in set(folds["model"])]
    pd.concat(blocks, ignore_index=True).to_csv(out / "comparison.csv", index=False)
    calibration(pd.read_csv(out / "loao_targets.csv")).to_csv(out / "calibration.csv",
                                                              index=False)


def _summary_rows(model, samples, diag):
    rows = []
    for site, v in samples.items():
        if site == "z":
            continue
        flat = np.asarray(v).reshape(len(v), -1)
        for j in range(flat.shape[1]):
            c = flat[:, j]
            rows.append({"model": model, "param": site, "index": j, "mean": c.mean(),
                         "sd": c.std(), "q05": np.quantile(c, 0.05),
                         "q50": np.quantile(c, 0.5), "q95": np.quantile(c, 0.95)})
    for key in ("max_rhat", "divergences"):
        x = float(diag[key])
        rows.append({"model": model, "param": f"diag/{key}", "index": 0, "mean": x, "sd": 0.0,
                     "q05": x, "q50": x, "q95": x})
    return rows


def step_full(models, d, out, cache, fit_kw, force):
    path = out / "posterior_summary.csv"
    for m in models:
        npz = cache / f"samples_{m}.npz"
        if _has(path, m) and npz.exists() and not force:
            log(f"  {m}: full fit done, skipping")
            continue
        samples, diag = fit(MODELS[m], d, seed=0, **fit_kw)
        cache.mkdir(parents=True, exist_ok=True)
        np.savez(npz, **samples)
        _append(path, pd.DataFrame(_summary_rows(m, samples, diag)), m)
        log(f"  {m}: full fit, max R-hat {diag['max_rhat']:.3f}, "
            f"divergences {diag['divergences']}")


def step_invert(out, cache, n_a):
    def load(m):
        return dict(np.load(cache / f"samples_{m}.npz"))
    rows = []
    for m in ("m3a", "m3b"):
        if (cache / f"samples_{m}.npz").exists():
            rows.append(invert.dose_for_level(m, load(m), level=0.2, n_a=n_a))
    if (cache / "samples_m3b.npz").exists():
        rows.append(invert.dose_at_plateau_fraction(load("m3b"), frac=0.9, n_a=n_a))
    pd.DataFrame(rows).to_csv(out / "dose_for_delivery.csv", index=False)


def _exclude_pairs(pairs: list[str]) -> list[tuple[int, int]]:
    """Parse ``['1:1', '2:3']`` → ``[(1, 1), (2, 3)]``."""
    out: list[tuple[int, int]] = []
    for p in pairs:
        m_str, t_str = p.split(":")
        out.append((int(m_str), int(t_str)))
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--models", nargs="*", default=ORDER)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--out", default=None, help="results directory (default burst_ml/results)")
    ap.add_argument("--exclude", nargs="*", default=[],
                      help="Mouse:Target pairs to drop from the design (e.g. 1:1 2:3)")
    a = ap.parse_args(argv)
    bad = [m for m in a.models if m not in MODELS]
    if bad:
        ap.error(f"unknown models: {bad}")
    if a.exclude and not a.out:
        ap.error("--exclude needs --out, so the pre-registered results and samples aren't touched")

    out = Path(a.out) if a.out else (CACHE / "quick" if a.quick else RESULTS)
    cache = out if (a.quick or a.exclude) else CACHE
    out.mkdir(parents=True, exist_ok=True)
    fit_kw = dict(chains=2, warmup=300, draws=300) if a.quick else dict(chains=4, warmup=1000,
                                                                        draws=1000)
    n_a = 50 if a.quick else 200

    d = load_design()
    if a.exclude:
        pairs = _exclude_pairs(a.exclude)
        keep = np.ones(len(d["k"]), dtype=bool)
        for m, t in pairs:
            keep = keep & ~((d["mouse"] == m) & (d["target"] == t))
        d = subset(d, np.flatnonzero(keep))
        log(f"excluded {len(pairs)} target(s): {pairs}")
    log(f"{len(d['k'])} targets, {d['n_animals']} animals; {fit_kw}, n_a {n_a}; out {out}")
    log("1. LOAO")
    step_loao(a.models, d, out, fit_kw, n_a, a.force)
    log("2. comparison")
    step_compare(out)
    log("3. full-data fits")
    step_full(a.models, d, out, cache, fit_kw, a.force)
    log("4. dose inversion")
    step_invert(out, cache, n_a)
    log("done")


if __name__ == "__main__":
    main()
