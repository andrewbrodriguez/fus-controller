"""Run the pre-registered analysis and write CSVs to ``burst_ml/results/``.

From ``burst_ml/``::

    .venv/bin/python -m bml.run                      # every model, LOAO + full-data fits
    .venv/bin/python -m bml.run --models m3a m4 --quick

Steps, each skipped if its output already holds that model (``--force`` redoes it):

1. **LOAO** for each model: ``loao.loao`` → append to ``loao_folds.csv`` and
   ``loao_targets.csv`` (saved after every model, so a killed run resumes).
2. **Comparison**: ``loao.compare`` → ``comparison.csv`` (vs m3a, and a second block vs m1
   with column ``base``); ``loao.calibration`` → ``calibration.csv``.
3. **Full-data fits** of every model on all 35 targets → ``posterior_summary.csv``: model,
   param, index (for vector params, else 0), mean, sd, q05, q50, q95; plus rows
   ``diag/max_rhat`` and ``diag/divergences`` per model. Skip ``z``.
4. **Dose inversion** from the full-data m3a and m3b samples → ``dose_for_delivery.csv``:
   ``invert.dose_for_level`` at level 0.2 for both, and ``invert.dose_at_plateau_fraction``
   at 0.9 for m3b.

``--quick`` = 2 chains, 300 warm-up, 300 draws, n_a 50 (for smoke tests only; never report
quick numbers). Defaults are the README's: 4 chains, 1000 + 1000, n_a 200, seed 0.

Call ``numpyro.set_host_device_count(4)`` before anything imports jax.
"""

from __future__ import annotations

import argparse

ORDER = ["m1", "m2", "m3a", "m3b", "m4", "m5a", "m5b"]


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--models", nargs="*", default=ORDER)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--out", default=None, help="results directory (default burst_ml/results)")
    a = ap.parse_args(argv)
    raise NotImplementedError


if __name__ == "__main__":
    main()
