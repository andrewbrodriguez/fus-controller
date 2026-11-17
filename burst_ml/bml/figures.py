"""Aggregate figures from ``burst_ml/results/*.csv`` and ``cache/samples_*.npz``.

From ``burst_ml/``: ``.venv/bin/python -m bml.figures``. Writes PNGs (dpi 150) to
``burst_ml/results/``. Every figure pools all six mice. Use matplotlib's Agg backend,
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


def main() -> None:
    raise NotImplementedError


if __name__ == "__main__":
    main()
