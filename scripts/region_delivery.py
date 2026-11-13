"""Delivery by brain region, with dose taken out first.

    .venv/bin/python scripts/region_delivery.py      # after scripts/slices_dose_delivery.py

The lab suggested (2026-11-09), from looking at the IF images, that thalamic targets get much
stronger delivery and that right anterior striatum (T1) gets less than left (T4). Regions come
from the lab's target mask (``docs/lab/README.md``): T1/T4 anterior striatum, T2/T5 lateral
striatum and cortex, T3/T6 thalamus; T1-T3 on the animal's right.

Region and dose are tangled: whichever region happened to get more dose will look stronger.
So each target's delivery is first put on a logit scale and fitted against log dose with one
offset per mouse and one common slope; the **residual** is delivery above or below what that
mouse's dose-response predicts. Regions and left/right pairs are compared on the residuals,
and the p-values come from shuffling target labels within each mouse (doses stay with their
targets), so animal differences can't drive them.

Targets used: one sonication each. Mouse 1's T1 (three sonications) and T3 (no-FUS control)
are left out of the fit and the tests, and drawn hollow where shown.

Caveats (also in the README this writes):
- The lab's suggestions came from viewing these same animals, so agreement here is not an
  independent replication. Only Mice 6+ can test them cleanly.
- "Region" is the planned target, not the tissue in the section: sections sit ~1.7-2.7 mm
  deep (``fus.reconstruct``), above most focal centres (2.9-3.7 mm).
- n is 4-5 mice per contrast.

Inputs: ``results/histology/slices/{cumulative_dose_delivery,slices_dose_delivery}.csv`` and
``results/acoustic_eda/recordings.csv``. Writes ``results/histology/region/``:
``region_delivery.csv`` (per target), ``region_tests.csv``, ``dose_by_position.csv``,
``region_delivery.png`` and ``README.md``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import slices_dose_delivery as sdd  # noqa: E402  (shared palette and styling)

ROOT = sdd.ROOT
SLICES = sdd.OUT
OUT = ROOT / "results/histology/region"

REGION = {1: "anterior striatum", 4: "anterior striatum",
          2: "lateral striatum / cortex", 5: "lateral striatum / cortex",
          3: "thalamus", 6: "thalamus"}
REGION_ORDER = ["anterior striatum", "lateral striatum / cortex", "thalamus"]
SIDE = {t: "right" if t <= 3 else "left" for t in range(1, 7)}
PAIRS = [(1, 4), (2, 5), (3, 6)]  # right vs left, same region
MEASURES = {"cells": ("gfp_fraction", "B: fraction of neurons GFP+"),
            "pixels": ("pixel_coverage", "A: GFP+ pixel coverage")}
N_PERM = 20000
EPS = 1e-3


def table() -> pd.DataFrame:
    """One row per target: dose, both cumulative delivery measures, region and side."""
    c = pd.read_csv(SLICES / "cumulative_dose_delivery.csv")
    s = pd.read_csv(SLICES / "slices_dose_delivery.csv")
    s = s[s.dropped.isna() | (s.dropped == "")].copy()
    s["gfp_mm2"] = s.pixel_coverage * s.roi_tissue_mm2
    px = (s.groupby(["mouse", "target"])
            .agg(gfp_mm2=("gfp_mm2", "sum"), roi_tissue_mm2=("roi_tissue_mm2", "sum"))
            .reset_index())
    px["pixel_coverage"] = px.gfp_mm2 / px.roi_tissue_mm2
    d = c.merge(px[["mouse", "target", "pixel_coverage", "roi_tissue_mm2"]], on=["mouse", "target"])
    d["region"] = d.target.map(REGION)
    d["side"] = d.target.map(SIDE)
    d["in_fit"] = d.n_sonications == 1
    return d


def add_residuals(d: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """logit(delivery) ~ log(dose) + mouse; residual per target, for each measure."""
    d = d.copy()
    slopes = {}
    fit = d[d.in_fit]
    for name, (col, _) in MEASURES.items():
        y = np.log(fit[col].clip(EPS, 1 - EPS) / (1 - fit[col].clip(EPS, 1 - EPS)))
        X = pd.get_dummies(fit.mouse.astype(str)).astype(float)
        X["log_dose"] = np.log(fit.dose)
        b = np.linalg.lstsq(X.values, y.values, rcond=None)[0]
        slopes[name] = b[-1]
        d[f"residual_{name}"] = np.nan
        d.loc[fit.index, f"residual_{name}"] = y - X.values @ b
    return d, slopes


def _contrasts(f: pd.DataFrame, res: str) -> dict:
    """Every tested contrast, from residuals of the in-fit targets."""
    out = {}
    for region in REGION_ORDER:
        inside = f.region == region
        out[f"{region} vs other regions"] = (f.loc[inside, res].mean() - f.loc[~inside, res].mean(), None)
    p = f.pivot(index="mouse", columns="target", values=res)
    for a, b in PAIRS:
        diff = (p[a] - p[b]).dropna() if a in p and b in p else pd.Series(dtype=float)
        out[f"T{a} (right) - T{b} (left)"] = (diff.mean(), diff)
    return out


def tests(d: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    rows = []
    for name in MEASURES:
        res = f"residual_{name}"
        f = d[d.in_fit][["mouse", "target", res]].copy()
        f["region"] = f.target.map(REGION)
        obs = _contrasts(f, res)
        hits = dict.fromkeys(obs, 0)
        for _ in range(N_PERM):
            g = f.copy()
            g["target"] = g.groupby("mouse").target.transform(lambda t: rng.permutation(t.values))
            g["region"] = g.target.map(REGION)
            for k, (v, _) in _contrasts(g, res).items():
                hits[k] += abs(v) >= abs(obs[k][0]) - 1e-12
        for k, (v, per_mouse) in obs.items():
            rows.append({
                "measure": name, "contrast": k, "mean_residual_diff": round(v, 3),
                "perm_p": round(hits[k] / N_PERM, 4),
                "mice_negative": None if per_mouse is None else int((per_mouse < 0).sum()),
                "n_mice": None if per_mouse is None else len(per_mouse),
                "per_mouse": "" if per_mouse is None else
                    "; ".join(f"M{m} {x:+.2f}" for m, x in per_mouse.items()),
            })
    t = pd.DataFrame(rows)
    t[["mice_negative", "n_mice"]] = t[["mice_negative", "n_mice"]].astype("Int64")
    return t


def dose_by_position() -> tuple[pd.DataFrame, float, pd.DataFrame]:
    """Dose per position over every mouse from 2 on (TargetN = position N), vs the imaged ones."""
    r = pd.read_csv(ROOT / "results/acoustic_eda/recordings.csv")
    r = r[r.mouse >= 2]
    imaged = r[r.mouse.isin([m for m in sdd.MICE if m >= 2])]
    t = pd.DataFrame({
        "region": [REGION[p] for p in range(1, 7)],
        "side": [SIDE[p] for p in range(1, 7)],
        "mean_dose_all": r.groupby("position").cum_h2.mean().values,
        "n_all": r.groupby("position").cum_h2.size().values,
        "mean_dose_imaged": imaged.groupby("position").cum_h2.mean().values,
    }, index=pd.Index(range(1, 7), name="position")).round(3)
    p = stats.kruskal(*[g.cum_h2 for _, g in r.groupby("position")]).pvalue
    return t.reset_index(), p, r


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def figure(d: pd.DataFrame, tt: pd.DataFrame, rec: pd.DataFrame, kw_p: float):
    sdd._style()
    C = sdd
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.9), gridspec_kw={"width_ratios": [1.35, 0.8, 1.0]})
    fig.subplots_adjust(left=0.055, right=0.985, top=0.78, bottom=0.17, wspace=0.32)
    rng = np.random.default_rng(1)

    # (a) residual by target, grouped by region
    ax = axes[0]
    C._clean(ax, "y")
    order = [1, 4, 2, 5, 3, 6]
    xpos = {t: i + (i // 2) * 0.6 for i, t in enumerate(order)}
    ax.axhline(0, color=C.AXIS, linewidth=1.2, zorder=1)
    f = d[d.in_fit]
    for t in order:
        v = f[f.target == t]
        jit = rng.uniform(-0.16, 0.16, len(v))
        ax.scatter(xpos[t] + jit, v.residual_cells, s=70, c=[C.MOUSE_COLOR[m] for m in v.mouse],
                   edgecolors=C.SURFACE, linewidths=1.5, zorder=3)
        ax.plot([xpos[t] - 0.28, xpos[t] + 0.28], [v.residual_cells.mean()] * 2, color=C.INK,
                linewidth=2, solid_capstyle="round", zorder=4)
    ax.set_xticks([xpos[t] for t in order])
    ax.set_xticklabels([f"T{t}\n{SIDE[t]}" for t in order], fontsize=9, color=C.INK_2)
    lo, hi = ax.get_ylim()
    for k, reg in enumerate(REGION_ORDER):
        mid = (xpos[order[2 * k]] + xpos[order[2 * k + 1]]) / 2
        ax.text(mid, hi + 0.05 * (hi - lo), reg, ha="center", va="bottom", fontsize=9.5, color=C.INK)
    ax.set_ylabel("Delivery vs dose-predicted (logit residual)", fontsize=10, color=C.INK)
    ax.set_title("a   Dose-adjusted delivery by target", loc="left", fontsize=11, color=C.INK, pad=26)

    # (b) T1 vs T4, raw fraction, one line per mouse
    ax = axes[1]
    C._clean(ax, "y")
    for m in sdd.MICE:
        v = d[(d.mouse == m) & d.target.isin([1, 4])].set_index("target")
        if len(v) < 2:
            continue
        hollow = not v.loc[1, "in_fit"]
        col = C.MOUSE_COLOR[m]
        ax.plot([0, 1], [v.loc[1, "gfp_fraction"], v.loc[4, "gfp_fraction"]], color=col, linewidth=2,
                linestyle="--" if hollow else "-", zorder=2)
        for x, t in ((0, 1), (1, 4)):
            ax.scatter(x, v.loc[t, "gfp_fraction"], s=80, zorder=3, linewidths=1.8,
                       facecolors=C.SURFACE if hollow else col, edgecolors=col)
        ax.annotate(f"M{m}", (1, v.loc[4, "gfp_fraction"]), xytext=(9, 0), textcoords="offset points",
                    va="center", fontsize=8.5, color=C.INK_2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["T1\nright anterior", "T4\nleft anterior"], fontsize=9, color=C.INK_2)
    ax.set_xlim(-0.3, 1.35)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    ax.set_ylabel("Neurons GFP+ (Σ GFP+ ÷ Σ NeuN)", fontsize=10, color=C.INK)
    r14 = tt[(tt.measure == "cells") & (tt.contrast.str.startswith("T1"))].iloc[0]
    ax.set_title("b   Anterior striatum, right vs left", loc="left", fontsize=11, color=C.INK, pad=26)
    ax.text(0.0, 1.02, f"dose-adjusted: T1 lower in {r14.mice_negative}/{r14.n_mice} mice, "
            f"perm p = {r14.perm_p:.3f}", transform=ax.transAxes, fontsize=8.5, color=C.INK_2)

    # (c) dose per position: every mouse vs the imaged ones
    ax = axes[2]
    C._clean(ax, "y")
    for p in range(1, 7):
        v = rec[rec.position == p]
        ax.scatter(p + rng.uniform(-0.14, 0.14, len(v)), v.cum_h2, s=22, color=C.GRID,
                   edgecolors=C.MUTED, linewidths=0.5, zorder=2)
        ax.plot([p - 0.3, p + 0.3], [v.cum_h2.mean()] * 2, color=C.MUTED, linewidth=2, zorder=3)
        im = v[v.mouse.isin([m for m in sdd.MICE if m >= 2])]
        ax.scatter([p + 0.32] * len(im), im.cum_h2, s=46, c=[C.MOUSE_COLOR[m] for m in im.mouse],
                   edgecolors=C.SURFACE, linewidths=1.2, zorder=4)
    ax.set_xticks(range(1, 7))
    ax.set_xticklabels([f"T{p}" for p in range(1, 7)], fontsize=9, color=C.INK_2)
    ax.set_ylabel("Acoustic dose (cumulative 2nd harmonic)", fontsize=10, color=C.INK)
    ax.set_title("c   Dose by target position", loc="left", fontsize=11, color=C.INK, pad=26)
    ax.text(0.0, 1.02, f"grey: all {rec.mouse.nunique()} mice (Kruskal-Wallis p = {kw_p:.2f}); "
            "colour: imaged mice", transform=ax.transAxes, fontsize=8.5, color=C.INK_2)

    handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=8, color=C.MOUSE_COLOR[m],
                          markeredgecolor=C.SURFACE, label=f"Mouse {m}") for m in sdd.MICE]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.985, 0.995), ncol=5,
               frameon=False, fontsize=9, handletextpad=0.3, columnspacing=1.2)
    C._header(fig, "Delivery by brain region, with dose taken out first",
              "Residual = logit(fraction of neurons GFP+) minus a fit on log dose with one offset per mouse. "
              "Bars = mean. Thalamus looks strongest raw because it got more dose in these mice.")
    C._footer(fig, "Pipeline B, summed over each mouse's kept slices. Mouse 1 T1 (three sonications, dashed in b) "
              "and T3 (no-FUS control) are outside the fit. p from shuffling targets within each mouse. "
              "The lab suggested both questions after viewing these same animals.")
    return fig


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------

def write_readme(d, tt, slopes, pos, kw_p):
    f = d[d.in_fit]
    reg = (f.groupby("region")
             .agg(targets=("dose", "size"), mean_dose=("dose", "mean"),
                  raw_fraction=("gfp_fraction", "mean"), raw_coverage=("pixel_coverage", "mean"),
                  residual_cells=("residual_cells", "mean"), residual_pixels=("residual_pixels", "mean"))
             .reindex(REGION_ORDER).round(3))

    def row(measure, start):
        return tt[(tt.measure == measure) & tt.contrast.str.startswith(start)].iloc[0]

    lines = [
        "# Delivery by brain region",
        "",
        "Generated by `scripts/region_delivery.py`; do not edit by hand. Mice "
        f"{sdd.MICE[0]}-{sdd.MICE[-1]}; {len(f)} single-sonication targets "
        "(Mouse 1's T1, three sonications, and T3, the no-FUS control, are left out).",
        "",
        "Regions are the planned targets from the lab's mask: T1/T4 anterior striatum, T2/T5 "
        "lateral striatum and cortex, T3/T6 thalamus; T1-T3 on the animal's right.",
        "",
        "**Method.** Delivery is put on a logit scale and fitted against log dose with one offset "
        "per mouse and a common slope (cells: "
        f"{slopes['cells']:.2f}, pixels: {slopes['pixels']:.2f}). The residual is delivery above "
        "or below what that mouse's dose-response predicts. p-values shuffle target labels within "
        f"each mouse ({N_PERM:,} permutations), so differences between animals can't drive them.",
        "",
        "## By region",
        "",
        "| Region | Targets | Mean dose | Raw GFP+ neurons | Raw pixel coverage | Residual (cells) | Residual (pixels) |",
        "|---|---|---|---|---|---|---|",
        *[f"| {r} | {int(x.targets)} | {x.mean_dose:.2f} | {x.raw_fraction:.1%} | {x.raw_coverage:.1%} | "
          f"{x.residual_cells:+.2f} | {x.residual_pixels:+.2f} |" for r, x in reg.iterrows()],
        "",
        "## Tests",
        "",
        "| Contrast | Measure | Mean residual difference | Right side lower in | perm p | Per mouse |",
        "|---|---|---|---|---|---|",
        *[f"| {x.contrast} | {x.measure} | {x.mean_residual_diff:+.2f} | "
          f"{'' if pd.isna(x.n_mice) else f'{x.mice_negative}/{x.n_mice}'} | {x.perm_p:.3f} | "
          f"{x.per_mouse} |" for x in tt.itertuples()],
        "",
        "## Dose by position",
        "",
        f"Over mice 2-24 dose doesn't differ by position (Kruskal-Wallis p = {kw_p:.2f}), but in the "
        "imaged mice it does, by chance:",
        "",
        "| Position | Region | Side | Mean dose, mice 2-24 | Mean dose, imaged mice 2-5 |",
        "|---|---|---|---|---|",
        *[f"| T{x.position} | {x.region} | {x.side} | {x.mean_dose_all:.2f} (n {x.n_all}) | "
          f"{x.mean_dose_imaged:.2f} |" for x in pos.itertuples()],
        "",
        "## Reading",
        "",
        f"- **Thalamus is strongest raw, not after dose.** Its targets got the most dose here "
        f"(mean {reg.loc['thalamus', 'mean_dose']:.2f} vs {reg.loc['anterior striatum', 'mean_dose']:.2f} "
        f"and {reg.loc['lateral striatum / cortex', 'mean_dose']:.2f}); dose-adjusted, thalamus vs the "
        f"rest is {row('cells', 'thalamus').mean_residual_diff:+.2f} (p {row('cells', 'thalamus').perm_p:.2f}).",
        f"- **T1 (right anterior striatum) is lower than T4 (left) at matched dose**, in "
        f"{row('cells', 'T1').mice_negative}/{row('cells', 'T1').n_mice} mice on cells "
        f"(p {row('cells', 'T1').perm_p:.3f}) and {row('pixels', 'T1').mice_negative}/"
        f"{row('pixels', 'T1').n_mice} on pixels (p {row('pixels', 'T1').perm_p:.3f}). It is not a "
        "whole-right-side effect: the other two right-left pairs go both ways.",
        "",
        "## Caveats",
        "",
        "- The lab raised both questions after looking at these same animals' images, so agreement "
        "here is not independent. Mice 6+ are the clean test.",
        "- \"Region\" is the planned target. Sections sit ~1.7-2.7 mm deep (`fus.reconstruct`), "
        "above most focal centres (2.9-3.7 mm; thalamic ones deepest), so the tissue in a T3/T6 ROI "
        "may not be thalamus.",
        "- From Mouse 2 on, position = firing order: T1 is always fired first and T4 fourth, so the "
        "T1-T4 gap could be firing order rather than side. Their doses in the imaged mice are "
        f"about equal ({pos.set_index('position').mean_dose_imaged[1]:.2f} vs "
        f"{pos.set_index('position').mean_dose_imaged[4]:.2f} mean).",
        "- 4-5 mice per contrast; three region contrasts and three pairs were tested per measure, "
        "uncorrected.",
        "",
    ]
    (OUT / "README.md").write_text("\n".join(lines))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d, slopes = add_residuals(table())
    tt = tests(d)
    pos, kw_p, rec = dose_by_position()
    cols = ["mouse", "target", "region", "side", "dose", "n_sonications", "in_fit", "gfp_fraction",
            "pixel_coverage", "residual_cells", "residual_pixels"]
    d[cols].round(4).to_csv(OUT / "region_delivery.csv", index=False)
    tt.to_csv(OUT / "region_tests.csv", index=False)
    pos.to_csv(OUT / "dose_by_position.csv", index=False)
    fig = figure(d, tt, rec, kw_p)
    fig.savefig(OUT / "region_delivery.png", dpi=170)
    plt.close(fig)
    write_readme(d, tt, slopes, pos, kw_p)
    print(tt.to_string(index=False))


if __name__ == "__main__":
    main()
