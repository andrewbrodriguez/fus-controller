"""Dose vs delivery: cumulative across slices (headline), per slice, and pooled.

    .venv/bin/python scripts/slices_dose_delivery.py

Runs once per delivery measure, on the same ROIs:

* **pixels** (pipeline A) -- anti-GFP coverage, ``results/histology/mouseNN_template_slots.csv``
* **cells** (pipeline B) -- fraction of NeuN cells tagged GFP+,
  ``results/histology/cells/mouseNN_cell_rois.csv`` (``python -m fus.cells``)

Other inputs: the section sheet (exclusions and review verdicts), and each target's acoustic dose
(cumulative 2nd harmonic, from the recordings; Mouse 1 target 1 is the sum of its
three sonications and target 3 is the no-FUS control at 0).

Slices are dropped if excluded or marked wrong in review. Of the rest, a slice is
**low-signal** when its mean *pixel* coverage over the six targets is below
`LOW_SIGNAL_MEAN` (the same slices for both measures) -- most likely cut outside most targets' focal column. That
cutoff was chosen *after* looking at the per-slice table (2026-09-30), so every
aggregate is reported with and without low-signal slices.

**Cumulative across slices (the headline, per N. Todd, 2026-10-02).** Each target's
delivery is summed over every slice of its animal, giving one number per target that
approximates the transduced volume: total GFP+ cells (B) and total GFP+ area in mm²
(A: pixel coverage × ROI tissue area). Low-signal slices stay in, since a slice that
misses a target's focal column simply contributes little; no post hoc cutoff applies.
Raw totals aren't comparable between animals (Mouse 1 has 4 slices, Mouse 2 has 12, and
staining and section depth differ), so each target's total is divided by its animal's mean
target total: 1.0 = that animal's average target. This cancels the slice count and anything
else that scales a whole animal, and puts both mice on one axis for a pooled fit (12 targets).

Standardising within a slice (z-score over its six targets) removes what a slice
shares -- depth, staining, exposure -- and leaves how targets rank against each
other, which is what dose should explain. Statistics treat slices, not
slice x target points, as the replicates: within each slice, Spearman rho of
coverage vs dose; then how many slices are positive.

Writes to ``results/histology/slices/``:

  cumulative_dose_delivery.csv  one row per mouse x target, delivery summed over slices
  cumulative_pooled.png         both mice on one axis, relative to each animal's mean target
  cumulative.png                raw totals, one panel per mouse and pipeline
  slices_dose_delivery.csv      long table: one row per slice x target, both measures
  heatmap.png, per_slice.png, aggregate_standardized.png          pipeline A (pixels)
  heatmap_cells.png, per_slice_cells.png, aggregate_standardized_cells.png   pipeline B
  pixels_vs_cells.png           B against A for every slice x target
  README.md                     what each figure shows, with the headline numbers
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from fus import extract, orientation

ROOT = Path(__file__).resolve().parents[1]
ACOUSTIC = ROOT / "data/acoustic/20260611"
OUT = ROOT / "results/histology/slices"
DOSES = OUT / "doses.csv"
MICE = (1, 2)

LOW_SIGNAL_MEAN = 0.2  # post hoc -- see module docstring

MEASURES = {
    "pixels": {"suffix": "", "name": "GFP coverage", "pipeline": "A",
               "what": "fraction of the ROI's pixels that are anti-GFP positive"},
    "cells": {"suffix": "_cells", "name": "Fraction of cells GFP+", "pipeline": "B",
              "what": "fraction of NeuN cells in the ROI tagged GFP+ (mean anti-GFP > 4.05x the slice's background)"},
}
M = MEASURES["pixels"]  # the measure being drawn; set in main()

# Recordings per target. Mouse 1 was the operator deviation (docs/mouse01-dose-delivery.md).
RECORDINGS = {
    1: {1: ["Target1", "Target2", "Target3"], 2: ["Target2_Repeat"], 3: [],
        4: ["Target4"], 5: ["Target5"], 6: ["Target6"]},
    2: {t: [f"Target{t}"] for t in range(1, 7)},
}

# Reference palette (dataviz skill), light mode: categorical slots 1-2, blue ramp.
SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
MOUSE_COLOR = {1: "#2a78d6", 2: "#eb6834"}
BLUES = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


def doses() -> pd.DataFrame:
    """Cumulative 2nd harmonic per mouse x target; cached, the .mat files are large."""
    if DOSES.exists():
        return pd.read_csv(DOSES)
    rows = []
    for mouse, targets in RECORDINGS.items():
        for target, files in targets.items():
            recs = [extract(ACOUSTIC / f"Mouse_Cntr_{mouse:02d}_{f}.mat").to_dict() for f in files]
            rows.append({"mouse": mouse, "target": target,
                         "dose": sum(r["cum_2nd_harmonic_xlsx"] for r in recs),
                         "n_sonications": len(files)})
    d = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    d.to_csv(DOSES, index=False)
    return d


def table(measure: str = "pixels") -> pd.DataFrame:
    """One row per slice x target; ``coverage`` holds the chosen measure."""
    sheet = orientation.read_sheet()
    parts = []
    for mouse in MICE:
        d = pd.read_csv(ROOT / f"results/histology/mouse{mouse:02d}_template_slots.csv")
        d = d[d.channel == "TRITC"][["section", "target", "coverage", "placement"]].copy()
        d["mouse"] = mouse
        cells = ROOT / f"results/histology/cells/mouse{mouse:02d}_cell_rois.csv"
        if cells.exists():
            b = pd.read_csv(cells)[["section", "target", "fraction_gfp_pos", "n_cells",
                                    "n_gfp_pos", "roi_tissue_mm2"]]
            d = d.merge(b.rename(columns={"fraction_gfp_pos": "cell_fraction"}),
                        on=["section", "target"], how="left")
        else:
            d["cell_fraction"], d["n_cells"] = np.nan, np.nan
            d["n_gfp_pos"], d["roi_tissue_mm2"] = np.nan, np.nan
        parts.append(d)
    d = pd.concat(parts, ignore_index=True).rename(columns={"coverage": "pixel_coverage"})
    # Low signal is judged on pixel coverage, so both measures use the same slices.
    d["low_signal"] = d.groupby(["mouse", "section"]).pixel_coverage.transform("mean") < LOW_SIGNAL_MEAN
    d["coverage"] = d.pixel_coverage if measure == "pixels" else d.cell_fraction
    d = d[d.groupby(["mouse", "section"]).coverage.transform("count") == 6]  # sections measured

    status = {}
    for (mouse, section), a in sheet.items():
        status[(mouse, section)] = ("excluded" if a.exclude else
                                    "marked wrong" if a.review == "wrong" else "")
    d["dropped"] = [status.get((m, s), "") for m, s in zip(d.mouse, d.section)]
    d = d.merge(doses(), on=["mouse", "target"])

    d["slice_mean"] = d.groupby(["mouse", "section"]).coverage.transform("mean")
    g = d.groupby(["mouse", "section"]).coverage
    d["z"] = (d.coverage - g.transform("mean")) / g.transform("std").replace(0, np.nan)
    d["label"] = [slice_name(m, s) for m, s in zip(d.mouse, d.section)]
    return d.sort_values(["mouse", "section", "target"]).reset_index(drop=True)


def slice_name(mouse: int, section: str) -> str:
    """section_s4 -> "s4"; slide01_s3 -> "slide 1 · s3"."""
    if section.startswith("slide"):
        slide, sec = section.split("_")
        return f"slide {int(slide[5:])} · {sec}"
    return section.removeprefix("section_")


def cumulative(d: pd.DataFrame) -> pd.DataFrame:
    """One row per mouse x target: delivery summed over every kept slice of that animal."""
    use = d[(d.dropped == "") & d.n_gfp_pos.notna()].copy()
    use["gfp_area_mm2"] = use.pixel_coverage * use.roi_tissue_mm2
    c = (use.groupby(["mouse", "target"])
            .agg(dose=("dose", "first"), n_sonications=("n_sonications", "first"),
                 n_slices=("section", "nunique"), gfp_area_mm2=("gfp_area_mm2", "sum"),
                 gfp_cells=("n_gfp_pos", "sum"), cells=("n_cells", "sum"))
            .reset_index())
    for col in CUMULATIVE:
        c[f"{col}_rel"] = c[col] / c.groupby("mouse")[col].transform("mean")
    return c


CUMULATIVE = {"gfp_area_mm2": ("A", "GFP+ area, summed over slices (mm²)"),
              "gfp_cells": ("B", "GFP+ cells, summed over slices")}


def cumulative_fits(c: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col, (pipe, _) in CUMULATIVE.items():
        for mouse in MICE:
            m = c[c.mouse == mouse]
            fit = stats.linregress(m.dose, m[col])
            rows.append({"pipeline": pipe, "measure": col, "mouse": mouse, "n": len(m),
                         "slope": fit.slope, "intercept": fit.intercept, "r": fit.rvalue,
                         "p": fit.pvalue, "rho": stats.spearmanr(m.dose, m[col]).statistic})
    return pd.DataFrame(rows)


def pooled_fits(c: pd.DataFrame) -> pd.DataFrame:
    """Both mice together on the relative scale (each animal's mean target = 1)."""
    rows = []
    for col, (pipe, _) in CUMULATIVE.items():
        for label, m in (("both", c), ("both, no control", c[~((c.mouse == 1) & (c.target == 3))])):
            fit = stats.linregress(m.dose, m[f"{col}_rel"])
            rows.append({"pipeline": pipe, "measure": col, "targets": label, "n": len(m),
                         "slope": fit.slope, "intercept": fit.intercept, "r": fit.rvalue,
                         "p": fit.pvalue, "rho": stats.spearmanr(m.dose, m[f"{col}_rel"]).statistic})
    return pd.DataFrame(rows)


def per_slice_rho(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (mouse, section), g in d.groupby(["mouse", "section"]):
        rho = stats.spearmanr(g.dose, g.coverage).statistic
        rows.append({"mouse": mouse, "section": section, "label": g.label.iloc[0],
                     "rho": rho, "low_signal": g.low_signal.iloc[0]})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def _style():
    plt.rcParams.update({
        "font.family": ["Helvetica Neue", "Arial", "DejaVu Sans"],
        "font.size": 9.5, "text.color": INK, "axes.labelcolor": INK_2,
        "axes.titlecolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.edgecolor": AXIS, "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    })


def _clean(ax, grid_axis="both"):
    ax.set_facecolor(SURFACE)
    ax.grid(True, axis=grid_axis, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(length=0)


def _header(fig, title, subtitle=""):
    fig.text(0.012, 0.985, title, ha="left", va="top", fontsize=14, weight="bold", color=INK)
    if subtitle:
        fig.text(0.012, 0.985 - 0.035 * (8 / fig.get_figheight()), subtitle, ha="left", va="top",
             fontsize=9.5, color=INK_2)


def _footer(fig, text):
    fig.text(0.012, 0.008, text, ha="left", va="bottom", fontsize=8, color=MUTED)


def _dose_label(mouse, target, dose):
    extra = " (×3)" if mouse == 1 and target == 1 else " (no FUS)" if mouse == 1 and target == 3 else ""
    return f"T{target}\n{dose:.2f}{extra}"


def heatmap(d: pd.DataFrame):
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("blues", BLUES)
    n_rows = {m: d[d.mouse == m].section.nunique() for m in MICE}
    fig = plt.figure(figsize=(8.2, 2.2 + 0.34 * sum(n_rows.values()) + 0.9))
    gs = fig.add_gridspec(len(MICE), 2, height_ratios=[n_rows[m] for m in MICE],
                          width_ratios=[40, 1], hspace=0.75, wspace=0.04,
                          left=0.24, right=0.92, top=0.84, bottom=0.07)
    for row, mouse in enumerate(MICE):
        ax = fig.add_subplot(gs[row, 0])
        m = d[d.mouse == mouse]
        dose = m.groupby("target").dose.first().sort_values()
        order = dose.index.tolist()
        sections = sorted(m.section.unique())
        grid = m.pivot(index="section", columns="target", values="coverage").loc[sections, order]
        flags = m.groupby("section").agg(low=("low_signal", "first"), drop=("dropped", "first"))
        im = ax.imshow(grid.to_numpy(), cmap=cmap, vmin=0, vmax=1, aspect="auto")
        for (i, j), v in np.ndenumerate(grid.to_numpy()):
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8.5,
                    color="#ffffff" if v > 0.5 else INK)
        ax.set_xticks(np.arange(-0.5, len(order)), minor=True)
        ax.set_yticks(np.arange(-0.5, len(sections)), minor=True)
        ax.grid(which="minor", color=SURFACE, linewidth=2)
        # a wider gap between slides
        slides = [s.split("_")[0] for s in sections]
        for i in range(1, len(slides)):
            if slides[i] != slides[i - 1]:
                ax.axhline(i - 0.5, color=SURFACE, linewidth=6)
        ax.tick_params(which="both", length=0)
        ax.xaxis.tick_top()
        ax.set_xticks(range(len(order)), [_dose_label(mouse, t, dose[t]) for t in order],
                      fontsize=8.5, color=INK_2)
        labels = []
        for sec in sections:
            note = "  dropped" if flags.loc[sec, "drop"] else "  low signal" if flags.loc[sec, "low"] else ""
            labels.append(slice_name(mouse, sec) + note)
        ax.set_yticks(range(len(sections)), labels)
        for sec, tick in zip(sections, ax.get_yticklabels()):
            tick.set_color(MUTED if flags.loc[sec, "low"] or flags.loc[sec, "drop"] else INK)
        for side in ax.spines.values():
            side.set_visible(False)
        ax.set_title(f"Mouse {mouse}", loc="left", fontsize=11, weight="bold", color=INK, pad=34,
                     x=-0.27)
        if row == 0:
            cax = fig.add_subplot(gs[:, 1])
            cb = fig.colorbar(im, cax=cax)
            cb.outline.set_visible(False)
            cb.ax.tick_params(length=0, labelsize=8, colors=MUTED)
            cb.ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
            cb.set_label(M["name"], color=INK_2)
    _header(fig, f"{M['name']} per slice and target (pipeline {M['pipeline']})",
            f"{M['what'][0].upper() + M['what'][1:]}, 1.74 mm ROIs.\nColumns: targets ordered by "
            "acoustic dose (cumulative 2nd harmonic), low → high.")
    _footer(fig, f"Low signal: slice mean pixel coverage < {LOW_SIGNAL_MEAN} (post hoc). "
            "Mouse 1 T1 got three sonications (doses summed).")
    return fig


def per_slice(d: pd.DataFrame, rho: pd.DataFrame):
    rho = rho.set_index("section")
    rows = [("Mouse 1", sorted(d[d.mouse == 1].section.unique()))]
    m2 = sorted(d[d.mouse == 2].section.unique())
    for slide in sorted({s.split("_")[0] for s in m2}):
        rows.append((f"Mouse 2 · slide {int(slide[5:])}", [s for s in m2 if s.startswith(slide)]))
    ncol = max(len(r[1]) for r in rows)
    fig, axes = plt.subplots(len(rows), ncol, figsize=(2.75 * ncol + 1.2, 2.3 * len(rows) + 1.3),
                             sharex=True, sharey=True, squeeze=False)
    fig.subplots_adjust(left=0.11, right=0.985, top=0.88, bottom=0.08, hspace=0.42, wspace=0.12)
    xmax = d.dose.max() * 1.1
    for r, (row_name, sections) in enumerate(rows):
        axes[r, 0].annotate(row_name.replace(" · ", "\n"), xy=(-0.3, 0.5),
                            xycoords="axes fraction", ha="center", va="center", rotation=90,
                            fontsize=9.5, weight="bold", color=INK, linespacing=1.3)
        for c in range(ncol):
            ax = axes[r, c]
            if c >= len(sections):
                ax.set_visible(False)
                continue
            sec = sections[c]
            g = d[d.section == sec].sort_values("dose")
            mouse = int(g.mouse.iloc[0])
            faded = bool(g.low_signal.iloc[0]) or bool(g.dropped.iloc[0])
            col = MOUSE_COLOR[mouse]
            _clean(ax)
            fit = stats.linregress(g.dose, g.coverage)
            gx = np.array([0, g.dose.max()])
            ax.plot(gx, fit.intercept + fit.slope * gx, color=col, linewidth=1.5,
                    alpha=0.25 if faded else 0.5, zorder=2)
            ax.scatter(g.dose, g.coverage, s=120, zorder=3, linewidths=1.6,
                       facecolors=SURFACE if faded else col, edgecolors=col if faded else SURFACE)
            for p in g.itertuples():
                ax.annotate(str(p.target), (p.dose, p.coverage), ha="center",
                            va="center_baseline", fontsize=7.5, weight="bold", zorder=4,
                            color=col if faded else "#ffffff")
            ax.set_title(sec.split("_")[1] if mouse == 2 else sec.removeprefix("section_"),
                         loc="left", fontsize=9.5, color=MUTED if faded else INK, pad=4)
            ax.text(1.0, 1.02, f"ρ {rho.loc[sec, 'rho']:+.2f}" + ("  · low signal" if faded else ""),
                    transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=MUTED)
            ax.set_xlim(-0.1, xmax)
            ax.set_xticks(np.arange(0, xmax, 1.0))
            ax.set_ylim(-0.06, 1.04)
            ax.set_yticks([0, 0.5, 1])
            ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    for c in range(ncol):
        bottom = next(axes[r, c] for r in range(len(rows) - 1, -1, -1) if axes[r, c].get_visible())
        bottom.set_xlabel("dose (cum. 2nd harmonic)")
        bottom.xaxis.set_tick_params(labelbottom=True)
    _header(fig, f"{M['name']} vs dose, one panel per slice (pipeline {M['pipeline']})",
            "Number = target. Line: least squares within the slice. ρ: Spearman with dose. "
            "Hollow = low-signal slice.")
    _footer(fig, f"Y: {M['what']}. Mouse 1 T3 is the no-FUS control at dose 0; "
            "T1 = three sonications summed.")
    return fig


def aggregate(d: pd.DataFrame, rho: pd.DataFrame):
    use = d[~d.low_signal & (d.dropped == "")]
    r_ok = rho[~rho.low_signal]
    pos = int((r_ok.rho > 0).sum())
    fig, ax = plt.subplots(figsize=(8.4, 6.0))
    fig.subplots_adjust(left=0.11, right=0.97, top=0.84, bottom=0.12)
    _clean(ax)
    ax.axhline(0, color=AXIS, linewidth=1, zorder=1)
    rng = np.random.default_rng(0)
    for mouse in MICE:
        m = use[use.mouse == mouse]
        col = MOUSE_COLOR[mouse]
        ax.scatter(m.dose + rng.uniform(-0.018, 0.018, len(m)), m.z, s=16, color=col,
                   alpha=0.35, linewidths=0, zorder=2)
        means = m.groupby("target").agg(dose=("dose", "first"), z=("z", "mean"), se=("z", "sem"))
        ax.errorbar(means.dose, means.z, yerr=means.se, fmt="none", ecolor=col, elinewidth=1.8,
                    capsize=0, zorder=3)
        ax.scatter(means.dose, means.z, s=260, color=col, edgecolors=SURFACE, linewidths=2,
                   zorder=4, label=f"Mouse {mouse} · target mean ± SE · {m.section.nunique()} slices")
        for t, p in means.iterrows():
            ax.annotate(str(t), (p.dose, p.z), ha="center", va="center_baseline", fontsize=9.5,
                        weight="bold", color="#ffffff", zorder=5)
    ax.set_xlim(-0.12, use.dose.max() * 1.08)
    ax.set_xticks(np.arange(0, use.dose.max() * 1.08, 0.5))
    lo, hi = np.floor(use.z.min() * 2) / 2, np.ceil(use.z.max() * 2) / 2
    ax.set_ylim(lo - 0.1, hi + 0.1)
    ax.set_yticks(np.arange(lo, hi + 0.01, 0.5))
    ax.set_xlabel("Acoustic dose (cumulative 2nd harmonic)", fontsize=12.5, color=INK, labelpad=8)
    ax.set_ylabel(f"{M['name']}, z-scored within slice", fontsize=12.5, color=INK, labelpad=8)
    ax.tick_params(labelsize=10, labelcolor=INK_2)
    ax.legend(loc="upper left", bbox_to_anchor=(0, 1.1), ncol=2, frameon=False, fontsize=8.5,
              handletextpad=0.3, columnspacing=1.6)
    _header(fig, f"Dose vs delivery, standardized within each slice (pipeline {M['pipeline']})")
    return fig


def cumulative_figure(c: pd.DataFrame, fits: pd.DataFrame):
    """Rows: pipeline A, B. Columns: mouse. One dot per target, OLS trend line per panel."""
    fig, axes = plt.subplots(2, len(MICE), figsize=(9.6, 7.6), squeeze=False)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.86, bottom=0.1, hspace=0.42, wspace=0.22)
    xmax = c.dose.max() * 1.1
    for r, (col, (pipe, ylabel)) in enumerate(CUMULATIVE.items()):
        for k, mouse in enumerate(MICE):
            ax = axes[r, k]
            _clean(ax)
            m = c[c.mouse == mouse].sort_values("dose")
            f = fits[(fits.measure == col) & (fits.mouse == mouse)].iloc[0]
            colr = MOUSE_COLOR[mouse]
            gx = np.array([0, xmax])
            ax.plot(gx, f.intercept + f.slope * gx, color=colr, linewidth=2, alpha=0.55, zorder=2)
            ax.scatter(m.dose, m[col], s=200, color=colr, edgecolors=SURFACE, linewidths=2, zorder=3)
            for p in m.itertuples():
                ax.annotate(str(p.target), (p.dose, getattr(p, col)), ha="center",
                            va="center_baseline", fontsize=8.5, weight="bold", color="#ffffff", zorder=4)
            top = max(m[col].max(), (f.intercept + f.slope * xmax)) * 1.12
            ax.set_xlim(-0.1, xmax)
            ax.set_ylim(min(0, f.intercept) - 0.04 * top, top)
            ax.set_xticks(np.arange(0, xmax, 0.5))
            ax.set_title(f"Mouse {mouse} · {int(m.n_slices.max())} slices · pipeline {pipe}",
                         loc="left", fontsize=10, weight="bold", color=INK, pad=6)
            ax.text(0.02, 0.97, f"r = {f.r:.2f}  (p = {f.p:.2g})\nρ = {f.rho:+.2f}\n"
                    f"slope = {f.slope:,.3g} per unit dose", transform=ax.transAxes,
                    ha="left", va="top", fontsize=8.5, color=INK_2)
            if k == 0:
                ax.set_ylabel(ylabel, fontsize=10, color=INK, labelpad=6)
            if r == 1:
                ax.set_xlabel("Acoustic dose (cumulative 2nd harmonic)", fontsize=10, color=INK)
    _header(fig, "Cumulative delivery across slices vs acoustic dose",
            "Each dot is a target (number), delivery summed over every slice of that animal. "
            "Line: least squares, n = 6 per panel.")
    _footer(fig, "Mouse 1 T3 is the no-FUS control at dose 0; Mouse 1 T1 = its three sonications at position 1 summed (Target1-3 recordings, per N. Todd 9/18; not the sheet). "
            "Low-signal slices are included. Animals have different slice counts, so compare slopes "
            "within a panel, not heights across.")
    return fig


def pooled_figure(c: pd.DataFrame, pf: pd.DataFrame):
    """Both mice on one axis per pipeline: delivery relative to the animal's mean target."""
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.6), squeeze=False)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.8, bottom=0.14, wspace=0.18)
    xmax = c.dose.max() * 1.08
    ymax = max(c[f"{col}_rel"].max() for col in CUMULATIVE) * 1.12
    for k, (col, (pipe, ylabel)) in enumerate(CUMULATIVE.items()):
        ax = axes[0, k]
        _clean(ax)
        ax.axhline(1, color=AXIS, linewidth=1, linestyle=":", zorder=1)
        f = pf[(pf.measure == col) & (pf.targets == "both")].iloc[0]
        fn = pf[(pf.measure == col) & (pf.targets == "both, no control")].iloc[0]
        gx = np.array([0, xmax])
        ax.plot(gx, f.intercept + f.slope * gx, color=INK_2, linewidth=2, alpha=0.7, zorder=2,
                label="least squares, both mice")
        for mouse in MICE:
            m = c[c.mouse == mouse]
            colr = MOUSE_COLOR[mouse]
            ax.scatter(m.dose, m[f"{col}_rel"], s=200, color=colr, edgecolors=SURFACE,
                       linewidths=2, zorder=3, label=f"Mouse {mouse} ({int(m.n_slices.max())} slices)")
            for p in m.itertuples():
                ax.annotate(str(p.target), (p.dose, getattr(p, f"{col}_rel")), ha="center",
                            va="center_baseline", fontsize=8.5, weight="bold", color="#ffffff", zorder=4)
        ax.set_xlim(-0.1, xmax)
        ax.set_ylim(-0.08, ymax)
        ax.set_xticks(np.arange(0, xmax, 0.5))
        ax.set_title(f"Pipeline {pipe}: {ylabel.split(',')[0]}", loc="left", fontsize=10.5,
                     weight="bold", color=INK, pad=6)
        ax.text(0.98, 0.03, f"r = {f.r:.2f} (p = {f.p:.2g}), ρ = {f.rho:+.2f}, n = {f.n}\n"
                f"without Mouse 1 control: r = {fn.r:.2f} (p = {fn.p:.2g}), n = {fn.n}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=INK_2)
        ax.set_xlabel("Acoustic dose (cumulative 2nd harmonic)", fontsize=10, color=INK)
        if k == 0:
            ax.set_ylabel("Delivery summed over slices,\nrelative to the animal's mean target",
                          fontsize=10, color=INK, labelpad=6)
            ax.legend(loc="upper left", frameon=False, fontsize=8.5)
    _header(fig, "Cumulative delivery vs acoustic dose, both mice",
            "Each dot is a target (number): delivery summed over every slice of its animal, divided "
            "by that animal's mean target.\nDotted line: the animal's average target (1.0).")
    _footer(fig, "Mouse 1 T3 is the no-FUS control at dose 0; Mouse 1 T1 = its three sonications at position 1 summed (Target1-3 recordings, per N. Todd 9/18; not the sheet). "
            "Low-signal slices included.")
    return fig


def pixels_vs_cells(d: pd.DataFrame):
    """Pipeline B against pipeline A for every slice x target."""
    use = d[(d.dropped == "") & d.cell_fraction.notna()]
    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    fig.subplots_adjust(left=0.13, right=0.97, top=0.86, bottom=0.11)
    _clean(ax)
    ax.plot([0, 1], [0, 1], color=AXIS, linewidth=1.2, linestyle=":", zorder=1)
    for mouse in MICE:
        m = use[use.mouse == mouse]
        col = MOUSE_COLOR[mouse]
        ax.scatter(m.pixel_coverage, m.cell_fraction, s=34, zorder=3, linewidths=1.2,
                   facecolors=[SURFACE if lo else col for lo in m.low_signal],
                   edgecolors=[col if lo else SURFACE for lo in m.low_signal],
                   label=f"Mouse {mouse} ({m.section.nunique()} slices)")
    r = stats.pearsonr(use.pixel_coverage, use.cell_fraction).statistic
    ax.set(xlim=(-0.02, 1.02), ylim=(-0.02, 1.02))
    ax.set_xticks(np.arange(0, 1.01, 0.25))
    ax.set_yticks(np.arange(0, 1.01, 0.25))
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    ax.set_xlabel("A: GFP+ pixel coverage", fontsize=12, color=INK, labelpad=6)
    ax.set_ylabel("B: fraction of cells GFP+", fontsize=12, color=INK, labelpad=6)
    ax.tick_params(labelsize=10, labelcolor=INK_2)
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    _header(fig, "Pipeline B vs pipeline A, every slice × target",
            f"Pearson r = {r:.2f}, n = {len(use)}. Dotted: equal. Hollow = low-signal slice.")
    return fig


def write_readme(results: dict, c: pd.DataFrame, fits: pd.DataFrame, pf: pd.DataFrame) -> None:
    def line(sub):
        pos = int((sub.rho > 0).sum())
        return f"median ρ {sub.rho.median():+.2f}, positive in {pos}/{len(sub)}"

    rows, means_line = [], []
    for label, pick in (("Both mice, excluding low-signal slices", lambda r: r[~r.low_signal]),
                        ("Mouse 1", lambda r: r[~r.low_signal & (r.mouse == 1)]),
                        ("Mouse 2, excluding low-signal slices", lambda r: r[~r.low_signal & (r.mouse == 2)]),
                        ("Mouse 2, all slices", lambda r: r[r.mouse == 2])):
        rows.append(f"| {label} | " + " | ".join(line(pick(rho)) for _, rho in results.values()) + " |")
    for measure, (d, _) in results.items():
        use = d[~d.low_signal & (d.dropped == "")]
        means = use.groupby(["mouse", "target"]).agg(dose=("dose", "first"), z=("z", "mean"))
        by = {m: stats.spearmanr(means.loc[m].dose, means.loc[m].z).statistic for m in MICE}
        means_line.append(f"{MEASURES[measure]['pipeline']} ({measure}): Mouse 1 {by[1]:+.2f}, "
                          f"Mouse 2 {by[2]:+.2f}")
    d = results["pixels"][0]
    both = d[(d.dropped == "") & d.cell_fraction.notna()]
    r_ab = stats.pearsonr(both.pixel_coverage, both.cell_fraction).statistic
    n_low = int(d.groupby(["mouse", "section"]).low_signal.first().sum())
    dropped = sorted(set(d.loc[d.dropped != "", "section"]))
    header = " | ".join(f"Pipeline {MEASURES[m]['pipeline']}: {m}" for m in results)

    fit_rows = "\n".join(
        f"| {f.pipeline} | {f.measure} | Mouse {f.mouse} | {f.slope:,.3g} | {f.r:+.2f} | {f.p:.2g} | {f.rho:+.2f} |"
        for f in fits.itertuples())
    pooled_rows = "\n".join(
        f"| {f.pipeline} | {f.measure} | {f.targets} | {f.n} | {f.slope:.3g} | {f.r:+.2f} | {f.p:.2g} | {f.rho:+.2f} |"
        for f in pf.itertuples())
    m1 = c[(c.mouse == 1) & (c.target != 3)]
    no_ctrl = "; ".join(
        f"{CUMULATIVE[col][0]} r = {stats.pearsonr(m1.dose, m1[col]).statistic:+.2f}"
        for col in CUMULATIVE)

    text = f"""# Dose vs delivery: cumulative across slices, per slice, and pooled

Generated by `scripts/slices_dose_delivery.py`; re-run it after re-measuring or
re-reviewing. Dose = cumulative 2nd harmonic from the recordings. Delivery is
measured two ways in the same hand-finetuned 1.74 mm ROIs:

- **Pipeline A, pixels:** {MEASURES["pixels"]["what"]}.
- **Pipeline B, cells:** {MEASURES["cells"]["what"]}.

## Headline: cumulative delivery across slices

Following N. Todd (2026-10-02), each target's delivery is **summed over every slice** of its
animal instead of treating each slice as its own measurement: total GFP+ cells (pipeline B)
and total GFP+ area (pipeline A, coverage × ROI tissue area). All slices that aren't excluded
or marked wrong are summed, including the low-signal ones, so no post hoc cutoff enters.

Raw totals aren't comparable between animals (Mouse 1 has 4 slices and Mouse 2 has 12, and
staining and section depth differ), so each total is divided by **its animal's mean target
total**: 1.0 = that animal's average target. That cancels the slice count and anything else
that scales a whole animal, and puts both mice on one axis. Figure: **`cumulative_pooled.png`**;
table: `cumulative_dose_delivery.csv` (`*_rel` columns).

| Pipeline | Measure | Targets | n | Slope (relative units per unit dose) | Pearson r | p | Spearman ρ |
|---|---|---|---|---|---|---|---|
{pooled_rows}

The 12 points come from two animals, and the relative scale forces each animal's mean to 1,
so p-values are optimistic: they treat targets as independent.

Per animal, on raw totals (`cumulative.png`, n = 6 each):

| Pipeline | Measure | Animal | Slope per unit dose | Pearson r | p | Spearman ρ |
|---|---|---|---|---|---|---|
{fit_rows}

Mouse 1 without its no-FUS control (T3), 5 targets: {no_ctrl}.

## Per-slice figures (secondary)

Each pipeline has the same three figures; B's carry a `_cells` suffix.

**`heatmap.png`: the raw numbers.** One row per slice, one column per target,
with targets ordered by dose left to right. If delivery followed dose, rows
would get darker toward the right. Grey row labels are low-signal slices.

**`per_slice.png`: delivery vs dose inside each slice.** One panel per slice,
grouped by animal and slide. Each dot is a target, and the line is a
least-squares fit within that slice. ρ is the Spearman rank correlation with
dose in that slice (6 targets).

**`aggregate_standardized.png`: all slices pooled.** Each slice's six values
are z-scored against that slice's own mean and SD. This removes what the whole
slice shares (depth relative to the focal column, staining, exposure) and keeps
how its targets rank against each other. Small dots are slice × target; big dots
are each target's mean ± SE across slices, placed at that target's dose.

**`pixels_vs_cells.png`: the two pipelines against each other**, one dot per
slice × target. Pearson r = {r_ab:.2f} over {len(both)} ROIs.

## Headline numbers

| Within-slice ρ with dose | {header} |
|---|---|---|
{chr(10).join(rows)}

Across the six target means (z), ρ with dose: {"; ".join(means_line)}.

## Read with care

- **Low signal** means a slice mean *pixel* coverage below {LOW_SIGNAL_MEAN}, a cutoff chosen
  after looking at the data. It flags {n_low} slice(s), probably cut outside most targets'
  focal column, and the same slices are set aside for both pipelines.
  {"Dropped (excluded or marked wrong): " + ", ".join(dropped) + "." if dropped else "No slice is excluded or marked wrong in `data/section_orientation.csv`."}
- **Pipeline B's threshold** (4.05× the slice's background, its median anti-GFP over tissue) was
  set on one 3 mm crop of one section and is not yet checked against a hand count
  (`docs/gfp-cell-tagging.md`). Dividing by each slice's background is what keeps Mouse 1's
  no-FUS control near zero: its tissue is several times brighter than Mouse 2's.
- **Mouse 1's correlation leans on the no-FUS control** (T3, dose 0) and on T1, whose
  dose is three sonications summed.
- **Mouse 1 T1's dose (2.44) is the actual dose at position 1**: the `Target1`, `Target2` and
  `Target3` recordings (0.9462 + 0.5658 + 0.9293), all fired there per N. Todd (2026-09-18).
  It comes from the recordings, not the summary sheet, whose Mouse 1 rows 3–6 are rotated
  by one run (the sheet would give 2.83). Mouse 1 T2 is `Target2_Repeat`; T4–T6 are
  `Target4`–`Target6`.
- **Mouse 2 T2** (mid dose, low delivery) is the main exception. It's a right-side
  lateral target, and the head was rolled with the right side higher, so most
  slices probably miss its focal column.
- **Slices are not independent**: they come from two animals. The sign-test p-values
  on the figures describe consistency within these animals, not a population effect.
- **Doses aren't calibrated between animals**; the z-scores compare targets within a
  slice, not absolute delivery across mice.
"""
    (OUT / "README.md").write_text(text)


def main():
    global M
    _style()
    OUT.mkdir(parents=True, exist_ok=True)
    results = {}
    for measure, info in MEASURES.items():
        M = info
        d = table(measure)
        if d.empty:
            print(f"{measure}: no measurements yet, skipped")
            continue
        rho = per_slice_rho(d[d.dropped == ""])
        results[measure] = (d, rho)
        for name, fig in (("heatmap", heatmap(d)), ("per_slice", per_slice(d, rho)),
                          ("aggregate_standardized", aggregate(d, rho))):
            fig.savefig(OUT / f"{name}{info['suffix']}.png", dpi=170)
            plt.close(fig)
    full = table("pixels")
    c = cumulative(full)
    fits = cumulative_fits(c)
    c.round(4).to_csv(OUT / "cumulative_dose_delivery.csv", index=False)
    pf = pooled_fits(c)
    fig = cumulative_figure(c, fits)
    fig.savefig(OUT / "cumulative.png", dpi=170)
    plt.close(fig)
    fig = pooled_figure(c, pf)
    fig.savefig(OUT / "cumulative_pooled.png", dpi=170)
    plt.close(fig)
    print(pf.round(3).to_string(index=False))
    full.drop(columns=["coverage", "slice_mean", "z"]).round(4).to_csv(OUT / "slices_dose_delivery.csv", index=False)
    if "cells" in results:
        fig = pixels_vs_cells(results["pixels"][0])
        fig.savefig(OUT / "pixels_vs_cells.png", dpi=170)
        plt.close(fig)
        write_readme(results, c, fits, pf)
    print(f"wrote {', '.join(p.name for p in sorted(OUT.glob('*')))} to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
