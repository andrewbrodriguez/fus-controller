"""Dose vs delivery per slice, and pooled across slices and mice.

    .venv/bin/python scripts/slices_dose_delivery.py

Inputs: the finetuned/template measurements for each mouse
(``results/histology/mouseNN_template_slots.csv``, TRITC channel), the section
sheet (exclusions and review verdicts), and each target's acoustic dose
(cumulative 2nd harmonic, from the recordings; Mouse 1 target 1 is the sum of its
three sonications and target 3 is the no-FUS control at 0).

Slices are dropped if excluded or marked wrong in review. Of the rest, a slice is
**low-signal** when its mean coverage over the six targets is below
`LOW_SIGNAL_MEAN` -- most likely cut outside most targets' focal column. That
cutoff was chosen *after* looking at the per-slice table (2026-09-30), so every
aggregate is reported with and without low-signal slices.

Standardising within a slice (z-score over its six targets) removes what a slice
shares -- depth, staining, exposure -- and leaves how targets rank against each
other, which is what dose should explain. Statistics treat slices, not
slice x target points, as the replicates: within each slice, Spearman rho of
coverage vs dose; then how many slices are positive.

Writes to ``results/histology/slices/``:

  slices_dose_delivery.csv      long table: one row per slice x target
  heatmap.png                   coverage, slice x target, targets ordered by dose
  per_slice.png                 coverage vs dose, one panel per slice
  aggregate_standardized.png    within-slice z vs dose, pooled; per-target means
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


def table() -> pd.DataFrame:
    sheet = orientation.read_sheet()
    parts = []
    for mouse in MICE:
        d = pd.read_csv(ROOT / f"results/histology/mouse{mouse:02d}_template_slots.csv")
        d = d[d.channel == "TRITC"][["section", "target", "coverage", "placement"]].copy()
        d["mouse"] = mouse
        parts.append(d)
    d = pd.concat(parts, ignore_index=True)

    status = {}
    for (mouse, section), a in sheet.items():
        status[(mouse, section)] = ("excluded" if a.exclude else
                                    "marked wrong" if a.review == "wrong" else "")
    d["dropped"] = [status.get((m, s), "") for m, s in zip(d.mouse, d.section)]
    d = d.merge(doses(), on=["mouse", "target"])

    d["slice_mean"] = d.groupby(["mouse", "section"]).coverage.transform("mean")
    d["low_signal"] = d.slice_mean < LOW_SIGNAL_MEAN
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
            cb.set_label("GFP coverage", color=INK_2)
    _header(fig, "GFP coverage per slice and target",
            "Anti-GFP stain (TRITC) in 1.74 mm ROIs. Columns: targets ordered by acoustic dose "
            "(cumulative 2nd harmonic), low → high.")
    _footer(fig, f"Low signal: slice mean coverage < {LOW_SIGNAL_MEAN} (post hoc). "
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
    _header(fig, "Coverage vs dose, one panel per slice",
            "Number = target. Line: least squares within the slice. ρ: Spearman, coverage vs dose. "
            "Hollow = low-signal slice.")
    _footer(fig, "Y: anti-GFP coverage in the target ROI. Mouse 1 T3 is the no-FUS control at dose 0; "
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
    ax.set_ylabel("Coverage, z-scored within slice", fontsize=12.5, color=INK, labelpad=8)
    ax.tick_params(labelsize=10, labelcolor=INK_2)
    ax.legend(loc="upper left", bbox_to_anchor=(0, 1.1), ncol=2, frameon=False, fontsize=8.5,
              handletextpad=0.3, columnspacing=1.6)
    _header(fig, "Dose vs delivery, standardized within each slice")
    return fig


def write_readme(d: pd.DataFrame, rho: pd.DataFrame) -> None:
    use = d[~d.low_signal & (d.dropped == "")]
    r_ok = rho[~rho.low_signal]
    means = use.groupby(["mouse", "target"]).agg(dose=("dose", "first"), z=("z", "mean"))
    by_mouse = {m: stats.spearmanr(means.loc[m].dose, means.loc[m].z).statistic for m in MICE}
    n_low = int(rho.low_signal.sum())
    dropped = sorted(set(d.loc[d.dropped != "", "section"]))

    def line(sub):
        pos = int((sub.rho > 0).sum())
        return f"median ρ {sub.rho.median():+.2f}, positive in {pos}/{len(sub)}"

    text = f"""# Dose vs delivery, per slice and pooled

Generated by `scripts/slices_dose_delivery.py`; re-run it after re-measuring or
re-reviewing. Delivery = anti-GFP (TRITC) coverage in each target's 1.74 mm ROI,
finetuned by hand. Dose = cumulative 2nd harmonic from the recordings.

## Figures

**`heatmap.png`: the raw numbers.** One row per slice, one column per target,
with targets ordered by dose left to right. If delivery followed dose, rows
would get darker toward the right. Grey row labels are low-signal slices.

**`per_slice.png`: coverage vs dose inside each slice.** One panel per slice,
grouped by animal and slide. Each dot is a target, and the line is a
least-squares fit within that slice. ρ is the Spearman rank correlation of
coverage with dose in that slice (6 targets).

**`aggregate_standardized.png`: all slices pooled.** Each slice's six coverages
are z-scored against that slice's own mean and SD. This removes what the whole
slice shares (depth relative to the focal column, staining, exposure) and keeps
how its targets rank against each other. Small dots are slice × target; big dots
are each target's mean ± SE across slices, placed at that target's dose.

## Headline numbers

| | Within-slice ρ (coverage vs dose) |
|---|---|
| Both mice, excluding low-signal slices | {line(r_ok)} |
| Mouse 1 | {line(r_ok[r_ok.mouse == 1])} |
| Mouse 2, excluding low-signal slices | {line(r_ok[r_ok.mouse == 2])} |
| Mouse 2, all slices | {line(rho[rho.mouse == 2])} |

Across the six target means (z), ρ with dose is {by_mouse[1]:+.2f} for Mouse 1 and
{by_mouse[2]:+.2f} for Mouse 2.

## Read with care

- **Low signal** means a slice mean coverage below {LOW_SIGNAL_MEAN}, a cutoff chosen after
  looking at the data. It flags {n_low} slice(s), probably cut outside most targets'
  focal column. {"Dropped (excluded or marked wrong): " + ", ".join(dropped) + "." if dropped else "No slice is excluded or marked wrong in `data/section_orientation.csv`."}
- **Mouse 1's correlation leans on the no-FUS control** (T3, dose 0) and on T1, whose
  dose is three sonications summed.
- **Mouse 2 T2** (mid dose, low coverage) is the main exception. It's a right-side
  lateral target, and the head was rolled with the right side higher, so most
  slices probably miss its focal column.
- **Slices are not independent**: they come from two animals. The sign-test p-values
  on the figure describe consistency within these animals, not a population effect.
- **Doses aren't calibrated between animals**; the z-scores compare targets within a
  slice, not absolute coverage across mice.
"""
    (OUT / "README.md").write_text(text)


def main():
    _style()
    d = table()
    rho = per_slice_rho(d[d.dropped == ""])
    OUT.mkdir(parents=True, exist_ok=True)
    d.round(4).to_csv(OUT / "slices_dose_delivery.csv", index=False)
    for name, fig in (("heatmap", heatmap(d)), ("per_slice", per_slice(d, rho)),
                      ("aggregate_standardized", aggregate(d, rho))):
        fig.savefig(OUT / f"{name}.png", dpi=170)
        plt.close(fig)
    write_readme(d, rho)
    print(f"wrote {', '.join(p.name for p in sorted(OUT.glob('*')))} to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
