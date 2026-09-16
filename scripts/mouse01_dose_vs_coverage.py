"""Mouse 1: acoustic dose vs GFP coverage, under explicit assumptions.

    PYTHONPATH=src python scripts/mouse01_dose_vs_coverage.py

Reads the per-slot histology table (``python -m fus.histology measure``) and
the seven Mouse_Cntr_01 recordings, joins them per target, and plots coverage
against dose. Writes to ``results/histology/``:

  mouse01_dose_vs_coverage.{csv,png}   both mappings x both dose measures x both channels
  mouse01_dose_delivery_scatter.png    one correlation: measured dose vs TRITC
                                       coverage, deck mapping, least-squares line

This is exploratory. Every join below rests on an assumption that has not been
confirmed with the lab:

A1  Section orientation: s2 and s4 as scanned, s5 and s6 mirrored (the best
    agreement on targets 1, 2, 4, 5). s3 is excluded (plan fit failed).
A2  The empty site is target 3, per the slide deck's "No FUS control".
A3  Coverage per target is the mean over s2, s4, s5, s6; bars span min-max.
A4  Which recordings were fired at which position -- two versions, one per
    row of the figure:
      deck         pos 1 = Target1 + Target2 + Target3 (120 + 60 + 240 = the
                   deck's "420 bursts"); every other position follows its
                   filename, with Target2_Repeat at pos 2.
      spreadsheet  the xlsx `Brain Region` column: pos 1 = Target1 + Target2 +
                   Target6, pos 4 = Target3, pos 5 = Target4, pos 6 = Target5.
A5  Repeat sonications at one position add: doses are summed.
A6  The no-FUS control has zero dose.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
import numpy as np
import pandas as pd
from scipy import stats

from fus import extract
from fus.histology import slot_to_target

ROOT = Path(__file__).resolve().parents[1]
ACOUSTIC = ROOT / "data" / "acoustic" / "20260611"
SLOTS = ROOT / "results" / "histology" / "mouse01_slots.csv"
OUT = ROOT / "results" / "histology" / "mouse01_dose_vs_coverage"
SCATTER = ROOT / "results" / "histology" / "mouse01_dose_delivery_scatter.png"

MIRRORED = {"section_s2": False, "section_s4": False, "section_s5": True, "section_s6": True}

MAPPINGS = {
    "deck": {1: ["Target1", "Target2", "Target3"], 2: ["Target2_Repeat"], 3: [],
             4: ["Target4"], 5: ["Target5"], 6: ["Target6"]},
    "spreadsheet": {1: ["Target1", "Target2", "Target6"], 2: ["Target2_Repeat"], 3: [],
                    4: ["Target3"], 5: ["Target4"], 6: ["Target5"]},
}
ROW_TITLES = {
    "deck": "Mapping A — slide deck",
    "spreadsheet": "Mapping B — spreadsheet Brain Region",
}

# Reference palette (dataviz skill), light mode.
SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
SERIES = {"TRITC": "#2a78d6", "FITC": "#eb6834"}
SERIES_LABEL = {"TRITC": "TRITC (anti-GFP stain)", "FITC": "FITC (native GFP)"}


def section_coverage() -> pd.DataFrame:
    """One row per section x target x channel, with slots relabelled (A1, A2)."""
    d = pd.read_csv(SLOTS)
    d = d[d.section.isin(MIRRORED)].copy()
    d["target"] = [slot_to_target(s, MIRRORED[sec]) for s, sec in zip(d.slot, d.section)]
    return d[["section", "target", "channel", "coverage"]]


def coverage_by_target(sections: pd.DataFrame) -> pd.DataFrame:
    return (
        sections.groupby(["target", "channel"])["coverage"]
        .agg(["mean", "min", "max"])
        .unstack("channel")
    )


def recordings() -> dict[str, dict]:
    out = {}
    for path in sorted(ACOUSTIC.glob("Mouse_Cntr_01_Target*.mat")):
        r = extract(path).to_dict()
        out[path.stem.removeprefix("Mouse_Cntr_01_")] = r
    return out


def join(recs: dict[str, dict], cov: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for mapping, positions in MAPPINGS.items():
        for target, files in positions.items():
            rs = [recs[f] for f in files]
            row = {
                "mapping": mapping,
                "target": target,
                "recordings": " + ".join(files) or "none (no FUS)",
                "n_sonications": len(files),
                "cum_2nd_harmonic": sum(r["cum_2nd_harmonic_xlsx"] for r in rs),
                "bursts_x_goal": sum(r["n_bursts"] * r["harmonic_setpoint"] for r in rs),
                "interlock_fired": any(r["goal_was_reduced"] for r in rs),
            }
            for ch in ("TRITC", "FITC"):
                for stat in ("mean", "min", "max"):
                    row[f"{ch}_{stat}"] = cov.loc[target, (stat, ch)]
            rows.append(row)
    return pd.DataFrame(rows)


def _style() -> None:
    plt.rcParams.update({
        "font.family": ["system-ui", "-apple-system", "Helvetica Neue", "Arial", "sans-serif"],
        "font.size": 10,
        "text.color": INK,
        "axes.labelcolor": INK_2,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
    })


def _clean_axes(ax) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(length=0)


def plot(df: pd.DataFrame):
    _style()
    xs = [
        ("cum_2nd_harmonic", "Measured dose — cumulative 2nd harmonic (sheet units)"),
        ("bursts_x_goal", "Prescribed dose — N bursts × harmonic goal"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharey=True, sharex="col", facecolor=SURFACE)
    for i, mapping in enumerate(MAPPINGS):
        sub = df[df.mapping == mapping].sort_values("target")
        for j, (xcol, xlabel) in enumerate(xs):
            ax = axes[i, j]
            _clean_axes(ax)

            rho = {}
            # Channels share x; FITC is nudged 10 px right on screen so the two
            # series' bars and markers do not overlap.
            for ch, filled, dx, size in (("FITC", False, 10, 55), ("TRITC", True, 0, 190)):
                c = SERIES[ch]
                tr = offset_copy(ax.transData, fig=fig, x=dx, y=0, units="dots")
                x, y = sub[xcol].to_numpy(), sub[f"{ch}_mean"].to_numpy()
                err = np.vstack([y - sub[f"{ch}_min"], sub[f"{ch}_max"] - y])
                ax.errorbar(x, y, yerr=err, fmt="none", ecolor=c, elinewidth=1.4,
                            capsize=0, alpha=0.55, zorder=2, transform=tr)
                ax.scatter(x, y, s=size, zorder=3, linewidths=2, transform=tr,
                           facecolors=c if filled else SURFACE, edgecolors=c if not filled else SURFACE,
                           label=SERIES_LABEL[ch])
                rho[ch] = stats.spearmanr(x, y).statistic

            # Target number inside each TRITC marker; dagger beside an interlocked run
            for r in sub.itertuples():
                xy = (getattr(r, xcol), r.TRITC_mean)
                ax.annotate(str(r.target), xy, ha="center", va="center_baseline",
                            fontsize=8.5, weight="bold", color="#ffffff", zorder=4)
                if r.interlock_fired:
                    ax.annotate("†", xy, xytext=(-7, 6), textcoords="offset points",
                                ha="right", fontsize=10, color=INK_2, zorder=4)

            ax.text(0.02, 0.97,
                    f"Spearman ρ   TRITC {rho['TRITC']:+.2f}   FITC {rho['FITC']:+.2f}   (n = 6)",
                    transform=ax.transAxes, ha="left", va="top", fontsize=9, color=MUTED)
            ax.set_ylim(-0.03, 1.0)
            ax.set_xlim(left=-0.05 * df[xcol].max(), right=1.1 * df[xcol].max())
            ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
            if i == 1:
                ax.set_xlabel(xlabel)
            if j == 0:
                ax.set_ylabel("GFP area coverage in target ROI")
            ax.set_title(ROW_TITLES[mapping] if j == 0 else "", loc="left",
                         fontsize=11, color=INK, weight="bold")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles[::-1], labels[::-1], loc="upper right", ncol=2, frameon=False,
               bbox_to_anchor=(0.99, 0.965), fontsize=10)
    fig.suptitle("Mouse 1 — GFP coverage vs FUS dose", x=0.01, y=0.975, ha="left",
                 fontsize=14, weight="bold")
    fig.text(
        0.01, 0.005,
        "Exploratory — rests on unconfirmed assumptions (see script header). One animal; numbers are targets. "
        "Points: mean of sections s2, s4, s5, s6; bars: min–max.\n"
        "Target 1 received three sonications (doses summed); target 3 is the no-FUS control.  "
        "† = wideband interlock lowered the goal mid-run.  Mappings A and B disagree on which exposure went to targets 1, 4, 5, 6.",
        fontsize=8.5, color=MUTED, va="bottom",
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.94))
    return fig


def plot_scatter(df: pd.DataFrame, mapping: str = "deck",
                 xcol: str = "cum_2nd_harmonic", channel: str = "TRITC"):
    """Single dose-delivery correlation: one dot per target, least-squares line."""
    _style()
    sub = df[df.mapping == mapping].sort_values("target")
    x, y = sub[xcol].to_numpy(), sub[f"{channel}_mean"].to_numpy()

    fig, ax = plt.subplots(figsize=(7, 5.2), facecolor=SURFACE)
    _clean_axes(ax)

    fit = stats.linregress(x, y)
    gx = np.array([0, x.max() * 1.05])
    ax.plot(gx, fit.intercept + fit.slope * gx, color=MUTED, linewidth=1.5, zorder=2)

    ax.scatter(x, y, s=260, color=SERIES[channel], edgecolors=SURFACE, linewidths=2, zorder=3)
    for r in sub.itertuples():
        ax.annotate(str(r.target), (getattr(r, xcol), getattr(r, f"{channel}_mean")),
                    ha="center", va="center_baseline", fontsize=10, weight="bold",
                    color="#ffffff", zorder=4)

    ax.text(0.03, 0.97, f"r = {fit.rvalue:.2f},  n = {len(x)} targets",
            transform=ax.transAxes, ha="left", va="top", fontsize=11, color=INK_2)

    ax.set_xlim(-0.05 * x.max(), x.max() * 1.1)
    ax.set_ylim(-0.03, 1.0)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.set_xlabel("Acoustic dose (cumulative 2nd harmonic)")
    ax.set_ylabel("AAV delivery (GFP coverage)")
    ax.set_title("Mouse 1: acoustic dose vs AAV delivery", loc="left",
                 fontsize=13, weight="bold", pad=12)
    fig.tight_layout()
    return fig


def main() -> None:
    sections = section_coverage()
    df = join(recordings(), coverage_by_target(sections))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.round(4).to_csv(OUT.with_suffix(".csv"), index=False)
    plot(df).savefig(OUT.with_suffix(".png"), dpi=150, facecolor=SURFACE)
    plot_scatter(df).savefig(SCATTER, dpi=150, facecolor=SURFACE)
    print(df[["mapping", "target", "recordings", "cum_2nd_harmonic", "bursts_x_goal",
              "TRITC_mean", "FITC_mean"]].round(3).to_string(index=False))
    print(f"\nwrote {OUT.with_suffix('.png')}\nwrote {SCATTER}")


if __name__ == "__main__":
    main()
