"""The constant treatment (120 bursts at goal 0.75) in every imaged mouse, on one plot.

    .venv/bin/python scripts/reference_target.py

Every mouse got exactly one 120-burst, goal-0.75 sonication: the candidate internal standard
for normalising delivery across animals. This plots its delivery -- the fraction of NeuN cells
GFP+, summed over the mouse's kept slices (sum GFP+ / sum NeuN) -- against its acoustic dose,
with each slice's own value behind it.

Mouse 1 is drawn hollow: its reference shot went into position 1 together with two other
sonications (N. Todd, 2026-09-18), so that site's dose is all three summed.

Inputs: ``results/histology/slices/slices_dose_delivery.csv`` (run
``scripts/slices_dose_delivery.py`` first) and ``results/acoustic_eda/recordings.csv``.
Writes ``results/histology/slices/reference_target.{png,csv}`` (the reference alone) and
``reference_normalised.{png,csv}`` (every target of Mice 2-5 divided by its mouse's reference).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import slices_dose_delivery as sdd  # noqa: E402  (shared palette and styling)

ROOT = sdd.ROOT
OUT = sdd.OUT


def reference_positions() -> dict[int, int]:
    """Position of the 120-burst, goal 0.75-band sonication, per mouse."""
    rec = pd.read_csv(ROOT / "results/acoustic_eda/recordings.csv")
    ref = rec[(rec.programmed == 120) & (rec.goal_lo == 0.7) & (rec.goal_hi == 0.8)]
    return ref.groupby("mouse").position.first().to_dict()


def table() -> tuple[pd.DataFrame, pd.DataFrame]:
    d = pd.read_csv(OUT / "slices_dose_delivery.csv")
    d = d[d.dropped.isna() | (d.dropped == "")]
    pos = reference_positions()
    d = d[[t == pos.get(m) for m, t in zip(d.mouse, d.target)]].copy()
    d["fraction"] = d.n_gfp_pos / d.n_cells
    per_mouse = (d.groupby("mouse")
                   .agg(target=("target", "first"), dose=("dose", "first"),
                        n_sonications=("n_sonications", "first"), n_slices=("section", "nunique"),
                        gfp_cells=("n_gfp_pos", "sum"), cells=("n_cells", "sum"),
                        slice_min=("fraction", "min"), slice_max=("fraction", "max"))
                   .reset_index())
    per_mouse["fraction"] = per_mouse.gfp_cells / per_mouse.cells
    return d, per_mouse


def figure(d: pd.DataFrame, m: pd.DataFrame):
    sdd._style()
    fig, ax = plt.subplots(figsize=(8.6, 6.0))
    fig.subplots_adjust(left=0.1, right=0.97, top=0.8, bottom=0.12)
    sdd._clean(ax)
    rng = np.random.default_rng(0)
    for p in m.itertuples():
        col = sdd.MOUSE_COLOR[p.mouse]
        mixed = p.n_sonications > 1
        s = d[d.mouse == p.mouse]
        ax.scatter(s.dose + rng.uniform(-0.025, 0.025, len(s)), s.fraction, s=28, zorder=2,
                   linewidths=1.2, facecolors=[sdd.SURFACE if lo else col for lo in s.low_signal],
                   edgecolors=col, alpha=0.8)
        ax.scatter(p.dose, p.fraction, s=320, zorder=4, linewidths=2.4,
                   facecolors=sdd.SURFACE if mixed else col, edgecolors=col if mixed else sdd.SURFACE,
                   label=f"Mouse {p.mouse} · T{p.target} · {p.n_slices} slices"
                         + (" · 3 sonications at this site" if mixed else ""))
        ax.annotate(str(p.mouse), (p.dose, p.fraction), ha="center", va="center_baseline",
                    fontsize=10, weight="bold", zorder=5, color=col if mixed else "#ffffff")
        ax.annotate(f"{p.fraction:.1%}", (p.dose, p.fraction), xytext=(15, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=9, color=sdd.INK_2)
    ax.set_xlim(0, max(m.dose.max(), d.dose.max()) * 1.12)
    ax.set_ylim(-0.03, max(0.6, d.fraction.max() * 1.1))
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    ax.set_xlabel("Acoustic dose at the reference site (cumulative 2nd harmonic)", fontsize=11,
                  color=sdd.INK, labelpad=6)
    ax.set_ylabel("Neurons GFP+ in the reference ROI", fontsize=11, color=sdd.INK, labelpad=6)
    ax.tick_params(labelsize=9.5, labelcolor=sdd.INK_2)
    ax.legend(loc="upper left", bbox_to_anchor=(0, 1.17), ncol=2, frameon=False, fontsize=8.5,
              handletextpad=0.4, columnspacing=1.4)
    sdd._header(fig, "The constant treatment (120 bursts, goal 0.75) in each mouse")
    sdd._footer(fig, "Big dot: sum GFP+ / sum NeuN over the mouse's kept slices. Small dots: each slice "
                "(hollow = low-signal slice). Pipeline B, 1.74 mm ROI. Mouse 1 is hollow: its "
                "reference site also got two other sonications.")
    return fig


def normalised() -> pd.DataFrame:
    """Every target's sum GFP+ / sum NeuN divided by its own mouse's reference target.

    Mouse 1 is left out: its reference site also got two other sonications.
    """
    c = pd.read_csv(OUT / "cumulative_dose_delivery.csv")
    c["fraction"] = c.gfp_cells / c.cells
    pos = reference_positions()
    rows = []
    for mouse, g in c[c.mouse != 1].groupby("mouse"):
        ref = g[g.target == pos[mouse]].iloc[0]
        for p in g.itertuples():
            rows.append({"mouse": mouse, "target": p.target, "is_reference": p.target == pos[mouse],
                         "dose": p.dose, "fraction": p.fraction, "ref_fraction": ref.fraction,
                         "ref_dose": ref.dose, "normalised": p.fraction / ref.fraction})
    return pd.DataFrame(rows)


def normalised_figure(n: pd.DataFrame):
    from scipy import stats

    sdd._style()
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 6.0), squeeze=False)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.78, bottom=0.14, wspace=0.22)
    others = n[~n.is_reference]
    panels = [("fraction", "Not normalised: neurons GFP+ (sum GFP+ / sum NeuN)", False),
              ("normalised", "Normalised: ÷ the same mouse's reference target (log scale)", True)]
    for ax, (col, title, log) in zip(axes[0], panels):
        sdd._clean(ax)
        if log:
            ax.set_yscale("log")
            ax.axhline(1, color=sdd.AXIS, linewidth=1, linestyle=":", zorder=1)
        for mouse, g in n.groupby("mouse"):
            colr = sdd.MOUSE_COLOR[mouse]
            ref = g[g.is_reference]
            rest = g[~g.is_reference]
            ax.scatter(rest.dose, rest[col], s=200, color=colr, edgecolors=sdd.SURFACE,
                       linewidths=2, zorder=3,
                       label=f"Mouse {mouse} (reference T{int(ref.target.iloc[0])}: "
                             f"{ref.fraction.iloc[0]:.1%} GFP+)")
            ax.scatter(ref.dose, ref[col], s=260, marker="D", facecolors=sdd.SURFACE,
                       edgecolors=colr, linewidths=2.2, zorder=4)
            for p in g.itertuples():
                ax.annotate(str(p.target), (p.dose, getattr(p, col)), ha="center",
                            va="center_baseline", fontsize=8.5, weight="bold", zorder=5,
                            color=colr if p.is_reference else "#ffffff")
        x, y = others.dose, others[col]
        r = stats.pearsonr(x, np.log(y) if log else y).statistic
        rho = stats.spearmanr(x, y).statistic
        ax.text(0.98, 0.03, f"{'r on log scale' if log else 'r'} = {r:+.2f}, ρ = {rho:+.2f}\n"
                f"n = {len(others)} targets (references left out)",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=sdd.INK_2)
        if not log:
            ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
            ax.set_ylim(-0.03, n.fraction.max() * 1.12)
        else:
            ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}×"))
        ax.set_xlim(0, n.dose.max() * 1.08)
        ax.set_xlabel("Acoustic dose (cumulative 2nd harmonic)", fontsize=10, color=sdd.INK)
        ax.set_title(title, loc="left", fontsize=10.5, weight="bold", color=sdd.INK, pad=6)
    axes[0, 0].set_ylabel("Delivery summed over slices", fontsize=10, color=sdd.INK)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.005, 0.9), ncol=4,
               frameon=False, fontsize=8.5)
    sdd._header(fig, "All targets, normalised by each mouse's constant treatment (120 bursts, goal 0.75)")
    sdd._footer(fig, "Number = target; hollow diamond = that mouse's reference target (1× by construction). "
                "Mice 2-5; Mouse 1 is left out because its reference site also got two other sonications. "
                "Pipeline B, 1.74 mm ROIs.")
    return fig


def main():
    d, m = table()
    m.round(4).to_csv(OUT / "reference_target.csv", index=False)
    fig = figure(d, m)
    fig.savefig(OUT / "reference_target.png", dpi=170)
    plt.close(fig)
    n = normalised()
    n.round(4).to_csv(OUT / "reference_normalised.csv", index=False)
    fig = normalised_figure(n)
    fig.savefig(OUT / "reference_normalised.png", dpi=170)
    plt.close(fig)
    print(n.pivot(index="mouse", columns="target", values="normalised").round(2).to_string())
    print(m.assign(fraction=m.fraction.round(3), slice_min=m.slice_min.round(3),
                   slice_max=m.slice_max.round(3)).to_string(index=False))
    print(f"wrote {OUT.relative_to(ROOT)}/reference_target.png and reference_target.csv")


if __name__ == "__main__":
    main()
