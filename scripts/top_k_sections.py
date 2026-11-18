"""Dose vs delivery using only each mouse's top-k sections, by the hand quality ranking.

    .venv/bin/python scripts/top_k_sections.py

For k = 1, 2, ... each mouse contributes only its k best kept sections, ranked in
``data/section_quality.csv`` (``scripts/rank_sections.py``). Delivery is then summed over those
sections exactly as in ``slices_dose_delivery.py`` (``cumulative``), and dose vs delivery is
scored:

* **pooled** over all mice and 36 targets: pipeline A, GFP+ area ÷ the animal's mean target;
  pipeline B, ΣGFP+ cells ÷ ΣNeuN cells
* **within each mouse** (6 targets): A, GFP+ area summed; B, ΣGFP+ ÷ ΣNeuN

Excluded and marked-wrong sections are never used, whatever their rank. Kept sections missing
from the ranking go last. A mouse with fewer than k sections uses all of them, so the curve
past k = 4 is driven by Mice 2 and 6 alone, and past k = 7 by Mouse 2 alone. At the largest k
every kept section is in, which reproduces ``results/histology/slices/README.md``.

The 90% intervals on r use Fisher's z with targets treated as independent. They are not
(six targets share each mouse), so the intervals are too narrow, if anything.

Writes ``results/histology/slices/top_k_sections.csv`` and ``top_k_sections.png``.
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
import slices_dose_delivery as sdd  # noqa: E402  (mouse list, cumulative sums, palette)

ROOT = sdd.ROOT
QUALITY = ROOT / "data/section_quality.csv"
OUT = sdd.OUT

POOLED = {"A": ("gfp_area_mm2_rel", "A · GFP+ pixels (area ÷ animal's mean target)"),
          "B": ("gfp_fraction", "B · GFP+ neurons (ΣGFP+ ÷ ΣNeuN)")}
WITHIN = {"A": "gfp_area_mm2", "B": "gfp_fraction"}
PIPE_STYLE = {"A": dict(color=sdd.INK, linestyle="-", marker="o"),
              "B": dict(color=sdd.INK_2, linestyle="--", marker="s")}


def ranked_sections(d: pd.DataFrame) -> dict[int, list[str]]:
    """Each mouse's kept sections, best first."""
    rank = pd.read_csv(QUALITY).set_index(["mouse", "section"])["rank"].to_dict()
    kept = d[d.dropped == ""].groupby("mouse").section.unique()
    return {m: sorted(kept[m], key=lambda s: (rank.get((m, s), np.inf), s)) for m in sdd.MICE}


def fisher_ci(r: float, n: int, level: float = 0.90) -> tuple[float, float]:
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    q = stats.norm.ppf(0.5 + level / 2)
    return float(np.tanh(z - q * se)), float(np.tanh(z + q * se))


def trace(d: pd.DataFrame, order: dict[int, list[str]]) -> pd.DataFrame:
    rows = []
    for k in range(1, max(len(v) for v in order.values()) + 1):
        use = {(m, s) for m, secs in order.items() for s in secs[:k]}
        c = sdd.cumulative(d[[(m, s) in use for m, s in zip(d.mouse, d.section)]])
        n_sec = {m: min(k, len(order[m])) for m in sdd.MICE}
        for pipe, (col, _) in POOLED.items():
            r = stats.pearsonr(c.dose, c[col]).statistic
            lo, hi = fisher_ci(r, len(c))
            rows.append({"k": k, "scope": "pooled", "mouse": "", "pipeline": pipe, "measure": col,
                         "n_targets": len(c), "n_sections": sum(n_sec.values()), "r": r,
                         "r_lo90": lo, "r_hi90": hi,
                         "rho": stats.spearmanr(c.dose, c[col]).statistic})
        for m in sdd.MICE:
            cm = c[c.mouse == m]
            for pipe, col in WITHIN.items():
                rows.append({"k": k, "scope": "within mouse", "mouse": m, "pipeline": pipe,
                             "measure": col, "n_targets": len(cm), "n_sections": n_sec[m],
                             "r": stats.pearsonr(cm.dose, cm[col]).statistic,
                             "r_lo90": np.nan, "r_hi90": np.nan,
                             "rho": stats.spearmanr(cm.dose, cm[col]).statistic})
    return pd.DataFrame(rows)


def figure(t: pd.DataFrame, order: dict[int, list[str]]):
    sdd._style()
    ks = sorted(t.k.unique())
    pooled = t[t.scope == "pooled"]
    within = t[t.scope == "within mouse"]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.6), sharey=True,
                             gridspec_kw={"width_ratios": [1.25, 1, 1]})
    fig.subplots_adjust(left=0.05, right=0.99, top=0.8, bottom=0.2, wspace=0.08)
    totals = pooled[pooled.pipeline == "A"].set_index("k").n_sections

    for ax in axes:
        sdd._clean(ax)
        ax.axhline(0, color=sdd.AXIS, linewidth=1, zorder=1)
        ax.set_xticks(ks)
        ax.set_xlim(0.5, ks[-1] + 0.5)
        ax.tick_params(labelsize=9.5, labelcolor=sdd.INK_2)
    axes[0].set_xticklabels([f"{k}\n({totals[k]})" for k in ks])
    axes[0].set_xlabel("Top k sections per mouse\n(sections used, all mice)", fontsize=10.5,
                       color=sdd.INK, labelpad=6)

    ax = axes[0]
    end = pooled[pooled.k == ks[-1]].set_index("pipeline").r
    nudge = {p: (7 if end[p] >= end.drop(p).max() else -7) for p in end.index}  # keep labels apart
    for i, (pipe, (_, label)) in enumerate(POOLED.items()):
        p = pooled[pooled.pipeline == pipe].sort_values("k")
        x = p.k + (i - 0.5) * 0.16
        ax.errorbar(x, p.r, yerr=[p.r - p.r_lo90, p.r_hi90 - p.r], fmt="none",
                    ecolor=PIPE_STYLE[pipe]["color"], elinewidth=1.2, alpha=0.45, zorder=2)
        ax.plot(x, p.r, linewidth=2, markersize=7, markeredgecolor=sdd.SURFACE,
                markeredgewidth=1.5, zorder=3, label=label, **PIPE_STYLE[pipe])
        last = p.iloc[-1]
        ax.annotate(f"{pipe}  r = {last.r:+.2f}", (x.iloc[-1], last.r), xytext=(6, nudge[pipe]),
                    textcoords="offset points", va="center", fontsize=9, color=sdd.INK_2)
    ax.set_title("Pooled over all mice (36 targets)", loc="left", fontsize=11, color=sdd.INK)
    ax.set_ylabel("Pearson r, delivery vs acoustic dose", fontsize=10.5, color=sdd.INK, labelpad=6)
    ax.set_xlim(0.5, ks[-1] + 1.6)
    ax.legend(loc="lower left", frameon=False, fontsize=8.8, handlelength=2.6)

    for ax, pipe, title in ((axes[1], "A", "Within each mouse · A, GFP+ pixels (area)"),
                            (axes[2], "B", "Within each mouse · B, ΣGFP+ ÷ ΣNeuN")):
        for m in sdd.MICE:
            w = within[(within.pipeline == pipe) & (within.mouse == m)].sort_values("k")
            w = w[w.k <= len(order[m])]  # stop where the mouse runs out of sections
            ax.plot(w.k, w.r, color=sdd.MOUSE_COLOR[m], linewidth=2, marker="o", markersize=6,
                    markeredgecolor=sdd.SURFACE, markeredgewidth=1.2, zorder=3,
                    label=f"Mouse {m} ({len(order[m])})")
            ax.annotate(f"M{m}", (w.k.iloc[-1], w.r.iloc[-1]), xytext=(5, 0),
                        textcoords="offset points", va="center", fontsize=8.5, color=sdd.INK_2)
        ax.set_title(title, loc="left", fontsize=11, color=sdd.INK)
        ax.set_xlabel("Top k sections of that mouse", fontsize=10.5, color=sdd.INK, labelpad=6)
    axes[1].legend(loc="lower right", ncol=2, frameon=False, fontsize=8.5,
                   title="sections kept", title_fontsize=8.5)
    lo = min(-0.1, np.floor(np.nanmin(t[["r", "r_lo90"]].to_numpy()) * 10) / 10)
    axes[0].set_ylim(lo, 1.02)
    axes[0].set_yticks(np.round(np.arange(0, 1.001, 0.2), 1))

    sdd._header(fig, "Dose vs delivery using only each mouse's best sections",
                "Sections ranked by hand for quality (data/section_quality.csv); k = how many of "
                "each mouse's best kept sections are summed. Rightmost point = every kept section.")
    sdd._footer(fig, "A mouse with fewer than k sections uses all of them, so past k = 4 only Mice 2 "
                "and 6 change, and past k = 7 only Mouse 2. Bars: 90% interval on pooled r (Fisher z, "
                "targets treated as independent, which they are not). Within a mouse, r is over 6 "
                "targets. Excluded sections are never used.")
    return fig


def main() -> int:
    d = sdd.table("pixels")
    order = ranked_sections(d)
    for m, secs in order.items():
        print(f"Mouse {m}: " + ", ".join(secs))
    t = trace(d, order)
    t.round(4).to_csv(OUT / "top_k_sections.csv", index=False)
    fig = figure(t, order)
    fig.savefig(OUT / "top_k_sections.png", dpi=170)
    plt.close(fig)

    p = t[t.scope == "pooled"].pivot(index="k", columns="pipeline", values="r")
    p.insert(0, "sections", t[(t.scope == "pooled") & (t.pipeline == "A")].set_index("k").n_sections)
    print("\nPooled r by k (36 targets):")
    print(p.round(2).to_string())
    w = t[t.scope == "within mouse"].pivot_table(index="k", columns=["pipeline", "mouse"], values="r")
    print("\nWithin-mouse r by k:")
    print(w.round(2).to_string())
    print(f"\nwrote top_k_sections.csv, top_k_sections.png to {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
