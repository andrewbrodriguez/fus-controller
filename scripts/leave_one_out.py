"""Leave-one-out sensitivity of the pooled dose-delivery result: each mouse, then each target.

    .venv/bin/python scripts/leave_one_out.py

Uses every kept section, summed per target (``slices_dose_delivery.cumulative``), and the two
pooled measures: A, GFP+ area ÷ the animal's mean target; B, ΣGFP+ ÷ ΣNeuN.

**Each mouse.** Drop it and recompute pooled Pearson r on the other 30 targets. A mouse
moves r more than any 6 points would if dropping it shifts r further than dropping 6 random
targets does (20,000 draws; ``p_random6`` is the share of draws that move r at least as far).
Also: fit the other five mice and predict the dropped one. ``resid_B_pp`` is its mean miss on
B in percentage points (a level effect; A is rescaled to each mouse's mean, so it has none),
and ``slope_ratio`` is its own dose slope over the others' (1 = same steepness).

**Each target** (each dose given). Drop it and recompute r. ``cooks_d`` is Cook's distance
in the pooled least-squares fit; above 4/n is the usual flag for an influential point.

p-values and Cook's distance treat targets as independent. They are not (six per mouse).

Writes ``results/histology/slices/leave_one_out.csv`` and ``leave_one_out.png``.
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
OUT = sdd.OUT
MEASURES = {"A": ("gfp_area_mm2_rel", "A · GFP+ pixels (area ÷ animal's mean)"),
            "B": ("gfp_fraction", "B · GFP+ neurons (ΣGFP+ ÷ ΣNeuN)")}
N_DRAWS = 20_000


def r_of(c: pd.DataFrame, col: str) -> float:
    return float(stats.pearsonr(c.dose, c[col]).statistic)


def cooks_distance(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    X = np.column_stack([np.ones_like(x), x])
    h = np.diag(X @ np.linalg.pinv(X.T @ X) @ X.T)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    e = y - X @ beta
    s2 = e @ e / (len(y) - 2)
    return e**2 / (2 * s2) * h / (1 - h) ** 2


def by_mouse(c: pd.DataFrame, rng) -> pd.DataFrame:
    rows = []
    random6 = {}
    for pipe, (col, _) in MEASURES.items():
        r_all = r_of(c, col)
        draws = np.empty(N_DRAWS)
        for i in range(N_DRAWS):
            keep = np.ones(len(c), bool)
            keep[rng.choice(len(c), 6, replace=False)] = False
            draws[i] = r_of(c[keep], col) - r_all
        random6[pipe] = draws
        for m in sdd.MICE:
            rest, one = c[c.mouse != m], c[c.mouse == m]
            dr = r_of(rest, col) - r_all
            fit = stats.linregress(rest.dose, rest[col])
            own = stats.linregress(one.dose, one[col])
            rows.append({
                "left_out": f"Mouse {m}", "kind": "mouse", "mouse": m, "target": "", "dose": np.nan,
                "pipeline": pipe, "r_all": r_all, "r_without": r_all + dr, "delta_r": dr,
                "p_random6": float(np.mean(np.abs(draws) >= abs(dr))),
                "resid_B_pp": (100 * float(np.mean(one[col] - (fit.intercept + fit.slope * one.dose)))
                               if pipe == "B" else np.nan),
                "slope_ratio": own.slope / fit.slope, "within_r": own.rvalue, "cooks_d": np.nan})
    return pd.DataFrame(rows), random6


def by_target(c: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pipe, (col, _) in MEASURES.items():
        r_all = r_of(c, col)
        cd = cooks_distance(c.dose.to_numpy(float), c[col].to_numpy(float))
        for i, t in enumerate(c.itertuples()):
            rows.append({"left_out": f"M{t.mouse} T{t.target}", "kind": "target", "mouse": t.mouse,
                         "target": t.target, "dose": t.dose, "pipeline": pipe, "r_all": r_all,
                         "r_without": r_of(c.drop(index=c.index[i]), col),
                         "delta_r": r_of(c.drop(index=c.index[i]), col) - r_all,
                         "p_random6": np.nan, "resid_B_pp": np.nan, "slope_ratio": np.nan,
                         "within_r": np.nan, "cooks_d": cd[i], "value": getattr(t, col)})
    return pd.DataFrame(rows)


def figure(mice: pd.DataFrame, targets: pd.DataFrame, random6: dict, n: int):
    sdd._style()
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.8),
                             gridspec_kw={"width_ratios": [1, 1.1, 1.1]})
    fig.subplots_adjust(left=0.07, right=0.99, top=0.8, bottom=0.17, wspace=0.28)

    ax = axes[0]
    sdd._clean(ax, grid_axis="x")
    ys = np.arange(len(sdd.MICE))[::-1]
    for i, (pipe, (_, label)) in enumerate(MEASURES.items()):
        sub = mice[mice.pipeline == pipe].set_index("mouse").loc[list(sdd.MICE)]
        r_all = sub.r_all.iloc[0]
        lo, hi = np.percentile(random6[pipe], [5, 95]) + r_all
        off = (0.5 - i) * 0.22
        ax.axvspan(lo, hi, color=sdd.GRID, alpha=0.6 if i == 0 else 0.35, zorder=0, linewidth=0)
        ax.axvline(r_all, color=[sdd.INK, sdd.INK_2][i], linestyle=["-", "--"][i], linewidth=1.2,
                   zorder=1)
        ax.scatter(sub.r_without, ys + off, s=70, zorder=3,
                   color=[sdd.MOUSE_COLOR[m] for m in sub.index], marker=["o", "s"][i],
                   edgecolors=sdd.SURFACE, linewidths=1.5, label=label)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"without Mouse {m}" for m in sdd.MICE], fontsize=9.5, color=sdd.INK_2)
    ax.set_xlabel("Pooled r on the other 30 targets", fontsize=10.5, color=sdd.INK, labelpad=6)
    ax.set_title("Drop each mouse", loc="left", fontsize=11, color=sdd.INK, pad=24)
    leg = ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False, fontsize=8.3,
                    ncol=2, handletextpad=0.2, columnspacing=0.8, borderaxespad=0.2)
    for h in leg.legend_handles:
        h.set_color(sdd.MUTED)

    for ax, pipe in ((axes[1], "A"), (axes[2], "B")):
        sdd._clean(ax)
        t = targets[targets.pipeline == pipe]
        ax.axhline(0, color=sdd.AXIS, linewidth=1, zorder=1)
        lo, hi = np.percentile(random6[pipe] / 6, [5, 95])  # rough per-point scale, for reference
        for m in sdd.MICE:
            tm = t[t.mouse == m]
            ax.scatter(tm.dose, tm.delta_r, s=np.where(tm.cooks_d > 4 / n, 120, 55),
                       color=sdd.MOUSE_COLOR[m], edgecolors=sdd.SURFACE, linewidths=1.5, zorder=3,
                       label=f"Mouse {m}")
        top = t.reindex(t.delta_r.abs().sort_values(ascending=False).index).head(4)
        for p in top.itertuples():
            ax.annotate(p.left_out, (p.dose, p.delta_r), xytext=(7, 0), textcoords="offset points",
                        va="center", fontsize=8.5, color=sdd.INK_2)
        ax.set_title(f"Drop each target · {MEASURES[pipe][1]}", loc="left", fontsize=11,
                     color=sdd.INK)
        ax.set_xlabel("Dose of the dropped target (cumulative 2nd harmonic)", fontsize=10.5,
                      color=sdd.INK, labelpad=6)
        ax.set_ylabel("Change in pooled r when dropped", fontsize=10.5, color=sdd.INK, labelpad=4)
        ax.set_xlim(-0.15, t.dose.max() * 1.18)
    axes[1].legend(loc="upper left", ncol=3, frameon=False, fontsize=8.3, handletextpad=0.2,
                   columnspacing=0.8)

    sdd._header(fig, "Is any mouse or any single dose driving the pooled result?",
                "Left: pooled r without each mouse; line = all 36 targets, band = middle 90% of r "
                "after dropping 6 random targets. Middle, right: change in r when one target is "
                "dropped; large dots have Cook's distance > 4/n.")
    sdd._footer(fig, "Every kept section is summed per target. Targets are treated as independent "
                "in the random-6 band and Cook's distance; they are not (six per mouse). Mouse 1 "
                "T3 is the no-FUS control (dose 0); Mouse 1 T1 is three sonications summed.")
    return fig


def main() -> int:
    c = sdd.cumulative(sdd.table("pixels")).reset_index(drop=True)
    rng = np.random.default_rng(0)
    mice, random6 = by_mouse(c, rng)
    targets = by_target(c)
    pd.concat([mice, targets], ignore_index=True).round(4).to_csv(OUT / "leave_one_out.csv",
                                                                   index=False)
    fig = figure(mice, targets, random6, len(c))
    fig.savefig(OUT / "leave_one_out.png", dpi=170)
    plt.close(fig)

    for pipe, (col, label) in MEASURES.items():
        m = mice[mice.pipeline == pipe]
        print(f"\n{label}: r = {m.r_all.iloc[0]:+.2f} on {len(c)} targets; "
              f"dropping 6 random targets moves it {np.percentile(random6[pipe], 5):+.2f} to "
              f"{np.percentile(random6[pipe], 95):+.2f} (middle 90%)")
        print(m[["left_out", "r_without", "delta_r", "p_random6", "within_r", "slope_ratio",
                 "resid_B_pp"]].round(3).to_string(index=False))
        t = targets[targets.pipeline == pipe]
        t = t.reindex(t.delta_r.abs().sort_values(ascending=False).index)
        print(f"  targets with Cook's D > 4/n = {4 / len(c):.2f}: "
              + (", ".join(f"{d} ({v:.2f})" for d, v in zip(t[t.cooks_d > 4 / len(c)].left_out,
                                                         t[t.cooks_d > 4 / len(c)].cooks_d))
                 or "none"))
        print(t.head(6)[["left_out", "dose", "value", "r_without", "delta_r", "cooks_d"]]
              .round(3).to_string(index=False))
    print(f"\nwrote leave_one_out.csv, leave_one_out.png to {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
