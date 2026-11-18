"""Within-mouse dose trends, how much mice differ, and whether that tracks the acoustics.

    .venv/bin/python scripts/mouse_effect.py      # run after slices_dose_delivery.py

Delivery is the fraction of neurons GFP+ (ΣGFP+ ÷ ΣNeuN over every kept section, pipeline B),
on a logit scale, against log dose (cumulative 2nd harmonic), over the 35 sonicated targets.
Mouse 1's no-FUS control (dose 0) is left out, since log dose is undefined there. This is the
same fit as ``region_delivery.py``: one offset per mouse and one shared slope, least squares.

* **Within each mouse:** r and ρ over its own targets, and its own slope. A test of whether
  separate slopes fit better than the shared one (F test).
* **Between mice:** each mouse's offset, shown as % GFP+ predicted at dose 1. The
  repeated-measures correlation (dose vs delivery after removing each mouse's mean) and the
  share of the variance that is between-mouse level.
* **Exploratory, six points:** does a mouse's level track its acoustics? Both covariates come
  from the one treatment every mouse got (120 bursts, goal band 0.70–0.80), read from the
  recordings: its mean drive voltage over sonicating bursts, and its dose.

GFP+ *area* isn't used for levels: it is summed over sections, and mice have 3 to 12, so its
level is partly section count. It appears only as a within-mouse r.

Writes ``results/histology/mouse_effect/``: ``within_mouse.csv``, ``mouse_level.csv``,
``within_mouse_trend.png``, ``mouse_level_vs_acoustics.png`` and ``README.md``.
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

from fus import archive  # noqa: E402
from fus.extract import extract  # noqa: E402

ROOT = sdd.ROOT
OUT = ROOT / "results/histology/mouse_effect"
EPS = 1e-3  # as in region_delivery.py
PCT_TICKS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8]
DOSE_TICKS = [0.3, 0.5, 1, 2, 3]


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def expit(z):
    return 1 / (1 + np.exp(-z))


# ---------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------


def data() -> pd.DataFrame:
    c = sdd.cumulative(sdd.table("pixels"))
    s = c[c.dose > 0].copy()
    s["x"] = np.log(s.dose)
    s["y"] = logit(s.gfp_fraction)
    s["y_area"] = np.log(s.gfp_area_mm2.clip(EPS))
    return s.reset_index(drop=True)


def shared_slope(s: pd.DataFrame) -> dict:
    """logit(fraction) ~ mouse offsets + slope * log dose, by least squares."""
    mice = sorted(s.mouse.unique())
    X = np.column_stack([(s.mouse == m).astype(float) for m in mice] + [s.x])
    b = np.linalg.lstsq(X, s.y, rcond=None)[0]
    resid = s.y - X @ b
    n, p = X.shape
    sigma2 = resid @ resid / (n - p)
    se = np.sqrt(np.diag(sigma2 * np.linalg.inv(X.T @ X)))
    q = stats.t.ppf(0.95, n - p)
    rss_shared = float(resid @ resid)
    rss_sep = sum(float(np.sum((g.y - np.polyval(np.polyfit(g.x, g.y, 1), g.x)) ** 2))
                  for _, g in s.groupby("mouse"))
    k = len(mice)
    F = ((rss_shared - rss_sep) / (k - 1)) / (rss_sep / (n - 2 * k))
    return {"slope": b[-1], "slope_lo": b[-1] - q * se[-1], "slope_hi": b[-1] + q * se[-1],
            "offset": pd.Series(b[:-1], index=mice), "F": F,
            "p_slopes": float(stats.f.sf(F, k - 1, n - 2 * k)), "n": n}


def rmcorr(s: pd.DataFrame, y: str) -> tuple[float, int, float]:
    """Repeated-measures correlation: pooled r after removing each mouse's means."""
    xc = s.x - s.groupby("mouse").x.transform("mean")
    yc = s[y] - s.groupby("mouse")[y].transform("mean")
    r = float(np.sum(xc * yc) / np.sqrt(np.sum(xc**2) * np.sum(yc**2)))
    df = len(s) - s.mouse.nunique() - 1
    t = r * np.sqrt(df / (1 - r**2))
    return r, df, float(2 * stats.t.sf(abs(t), df))


def within_table(s: pd.DataFrame, fit: dict) -> pd.DataFrame:
    rows = []
    for m, g in s.groupby("mouse"):
        rows.append({
            "mouse": m, "n_targets": len(g), "n_sections": int(g.n_slices.max()),
            "dose_min": g.dose.min(), "dose_max": g.dose.max(),
            "dose_fold": g.dose.max() / g.dose.min(),
            "gfp_min": g.gfp_fraction.min(), "gfp_max": g.gfp_fraction.max(),
            "r_cells": stats.pearsonr(g.x, g.y).statistic,
            "rho_cells": stats.spearmanr(g.dose, g.gfp_fraction).statistic,
            "r_pixels": stats.pearsonr(g.x, g.y_area).statistic,
            "own_slope": stats.linregress(g.x, g.y).slope,
            "level_logit": fit["offset"][m],
            "pct_at_dose_0_5": expit(fit["offset"][m] + fit["slope"] * np.log(0.5)),
            "pct_at_dose_1": expit(fit["offset"][m]),
            "pct_at_dose_2": expit(fit["offset"][m] + fit["slope"] * np.log(2)),
        })
    return pd.DataFrame(rows)


def shared_treatment() -> pd.DataFrame:
    """The 120-burst, goal 0.70-0.80 sonication of each mouse: mean drive voltage and dose."""
    rec = pd.read_csv(ROOT / "results/acoustic_eda/recordings.csv")
    ref = rec[(rec.programmed == 120) & (rec.goal_lo == 0.7) & (rec.goal_hi == 0.8)
              & rec.mouse.isin(sdd.MICE)].sort_values(["mouse", "order"]).groupby("mouse").first()
    rows = []
    for m, r in ref.iterrows():
        with archive.paths(int(m), [r.recording]) as paths:
            x = extract(paths[0]).to_dict()
        rows.append({"mouse": int(m), "recording": r.recording, "position": int(r.position),
                     "ref_mean_voltage": x["mean_voltage"], "ref_dose": x["cum_2nd_harmonic_xlsx"],
                     "ref_at_voltage_cap": bool(r.at_v_cap)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def _logit_axis(ax, lo, hi):
    ticks = [t for t in PCT_TICKS if lo <= t <= hi]
    ax.set_yticks(logit(np.array(ticks)))
    ax.set_yticklabels([f"{100 * t:g}%" for t in ticks])
    ax.set_ylim(logit(lo), logit(hi))


def trend_figure(s: pd.DataFrame, fit: dict, w: pd.DataFrame, rm: tuple):
    sdd._style()
    fig, ax = plt.subplots(figsize=(11.2, 6.8))
    fig.subplots_adjust(left=0.08, right=0.72, top=0.82, bottom=0.17)
    sdd._clean(ax)
    for m in sdd.MICE:
        g = s[s.mouse == m]
        col = sdd.MOUSE_COLOR[m]
        gx = np.linspace(g.x.min(), g.x.max(), 50)
        level = w.set_index("mouse").pct_at_dose_1[m]
        ax.plot(gx, fit["offset"][m] + fit["slope"] * gx, color=col, linewidth=2.2, zorder=2,
                label=f"Mouse {m} · {100 * level:.0f}% at dose 1")
        ax.scatter(g.x, g.y, s=75, color=col, edgecolors=sdd.SURFACE, linewidths=1.5, zorder=3)
        ax.annotate(f"M{m}", (gx[-1], fit["offset"][m] + fit["slope"] * gx[-1]), xytext=(9, 0),
                    textcoords="offset points", va="center", fontsize=9, color=sdd.INK_2)
    ax.set_xticks(np.log(DOSE_TICKS))
    ax.set_xticklabels([f"{t:g}" for t in DOSE_TICKS])
    ax.set_xlim(np.log(0.3), np.log(3.6))
    _logit_axis(ax, 0.005, 0.85)
    ax.set_xlabel("Acoustic dose (cumulative 2nd harmonic, log scale)", fontsize=11,
                  color=sdd.INK, labelpad=6)
    ax.set_ylabel("Neurons GFP+ in the target (Σ GFP+ ÷ Σ NeuN, logit scale)", fontsize=11,
                  color=sdd.INK, labelpad=6)
    ax.tick_params(labelsize=9.5, labelcolor=sdd.INK_2)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=8.8,
              handlelength=1.6)
    odds2 = 2 ** fit["slope"]
    lv = w.pct_at_dose_1
    ax.text(1.02, 0.42,
            f"Shared slope: doubling the dose\nmultiplies the odds by {odds2:.1f}\n"
            f"(90%: {2 ** fit['slope_lo']:.1f}–{2 ** fit['slope_hi']:.1f}).\n\n"
            f"Separate slopes per mouse\nfit no better (p {fit['p_slopes']:.2f}).\n\n"
            f"Within-mouse r = {rm[0]:+.2f}\n(p {rm[2]:.0e}).\n\n"
            f"At dose 1, mice range\n{100 * lv.min():.0f}% to {100 * lv.max():.0f}% GFP+.",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.8, color=sdd.INK_2,
            linespacing=1.3)
    sdd._header(fig, "Within each mouse, more dose means more delivery; mice sit at different levels",
                "One dot per target, one line per mouse: its own level, the slope shared by all six.")
    sdd._footer(fig, "35 sonicated targets; Mouse 1's no-FUS control (dose 0) is left out, and its T1 is three "
                "sonications summed.\nEvery kept section summed per target (pipeline B). Least squares on "
                "logit fraction vs log dose with one offset per mouse, as in region_delivery.py.")
    return fig


def covariate_figure(lv: pd.DataFrame, cors: pd.DataFrame):
    sdd._style()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.9), sharey=True)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.76, bottom=0.2, wspace=0.08)
    panels = (("ref_mean_voltage", "Mean drive voltage on the shared treatment (V)"),
              ("ref_dose", "Dose delivered by the shared treatment"))
    for ax, (col, label) in zip(axes, panels):
        sdd._clean(ax)
        for p in lv.itertuples():
            ax.scatter(getattr(p, col), p.level_logit, s=200, color=sdd.MOUSE_COLOR[p.mouse],
                       edgecolors=sdd.SURFACE, linewidths=1.8, zorder=3)
            ax.annotate(f"M{p.mouse}", (getattr(p, col), p.level_logit), xytext=(8, 0),
                        textcoords="offset points", va="center", fontsize=9, color=sdd.INK_2)
        c = cors.set_index("covariate").loc[col]
        ax.set_title(f"Pearson r {c.r:+.2f} (p {c.p:.2f}) · Spearman ρ {c.rho:+.2f}", loc="left",
                     fontsize=10, color=sdd.INK)
        ax.set_xlabel(label, fontsize=10.5, color=sdd.INK, labelpad=6)
        xv = lv[col]
        pad = (xv.max() - xv.min()) * 0.15
        ax.set_xlim(xv.min() - pad, xv.max() + pad * 1.6)
        ax.tick_params(labelsize=9.5, labelcolor=sdd.INK_2)
    _logit_axis(axes[0], 0.04, 0.5)
    axes[0].set_ylabel("Mouse level: % GFP+ predicted at dose 1", fontsize=10.5, color=sdd.INK,
                       labelpad=6)
    sdd._header(fig, "Does a mouse's level track its acoustics? (exploratory, six mice)",
                "Both covariates come from the one treatment every mouse got: 120 bursts, goal band "
                "0.70–0.80, read from the recordings.")
    sdd._footer(fig, "Six points: any pattern here is a hypothesis for new animals, not a finding. "
                "Level = the mouse's offset in the shared-slope fit (within_mouse_trend.png).")
    return fig


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------


def write_readme(s, fit, w, rm_cells, rm_px, shares, lv, cors):
    odds2 = 2 ** fit["slope"]
    hi, lo = w.level_logit.max(), w.level_logit.min()
    lines = [
        "# Within-mouse trends and the mouse effect",
        "",
        "Generated by `scripts/mouse_effect.py`; do not edit by hand. Delivery is the fraction of "
        "neurons GFP+ (ΣGFP+ ÷ ΣNeuN, pipeline B) on a logit scale against log dose, over "
        f"{fit['n']} sonicated targets in {s.mouse.nunique()} mice (Mouse 1's no-FUS control left out).",
        "",
        "Figures: `within_mouse_trend.png`, `mouse_level_vs_acoustics.png` (exploratory).",
        "",
        "## Within each mouse",
        "",
        "| Mouse | Targets | Sections | Dose range | GFP+ neurons | r (cells) | ρ (cells) | r (pixels) | Own slope |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for p in w.itertuples():
        lines.append(f"| {p.mouse} | {p.n_targets} | {p.n_sections} | {p.dose_min:.2f}–{p.dose_max:.2f} "
                     f"({p.dose_fold:.1f}×) | {100 * p.gfp_min:.0f}–{100 * p.gfp_max:.0f}% | "
                     f"{p.r_cells:+.2f} | {p.rho_cells:+.2f} | {p.r_pixels:+.2f} | {p.own_slope:.2f} |")
    lines += [
        "",
        f"- **Shared slope {fit['slope']:.2f}** logit per log dose (90%: {fit['slope_lo']:.2f}–"
        f"{fit['slope_hi']:.2f}): doubling the dose multiplies the odds of a neuron being GFP+ by "
        f"{odds2:.1f}.",
        f"- Separate slopes per mouse fit no better than the shared one (F = {fit['F']:.2f}, "
        f"p {fit['p_slopes']:.2f}). With 5–6 targets per mouse, modest differences could hide.",
        f"- **Repeated-measures correlation** (dose vs delivery after removing each mouse's means): "
        f"cells r {rm_cells[0]:+.2f} (df {rm_cells[1]}, p {rm_cells[2]:.1e}); pixels r "
        f"{rm_px[0]:+.2f} (p {rm_px[2]:.1e}).",
        "- With six targets, a single mouse's r is very uncertain (r = 0.7 has a 90% interval of "
        "roughly 0.1–0.95).",
        "",
        "## Between mice",
        "",
        "Each mouse's level in the shared-slope fit, as % GFP+ predicted at a given dose:",
        "",
        "| Mouse | Dose 0.5 | Dose 1 | Dose 2 |",
        "|---|---|---|---|",
    ]
    for p in w.sort_values("level_logit", ascending=False).itertuples():
        lines.append(f"| {p.mouse} | {100 * p.pct_at_dose_0_5:.0f}% | **{100 * p.pct_at_dose_1:.0f}%** "
                     f"| {100 * p.pct_at_dose_2:.0f}% |")
    lines += [
        "",
        f"- Highest vs lowest mouse at the same dose: odds ratio {np.exp(hi - lo):.1f}; the lowest "
        f"would need {np.exp((hi - lo) / fit['slope']):.1f}× the dose to match the highest. SD of "
        f"the levels {w.level_logit.std(ddof=1):.2f} logit (unshrunk; `burst_ml`'s hierarchical "
        "estimate is smaller).",
        f"- Between-mouse level is {100 * shares['cells']:.0f}% of the variance in logit GFP+ "
        f"fraction and {100 * shares['pixels']:.0f}% in log GFP+ area (the latter partly section "
        "count).",
        "",
        "## Does the level track the acoustics? (exploratory, six mice)",
        "",
        "| Mouse | Shared-treatment recording | Mean drive voltage (V) | Its dose | At voltage cap | Level at dose 1 |",
        "|---|---|---|---|---|---|",
    ]
    for p in lv.itertuples():
        lines.append(f"| {p.mouse} | {p.recording} (T{p.position}) | {p.ref_mean_voltage:.4f} | "
                     f"{p.ref_dose:.3f} | {'yes' if p.ref_at_voltage_cap else 'no'} | "
                     f"{100 * p.pct_at_dose_1:.0f}% |")
    lines += ["", "| Covariate | Pearson r | p | Spearman ρ |", "|---|---|---|---|"]
    for p in cors.itertuples():
        lines.append(f"| {p.covariate} | {p.r:+.2f} | {p.p:.2f} | {p.rho:+.2f} |")
    lines += ["", "Six points: treat any pattern as a hypothesis to test on new animals."]
    (OUT / "README.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    s = data()
    fit = shared_slope(s)
    w = within_table(s, fit)
    rm_cells, rm_px = rmcorr(s, "y"), rmcorr(s, "y_area")
    shares = {"cells": s.groupby("mouse").y.transform("mean").var() / s.y.var(),
              "pixels": s.groupby("mouse").y_area.transform("mean").var() / s.y_area.var()}
    lv = shared_treatment().merge(w[["mouse", "level_logit", "pct_at_dose_1"]], on="mouse")
    cors = pd.DataFrame([{"covariate": col, "r": stats.pearsonr(lv[col], lv.level_logit).statistic,
                          "p": stats.pearsonr(lv[col], lv.level_logit).pvalue,
                          "rho": stats.spearmanr(lv[col], lv.level_logit).statistic}
                         for col in ("ref_mean_voltage", "ref_dose")])

    w.round(4).to_csv(OUT / "within_mouse.csv", index=False)
    lv.round(4).to_csv(OUT / "mouse_level.csv", index=False)
    for name, fig in (("within_mouse_trend.png", trend_figure(s, fit, w, rm_cells)),
                      ("mouse_level_vs_acoustics.png", covariate_figure(lv, cors))):
        fig.savefig(OUT / name, dpi=170)
        plt.close(fig)
    write_readme(s, fit, w, rm_cells, rm_px, shares, lv, cors)

    print(w[["mouse", "n_targets", "dose_fold", "r_cells", "rho_cells", "r_pixels", "own_slope",
             "pct_at_dose_1"]].round(2).to_string(index=False))
    print(f"\nshared slope {fit['slope']:.2f} (90% {fit['slope_lo']:.2f}-{fit['slope_hi']:.2f}); "
          f"doubling dose x{2 ** fit['slope']:.2f} odds; slopes differ? p {fit['p_slopes']:.2f}")
    print(f"rmcorr cells {rm_cells[0]:+.2f} p {rm_cells[2]:.1e}; pixels {rm_px[0]:+.2f} p {rm_px[2]:.1e}; "
          f"between-mouse share cells {shares['cells']:.0%}, pixels {shares['pixels']:.0%}")
    print("\n" + lv.round(4).to_string(index=False))
    print("\n" + cors.round(3).to_string(index=False))
    print(f"\nwrote {', '.join(sorted(p.name for p in OUT.iterdir()))} to {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
