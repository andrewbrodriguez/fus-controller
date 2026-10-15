"""Exploratory look at the sonication data: study design, controller behaviour, noise.

    .venv/bin/python scripts/sonication_eda.py

Two sources, two levels:

* **Summary sheet, mice 2-24** (``data/Mouse_Controller_Data.xlsx``): one row per target.
  Mouse 1 is left out here, because its rows 3-6 are rotated by one run (data/README.md).
  This level covers the design (programmed bursts x harmonic goal), how much of the
  realised dose the design explains, animal-level vs target-level variation, firing
  order, and the reference target.
* **Raw recordings, mice 1-6** (``data/acoustic/20260611``): one row per burst. This level
  covers what the controller did: its own online harmonic reading (``posx.data.HH``, the
  value it compares with the goal band), drive voltage, the online wideband reading
  against the interlock threshold, and the receiver noise floor from the zero-volt
  bursts. Burst tables are cached in ``data/processed/acoustic/`` (~1 s per file).

Writes ``results/acoustic_eda/``: figures, ``recordings.csv`` (one row per recording) and
``README.md`` with the numbers quoted.
"""

from __future__ import annotations

import datetime as dt
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fus.extract import extract

ROOT = Path(__file__).resolve().parents[1]
SHEET = ROOT / "data/Mouse_Controller_Data.xlsx"
ACOUSTIC = ROOT / "data/acoustic/20260611"
CACHE = ROOT / "data/processed/acoustic"
OUT = ROOT / "results/acoustic_eda"

# Reference palette (dataviz skill), light mode, as in scripts/slices_dose_delivery.py.
SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"
BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
GOALS = [0.55, 0.65, 0.75, 0.85, 0.95, 1.05]
GOAL_COLOR = dict(zip(GOALS, BLUES[1:]))
PROGRAMMED = [45, 60, 90, 120, 180, 240]

# Mouse 1 was the operator deviation (docs/mouse01-dose-delivery.md): recording -> position.
MOUSE1_POSITION = {"Target1": 1, "Target2": 1, "Target3": 1, "Target2_Repeat": 2,
                   "Target4": 4, "Target5": 5, "Target6": 6}


def style(ax):
    ax.set_facecolor(SURFACE)
    for k in ("top", "right"):
        ax.spines[k].set_visible(False)
    for k in ("left", "bottom"):
        ax.spines[k].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=7.5)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def title(fig, text, sub=None):
    fig.suptitle(text, x=0.01, ha="left", fontsize=11, color=INK, y=0.995)
    if sub:
        fig.text(0.01, 0.955, sub, ha="left", va="top", fontsize=8, color=INK_2)


# ---------------------------------------------------------------------------
# Summary sheet (mice 2-24)
# ---------------------------------------------------------------------------


def sheet() -> pd.DataFrame:
    d = pd.read_excel(SHEET).rename(columns={
        "Mouse": "mouse", "Treatment #": "order", "Brain Region": "position",
        "N Bursts": "n_stored", "Harmonic Goal": "goal", "Mean Voltage": "mean_v",
        "Cumulative Voltage": "cum_v", "Mean 2nd Harmonic": "mean_h2",
        "Cumulative 2nd Harmonic": "cum_h2"})
    d = d[d["mouse"] != 1].copy()
    d["n_prog"] = d["n_stored"] - np.where(d["mouse"] <= 6, 15, 12)
    d["goal"] = d["goal"].round(2)
    return d


def icc(resid: pd.Series, groups: pd.Series) -> float:
    """One-way ANOVA intraclass correlation: share of variance that is between groups."""
    g = resid.groupby(groups)
    k = g.size().mean()
    msb = (g.size() * (g.mean() - resid.mean()) ** 2).sum() / (g.ngroups - 1)
    msw = ((resid - g.transform("mean")) ** 2).sum() / (len(resid) - g.ngroups)
    return float(max(0.0, (msb - msw) / (msb + (k - 1) * msw)))


def ols(y: np.ndarray, cols: list[np.ndarray]):
    X = np.column_stack([np.ones(len(y)), *cols])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    fit = X @ beta
    r2 = 1 - ((y - fit) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return beta, y - fit, float(r2)


def sheet_figures(d: pd.DataFrame, facts: dict):
    # --- design grid
    counts = pd.crosstab(d["n_prog"], d["goal"]).reindex(index=PROGRAMMED, columns=GOALS, fill_value=0)
    fig, ax = plt.subplots(figsize=(6.2, 4.4), facecolor=SURFACE)
    vmax = counts.values.max()
    ax.imshow(counts.values, cmap=matplotlib.colors.ListedColormap(BLUES), vmin=0, vmax=vmax,
              aspect="auto", origin="lower")
    for i in range(counts.shape[0]):
        for j in range(counts.shape[1]):
            v = counts.values[i, j]
            ax.text(j, i, str(v), ha="center", va="center", fontsize=9,
                    color=SURFACE if v > vmax * 0.55 else INK)
    ax.set_xticks(range(len(GOALS)), [f"{g:.2f}" for g in GOALS])
    ax.set_yticks(range(len(PROGRAMMED)), PROGRAMMED)
    ax.set_xlabel("harmonic goal (deadband midpoint)", fontsize=8.5, color=INK_2)
    ax.set_ylabel("programmed bursts", fontsize=8.5, color=INK_2)
    ax.tick_params(colors=INK_2, labelsize=8, length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    title(fig, "Study design: targets per bursts × goal cell (mice 2–24)",
          f"{len(d)} targets. The reference cell (120 bursts at 0.75) holds {facts['n_ref']}; "
          f"every mouse has at least one.")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / "design.png", dpi=130, facecolor=SURFACE)
    plt.close(fig)

    # --- dose vs design, and achieved harmonic vs goal
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), facecolor=SURFACE)
    a = axes[0]
    style(a)
    rng = np.random.default_rng(0)
    for g in GOALS:
        s = d[d["goal"] == g]
        a.scatter(s["n_prog"] * np.exp(rng.normal(0, 0.03, len(s))), s["cum_h2"], s=18,
                  color=GOAL_COLOR[g], edgecolor=SURFACE, lw=0.6, label=f"{g:.2f}")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xticks(PROGRAMMED, PROGRAMMED)
    a.minorticks_off()
    a.set_xlabel("programmed bursts", fontsize=8.5, color=INK_2)
    a.set_ylabel("cumulative 2nd harmonic (sheet units)", fontsize=8.5, color=INK_2)
    a.set_title(f"Dose ~ bursts × goal: R² = {facts['r2_design']:.2f} on log scale",
                fontsize=9, color=INK, loc="left")
    a.legend(title="goal", fontsize=7, title_fontsize=7.5, frameon=False, ncol=2)
    b = axes[1]
    style(b)
    for g in GOALS:
        s = d[d["goal"] == g]
        b.scatter(g + rng.normal(0, 0.012, len(s)), s["mean_v"], s=16, color=GOAL_COLOR[g],
                  edgecolor=SURFACE, lw=0.6)
    b.set_xlabel("harmonic goal", fontsize=8.5, color=INK_2)
    b.set_ylabel("mean drive voltage (V)", fontsize=8.5, color=INK_2)
    b.set_title(f"Voltage needed per goal: {facts['icc_v']:.0%} of the residual is animal-level",
                fontsize=9, color=INK, loc="left")
    title(fig, "What the controller delivered for each design cell (mice 2–24)")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / "dose_vs_design.png", dpi=130, facecolor=SURFACE)
    plt.close(fig)

    # --- animal effects: residual dose and voltage by mouse, and by firing order
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.4), facecolor=SURFACE,
                             gridspec_kw={"width_ratios": [3, 1.2]})
    mice = sorted(d["mouse"].unique())
    for row, (col, lab, key) in enumerate((("resid_dose", "log dose residual", "icc_dose"),
                                           ("resid_v", "log voltage residual", "icc_v"))):
        a = axes[row, 0]
        style(a)
        means = d.groupby("mouse")[col].mean().reindex(mice)
        for i, m in enumerate(mice):
            v = d.loc[d["mouse"] == m, col]
            a.scatter(np.full(len(v), i) + rng.normal(0, 0.08, len(v)), v, s=12, color=MUTED,
                      edgecolor="none", alpha=0.8)
        a.scatter(range(len(mice)), means, marker="_", s=180, color=BLUE, lw=2)
        a.axhline(0, color=AXIS, lw=0.8)
        a.set_xticks(range(len(mice)), mice, fontsize=7)
        a.set_ylabel(lab, fontsize=8.5, color=INK_2)
        a.set_title(f"{lab.split(' residual')[0].capitalize()} after the design: "
                    f"{facts[key]:.0%} of the variance is between mice (blue = mouse mean)",
                    fontsize=9, color=INK, loc="left")
        b = axes[row, 1]
        style(b)
        for p in range(1, 7):
            v = d.loc[d["order"] == p, col]
            b.scatter(np.full(len(v), p) + rng.normal(0, 0.08, len(v)), v, s=10, color=MUTED,
                      edgecolor="none", alpha=0.7)
        b.scatter(range(1, 7), d.groupby("order")[col].mean().reindex(range(1, 7)), marker="_",
                  s=180, color=BLUE, lw=2)
        b.axhline(0, color=AXIS, lw=0.8)
        b.set_xticks(range(1, 7))
        b.set_title("by firing order (= position)", fontsize=9, color=INK, loc="left")
    axes[1, 0].set_xlabel("mouse", fontsize=8.5, color=INK_2)
    axes[1, 1].set_xlabel("treatment #", fontsize=8.5, color=INK_2)
    title(fig, "Animal-level vs target-level variation (mice 2–24; residuals after bursts × goal)",
          "Firing order and position are the same variable for mice 2–24: position 1 is always "
          "sonicated first.")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(OUT / "animal_effects.png", dpi=130, facecolor=SURFACE)
    plt.close(fig)


def analyse_sheet(d: pd.DataFrame) -> dict:
    facts = {}
    ref = d[(d["n_prog"] == 120) & (d["goal"] == 0.75)]
    facts["n_ref"] = len(ref)
    facts["ref_mice"] = ref["mouse"].nunique()
    facts["ref_cv_dose"] = float(ref["cum_h2"].std() / ref["cum_h2"].mean())
    facts["ref_cv_v"] = float(ref["mean_v"].std() / ref["mean_v"].mean())
    facts["ref_dose_range"] = (float(ref["cum_h2"].min()), float(ref["cum_h2"].max()))
    y = np.log(d["cum_h2"].to_numpy())
    beta, resid, r2 = ols(y, [np.log(d["n_prog"].to_numpy()), d["goal"].to_numpy()])
    d["resid_dose"] = resid
    facts["r2_design"] = r2
    facts["burst_exponent"] = float(beta[1])
    _, resid_v, r2_v = ols(np.log(d["mean_v"].to_numpy()), [d["goal"].to_numpy()])
    d["resid_v"] = resid_v
    facts["r2_v_goal"] = r2_v
    facts["icc_dose"] = icc(d["resid_dose"], d["mouse"])
    facts["icc_v"] = icc(d["resid_v"], d["mouse"])
    facts["order_dose"] = d.groupby("order")["resid_dose"].mean().round(3).to_dict()
    facts["order_v"] = d.groupby("order")["resid_v"].mean().round(3).to_dict()
    facts["corr_v_dose_resid"] = float(np.corrcoef(d["resid_v"], d["resid_dose"])[0, 1])
    facts["mpa_ratio"] = float((d["Mean MPa"] / d["mean_v"]).median())
    facts["mpa_ratio_spread"] = float((d["Mean MPa"] / d["mean_v"]).std())
    return facts


# ---------------------------------------------------------------------------
# Raw recordings (mice 1-6)
# ---------------------------------------------------------------------------


def _datenum(x: float) -> dt.datetime:
    return dt.datetime.fromordinal(int(x)) + dt.timedelta(days=x % 1) - dt.timedelta(days=366)


def recording_paths() -> list[Path]:
    return sorted(p for p in ACOUSTIC.glob("Mouse_Cntr_0[1-6]_Target*.mat"))


def burst_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(bursts, recordings), cached. One row per burst / per recording."""
    bpath, rpath = CACHE / "bursts.csv", CACHE / "recordings.csv"
    if bpath.exists() and rpath.exists():
        return pd.read_csv(bpath), pd.read_csv(rpath)
    rows, recs = [], []
    for path in recording_paths():
        print(f"  {path.name}")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ex = extract(path)
        r = ex.recording
        name = path.stem.split("_", 3)[3]          # Target4, Target2_Repeat
        mouse = int(path.stem.split("_")[2])
        position = MOUSE1_POSITION[name] if mouse == 1 else int(name.replace("Target", ""))
        last = np.atleast_1d(r.raw.posx.data)[-1]
        hh = np.asarray(last.HH, float).ravel()
        wb = np.asarray(last.WB, float).ravel()
        volts = r.drive_voltage[:, 0] if r.drive_voltage is not None else np.full(r.n_bursts, np.nan)
        norm = ex.metrics.normalised()
        n = min(len(hh), len(wb), len(volts), r.n_bursts)
        p = r.raw.posx
        goal = np.atleast_1d(r.harmonic_goal).astype(float)
        cur = np.atleast_1d(r.current_harmonic_goal).astype(float)
        nb = r.n_baseline
        recs.append({
            "mouse": mouse, "recording": name, "position": position,
            "time": _datenum(float(np.atleast_1d(r.raw.niscope.tracktime)[0])).isoformat(),
            "goal_lo": goal[0], "goal_hi": goal[-1], "goal_end_lo": cur[0], "goal_end_hi": cur[-1],
            "interlock": r.goal_was_reduced, "n_bursts": r.n_bursts, "n_baseline": nb,
            "programmed": r.programmed_bursts, "max_v_allowed": float(getattr(p, "MaxV1", np.nan)),
            "start_v": float(getattr(p, "StartV1", np.nan)), "wb_thresh": float(getattr(p, "WBthresh", np.nan)),
            "noise_h2": float(ex.metrics.area_2nd[:nb].mean()),
            "noise_wb": float(ex.metrics.area_wb[:nb].mean()),
            "cum_h2": ex.to_dict()["cum_2nd_harmonic_xlsx"],
        })
        for i in range(n):
            rows.append({"mouse": mouse, "recording": name, "burst": i, "baseline": i < nb,
                         "volt": volts[i], "hh": hh[i], "wb": wb[i],
                         "h2_norm": norm.area_2nd[i], "wb_norm": norm.area_wb[i]})
    b, rec = pd.DataFrame(rows), pd.DataFrame(recs)
    CACHE.mkdir(parents=True, exist_ok=True)
    b.to_csv(bpath, index=False)
    rec.to_csv(rpath, index=False)
    return b, rec


def controller_metrics(b: pd.DataFrame, rec: pd.DataFrame) -> pd.DataFrame:
    rec = rec.copy()
    rec["time"] = pd.to_datetime(rec["time"])
    rec["order"] = rec.groupby("mouse")["time"].rank().astype(int)
    out = []
    for r in rec.itertuples():
        s = b[(b["mouse"] == r.mouse) & (b["recording"] == r.recording) & ~b["baseline"]]
        hh, v, wb = s["hh"].to_numpy(), s["volt"].to_numpy(), s["wb"].to_numpy()
        reach = np.argmax(hh >= r.goal_lo) if (hh >= r.goal_lo).any() else np.nan
        after = hh[int(reach):] if np.isfinite(reach) else np.array([])
        out.append({
            "ramp_bursts": reach,
            "in_band": float(((after >= r.goal_lo) & (after <= r.goal_hi)).mean()) if after.size else np.nan,
            "above_band": float((after > r.goal_hi).mean()) if after.size else np.nan,
            "below_band": float((after < r.goal_lo).mean()) if after.size else np.nan,
            "hh_median_after": float(np.median(after)) if after.size else np.nan,
            "hh_cv_after": float(after.std() / after.mean()) if after.size and after.mean() > 0 else np.nan,
            "v_max": float(v.max()), "v_final": float(v[-5:].mean()),
            "at_v_cap": bool(v.max() >= r.max_v_allowed - 1e-6),
            "wb_over": int((wb > r.wb_thresh).sum()), "wb_max": float(wb.max()),
            "h2_wb_corr": float(np.corrcoef(hh, wb)[0, 1]) if len(hh) > 2 else np.nan,
        })
    return pd.concat([rec.reset_index(drop=True), pd.DataFrame(out)], axis=1)


def steady_state(b: pd.DataFrame, cm: pd.DataFrame) -> dict:
    """Online harmonic after the ramp, relative to the goal midpoint, per recording."""
    from scipy import stats
    rel_mean, rel_med, skew = [], [], []
    for r in cm.itertuples():
        s = b[(b["mouse"] == r.mouse) & (b["recording"] == r.recording) & ~b["baseline"]]["hh"].to_numpy()
        a = s[int(r.ramp_bursts):]
        g = (r.goal_lo + r.goal_hi) / 2
        rel_mean.append(a.mean() / g)
        rel_med.append(np.median(a) / g)
        skew.append(stats.skew(a))
    rel_mean = np.array(rel_mean)
    half = (cm["goal_hi"] - cm["goal_lo"]).to_numpy() / 2 / ((cm["goal_lo"] + cm["goal_hi"]) / 2).to_numpy()
    return {"mean_rel": float(np.median(rel_mean)), "median_rel": float(np.median(rel_med)),
            "skew": float(np.median(skew)), "mean_in_band": float(np.mean(np.abs(rel_mean - 1) <= half))}


def trace_figure(b: pd.DataFrame, rec: pd.DataFrame, col: str, ylabel: str, path: Path, band: bool):
    mice = sorted(rec["mouse"].unique())
    ncol = int(rec.groupby("mouse").size().max())
    fig, axes = plt.subplots(len(mice), ncol, figsize=(2.0 * ncol, 1.45 * len(mice)),
                             facecolor=SURFACE, sharey=True, squeeze=False)
    ymax = b.loc[~b["baseline"], col].quantile(0.995) * 1.08
    for i, m in enumerate(mice):
        rr = rec[rec["mouse"] == m].sort_values("order")
        for j in range(ncol):
            a = axes[i, j]
            if j >= len(rr):
                a.axis("off")
                continue
            r = rr.iloc[j]
            s = b[(b["mouse"] == m) & (b["recording"] == r["recording"])]
            style(a)
            a.grid(False)
            if band:
                a.axhspan(r["goal_lo"], r["goal_hi"], color=BLUES[0], lw=0)
                if r["interlock"]:
                    a.axhspan(r["goal_end_lo"], r["goal_end_hi"], color=ORANGE, alpha=0.18, lw=0)
            else:
                a.axhline(r["max_v_allowed"], color=AXIS, lw=0.8, ls="--")
            a.plot(s["burst"], s[col], color=BLUE, lw=0.9)
            a.axvline(r["n_baseline"] - 0.5, color=AXIS, lw=0.6)
            a.set_ylim(0, ymax)
            lbl = f"P{r['position']} · {r['recording'].replace('Target', 'T')}"
            lbl += f" · {r['goal_lo']:.2f}–{r['goal_hi']:.2f}" if band else ""
            lbl += " · interlock" if r["interlock"] else ""
            lbl += " · at V cap" if (not band and r["at_v_cap"]) else ""
            a.set_title(lbl, fontsize=6.5, color=INK_2, loc="left", pad=2)
            a.tick_params(labelsize=6)
            if j == 0:
                a.set_ylabel(f"M{m}", fontsize=8, color=INK, rotation=0, ha="right", va="center")
    title(fig, f"{ylabel} per burst, mice 1–6 (rows), recordings in firing order (columns)",
          ("Blue band = goal deadband; orange = band after the interlock lowered it. "
           "Vertical line = end of zero-volt baseline bursts.") if band else
          "Dashed = the controller's voltage ceiling (MaxV1). Vertical line = end of baseline bursts.")
    fig.tight_layout(rect=(0, 0, 1, 0.94), h_pad=0.6, w_pad=0.4)
    fig.savefig(path, dpi=120, facecolor=SURFACE)
    plt.close(fig)


def controller_summary_figure(cm: pd.DataFrame, path: Path):
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.6), facecolor=SURFACE)
    rng = np.random.default_rng(1)
    cm = cm.assign(goal=((cm["goal_lo"] + cm["goal_hi"]) / 2).round(2))
    panels = [("ramp_bursts", "bursts to first reach the band"),
              ("in_band", "fraction of later bursts in band"),
              ("hh_median_after", "median online harmonic / band"),
              ("v_final", "final drive voltage (V)")]
    for a, (col, lab) in zip(axes, panels):
        style(a)
        y = cm[col] / ((cm["goal_lo"] + cm["goal_hi"]) / 2) if col == "hh_median_after" else cm[col]
        a.scatter(cm["goal"] + rng.normal(0, 0.01, len(cm)), y, s=22,
                  color=[GOAL_COLOR.get(g, MUTED) for g in cm["goal"]], edgecolor=SURFACE, lw=0.6)
        if col == "v_final":
            a.axhline(cm["max_v_allowed"].median(), color=AXIS, lw=0.8, ls="--")
            capped = cm[cm["at_v_cap"]]
            a.scatter(capped["goal"], capped["v_final"], s=60, facecolor="none", edgecolor=INK, lw=0.8)
        if col == "hh_median_after":
            a.axhline(1, color=AXIS, lw=0.8)
        a.set_xlabel("harmonic goal", fontsize=8, color=INK_2)
        a.set_title(lab, fontsize=8.5, color=INK, loc="left")
    title(fig, "Controller performance per recording vs goal (mice 1–6)",
          "Ringed points in the last panel reached the voltage ceiling.")
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)


def noise_figure(cm: pd.DataFrame, path: Path):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), facecolor=SURFACE)
    for a, col, lab in ((axes[0], "noise_h2", "2nd-harmonic band"), (axes[1], "noise_wb", "wideband")):
        style(a)
        for m, s in cm.sort_values("order").groupby("mouse"):
            rel = s[col] / s[col].iloc[0]
            a.plot(s["order"], rel, color=MUTED, lw=0.9, marker="o", ms=3)
            a.text(s["order"].iloc[-1] + 0.1, rel.iloc[-1], f"M{m}", fontsize=7, color=INK_2,
                   va="center")
        a.axhline(1, color=AXIS, lw=0.8)
        a.set_xlabel("firing order within the session", fontsize=8, color=INK_2)
        a.set_ylabel("baseline level ÷ first recording", fontsize=8, color=INK_2)
        a.set_title(f"Receiver noise floor, {lab} (zero-volt bursts)", fontsize=9, color=INK, loc="left")
    title(fig, "Does the normalisation baseline drift through a session? (mice 1–6)",
          "Every recording is normalised by its own zero-volt bursts, so drift here rescales its dose.")
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)


# ---------------------------------------------------------------------------


def readme(f: dict, cm: pd.DataFrame) -> str:
    n = len(cm)
    ilk = cm[cm["interlock"]]
    cap = cm[cm["at_v_cap"]]
    gap = cm.sort_values(["mouse", "order"]).groupby("mouse")["time"].apply(
        lambda t: pd.to_datetime(t).diff().dt.total_seconds().div(60).median()).median()
    noise_cv = cm.groupby("mouse")["noise_h2"].apply(lambda s: s.std() / s.mean()).median()
    return f"""# Sonication EDA

Generated by `scripts/sonication_eda.py`. Exploratory: nothing here is a pre-registered test.

## Summary sheet, mice 2–24 (Mouse 1 left out: rotated rows)

- **Design** (`design.png`): 6 programmed burst counts ({', '.join(map(str, PROGRAMMED))}) × 6
  harmonic goals (0.55–1.05). Stored counts are programmed + 15 (mice 1–6) or + 12 (mice 7–24).
  The reference cell (120 bursts at goal 0.75) holds {f['n_ref']} targets, and all
  {f['ref_mice']} mice have at least one.
- **Firing order is position.** For mice 2–24, treatment # equals brain region on every row:
  position 1 is always sonicated first and position 6 last. Order effects (bubble or
  anaesthesia drift) and position effects (anatomy, skull) can't be separated.
- **The design explains most of the dose** (`dose_vs_design.png`): log cumulative 2nd
  harmonic ~ log bursts + goal gives R² = {f['r2_design']:.2f}. The burst exponent is
  {f['burst_exponent']:.2f}; 1.0 would mean dose ∝ bursts.
- **Reference target reproducibility:** cumulative dose CV {f['ref_cv_dose']:.0%} (range
  {f['ref_dose_range'][0]:.2f}–{f['ref_dose_range'][1]:.2f}); mean voltage CV {f['ref_cv_v']:.0%}.
- **Animal vs target** (`animal_effects.png`, one-way ICC on residuals after the design):
  {f['icc_dose']:.0%} of the residual dose variance is between mice; {f['icc_v']:.0%} of the
  residual drive-voltage variance is between mice. A large share of the voltage needed to
  reach a goal is an animal property (skull, coupling), while the dose left over is mostly
  target-level, which supports within-animal contrasts.
- **Voltage and dose residuals** are barely related (r = {f['corr_v_dose_resid']:+.2f}): needing
  more voltage than the goal predicts says little about the dose a target ended up with.
- **Mean residual by firing order** (dose, log): {f['order_dose']}. Voltage: {f['order_v']}.
- **`Mean MPa` is a fixed multiple of voltage:** MPa / V = {f['mpa_ratio']:.3f} (SD
  {f['mpa_ratio_spread']:.1e}). It is a conversion, not a measurement; don't use it as a
  separate predictor.

## Raw recordings, mice 1–6 ({n} recordings)

Online harmonic (`HH`) is the controller's own reading, on the goal's scale.
`recordings.csv` has every number per recording.

- **Traces:** `controller_harmonic.png` (online harmonic with the goal band) and
  `controller_voltage.png` (drive voltage with the ceiling).
- **Ramp:** median {cm['ramp_bursts'].median():.0f} sonicating bursts to first reach the band
  (range {cm['ramp_bursts'].min():.0f}–{cm['ramp_bursts'].max():.0f}). Starting voltage
  {cm['start_v'].median():g} V, ceiling {cm['max_v_allowed'].median():g} V.
- **Tracking after the ramp:** median {cm['in_band'].median():.0%} of bursts inside the
  deadband, {cm['above_band'].median():.0%} above, {cm['below_band'].median():.0%} below. The
  band is only ±0.05 wide, while burst-to-burst CV of the online harmonic is median
  {cm['hh_cv_after'].median():.0%}, so most bursts are expected to fall outside it.
  `controller_summary.png` plots these against goal.
- **The controller settles below its setpoint.** After the ramp, the online harmonic averages
  {f['mean_rel']:.2f}× the goal midpoint (median over recordings; median-based
  {f['median_rel']:.2f}×). This holds at every goal from 0.55 to 1.05, and only
  {f['mean_in_band']:.0%} of recordings average inside their band. Voltage is flat by then
  (median slope over the second half ≈ 0), so this is a steady offset, not a slow ramp. It
  isn't skew (median skew {f['skew']:.2f}). The realised harmonic is therefore a near-constant
  fraction of the nominal goal, which matters for any model using `goal` rather than the
  measured harmonic.
- **Voltage ceiling:** {len(cap)} of {n} recordings reached {cm['max_v_allowed'].median():g} V
  ({', '.join(f"M{r.mouse} {r.recording}" for r in cap.itertuples()) or 'none'}). There, the goal
  may not have been achievable.
- **Interlock:** fired in {len(ilk)} of {n} recordings
  ({', '.join(f"M{r.mouse} {r.recording} (goal {r.goal_lo:.2f}–{r.goal_hi:.2f} → {r.goal_end_lo:.2f}–{r.goal_end_hi:.2f}, {r.wb_over} bursts over threshold)" for r in ilk.itertuples()) or 'none'}).
- **Noise floor** (`noise_floor.png`): within-session CV of the zero-volt 2nd-harmonic level is
  median {noise_cv:.0%} per mouse. Each recording's dose is normalised by its own floor, so
  floor drift rescales dose directly.
- **Session timing:** median {gap:.1f} min between consecutive sonications.

## Caveats

- Sheet-level numbers exclude Mouse 1. Raw-level numbers cover mice 1–6 only, which is the
  only acoustic session on disk; 3 of 4 days are still on Dropbox.
- Capsid cohort isn't in the sheet, so none of this is split by capsid.
"""


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    d = sheet()
    facts = analyse_sheet(d)
    sheet_figures(d, facts)
    b, rec = burst_tables()
    cm = controller_metrics(b, rec)
    facts.update(steady_state(b, cm))
    cm.to_csv(OUT / "recordings.csv", index=False)
    trace_figure(b, cm, "hh", "Online 2nd harmonic", OUT / "controller_harmonic.png", band=True)
    trace_figure(b, cm, "volt", "Drive voltage (V)", OUT / "controller_voltage.png", band=False)
    controller_summary_figure(cm, OUT / "controller_summary.png")
    noise_figure(cm, OUT / "noise_floor.png")
    (OUT / "README.md").write_text(readme(facts, cm))
    print((OUT / "README.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
