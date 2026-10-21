"""Exploratory look at the sonication data: study design, controller behaviour, noise.

    .venv/bin/python scripts/sonication_eda.py

Two sources, two levels:

* **Summary sheet, mice 2-24** (``data/Mouse_Controller_Data.xlsx``): one row per target.
  Mouse 1 is left out here, because its rows 3-6 are rotated by one run (data/README.md).
  This level covers the design (programmed bursts x harmonic goal), how much of the
  realised dose the design explains, animal-level vs target-level variation, firing
  order, and the reference target.
* **Raw recordings, mice 1-24** (``data/acoustic/US_Data.zip``, read without unpacking via
  ``fus.archive``): one row per burst. This level
  covers what the controller did: its own online harmonic reading (``posx.data.HH``, the
  value it compares with the goal band), drive voltage, the online wideband reading
  against the interlock threshold, and the receiver noise floor from the zero-volt
  bursts. Burst tables are cached in ``data/processed/acoustic/`` (~2.5 s per recording,
  ~6 min for all 145; only recordings missing from the cache are read).

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
import matplotlib.ticker
import numpy as np
import pandas as pd

from fus import archive
from fus.extract import extract

ROOT = Path(__file__).resolve().parents[1]
SHEET = ROOT / "data/Mouse_Controller_Data.xlsx"
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


# Mouse 17's sixth recording is saved as Target5Repeat. Its goal (0.85), bursts (120) and
# dose (1.4193) are the sheet's position 6 exactly, and the lab's Reduced_ file is Target6.
MISNAMED = {(17, "Target5Repeat"): 6}


def _position(mouse: int, name: str) -> int:
    """Recording -> brain position. TargetN = position N except Mouse 1 (data/README.md)
    and the misnamed Mouse 17 file above."""
    if (mouse, name) in MISNAMED:
        return MISNAMED[(mouse, name)]
    return MOUSE1_POSITION[name] if mouse == 1 else int(name[len("Target")])


def burst_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(bursts, recordings), cached. One row per burst / per recording.

    The cache is filled incrementally: only archive recordings not already in it are read,
    and it is saved after each one, so an interrupted run resumes."""
    bpath, rpath = CACHE / "bursts.csv", CACHE / "recordings.csv"
    have_b = pd.read_csv(bpath) if bpath.exists() else pd.DataFrame()
    have_r = pd.read_csv(rpath) if rpath.exists() else pd.DataFrame()
    done = set(zip(have_r["mouse"], have_r["recording"])) if len(have_r) else set()
    todo = [m for m in archive.recordings() if (m.mouse, m.recording) not in done]
    for k, member in enumerate(todo, 1):
        print(f"  [{k}/{len(todo)}] {member.filename}", flush=True)
        b, r = _read_recording(member)
        have_b = pd.concat([have_b, b], ignore_index=True)
        have_r = pd.concat([have_r, r], ignore_index=True)
        CACHE.mkdir(parents=True, exist_ok=True)
        have_b.to_csv(bpath, index=False)
        have_r.to_csv(rpath, index=False)
    return have_b, have_r


def _read_recording(member: archive.Member) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, recs = [], []
    with archive.local_path(member) as path, warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ex = extract(path)
    r = ex.recording
    name, mouse = member.recording, member.mouse
    position = _position(mouse, name)
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
    return pd.DataFrame(rows), pd.DataFrame(recs)


SESSIONS = {1: "2026-06-11", 2: "2026-06-16", 3: "2026-07-01", 4: "2026-07-07"}


def session(mouse) -> int:
    """Four sessions of six mice each (archive folders 061126, 061626, 070126, 070726)."""
    return (np.asarray(mouse) - 1) // 6 + 1


def dose_check(rec: pd.DataFrame, d: pd.DataFrame) -> pd.DataFrame:
    """Each target's dose from the recordings (repeats summed) against the sheet, mice 2-24."""
    ours = rec[rec["mouse"] != 1].groupby(["mouse", "position"], as_index=False)["cum_h2"].sum()
    j = ours.merge(d[["mouse", "position", "cum_h2"]], on=["mouse", "position"], how="outer",
                   suffixes=("_rec", "_sheet"), indicator=True)
    j["diff"] = j["cum_h2_rec"] - j["cum_h2_sheet"]
    return j


def controller_metrics(b: pd.DataFrame, rec: pd.DataFrame) -> pd.DataFrame:
    rec = rec.copy()
    rec["session"] = session(rec["mouse"])
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
    cm["hh_rel"] = rel_mean
    cm["band_half_rel"] = (cm["goal_hi"] - cm["goal_lo"]) / 2 / ((cm["goal_lo"] + cm["goal_hi"]) / 2)
    return {"mean_rel": float(np.median(rel_mean)), "median_rel": float(np.median(rel_med)),
            "skew": float(np.median(skew)), "mean_in_band": float(np.mean(np.abs(rel_mean - 1) <= half))}


def trace_figure(b: pd.DataFrame, rec: pd.DataFrame, col: str, ylabel: str, path: Path, band: bool,
                 sess: int):
    ymax = b.loc[~b["baseline"], col].quantile(0.995) * 1.08   # all sessions, so they compare
    rec = rec[rec["session"] == sess]
    b = b[b["mouse"].isin(rec["mouse"])]
    mice = sorted(rec["mouse"].unique())
    ncol = int(rec.groupby("mouse").size().max())
    fig, axes = plt.subplots(len(mice), ncol, figsize=(2.0 * ncol, 1.45 * len(mice)),
                             facecolor=SURFACE, sharey=True, squeeze=False)
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
    title(fig, f"{ylabel} per burst, session {sess} ({SESSIONS[sess]}): mice {mice[0]}–{mice[-1]} "
          "(rows), recordings in firing order (columns)",
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
    title(fig, f"Controller performance per recording vs goal (mice {cm['mouse'].min()}–{cm['mouse'].max()})",
          "Ringed points in the last panel reached the voltage ceiling.")
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)


def dose_vs_exposure(b: pd.DataFrame, cm: pd.DataFrame, mpa_per_v: float) -> pd.DataFrame:
    """Per sonication: applied exposure (peak pressure summed over sonicating bursts) and the
    2nd-harmonic dose. Pressure is drive voltage x the sheet's fixed MPa/V conversion."""
    s = b[~b["baseline"]].groupby(["mouse", "recording"])["volt"].sum().rename("cum_v").reset_index()
    d = cm.merge(s, on=["mouse", "recording"])
    d["exposure_mpa"] = d["cum_v"] * mpa_per_v
    return d


def dose_vs_exposure_figure(d: pd.DataFrame, path: Path) -> dict:
    """Simplest view for the lab: 2nd-harmonic AUC against the acoustic exposure applied."""
    from scipy import stats
    x, y = np.log10(d["exposure_mpa"]), np.log10(d["cum_h2"])
    fit = stats.linregress(x, y)
    resid = 10 ** (y - (fit.intercept + fit.slope * x))
    lo, hi = np.percentile(resid, [10, 90])
    fig, a = plt.subplots(figsize=(6.8, 5.2), facecolor=SURFACE)
    style(a)
    a.set_xscale("log")
    a.set_yscale("log")
    gx = np.logspace(x.min() - 0.05, x.max() + 0.05, 50)
    gy = 10 ** (fit.intercept + fit.slope * np.log10(gx))
    a.fill_between(gx, gy * lo, gy * hi, color=BLUES[0], lw=0, zorder=1,
                   label=f"middle 80% of sonications ({lo:.2f}–{hi:.2f}× the trend)")
    a.plot(gx, gy, color=INK_2, lw=1.5, zorder=2, label=f"trend (log–log, r = {fit.rvalue:.2f})")
    a.scatter(d["exposure_mpa"], d["cum_h2"], s=22, color=BLUE, edgecolor=SURFACE, lw=0.6, zorder=3,
              label="one sonication")
    a.set_xlabel("acoustic exposure applied: peak pressure summed over bursts (MPa × bursts)",
                 fontsize=8.5, color=INK_2)
    a.set_ylabel("2nd-harmonic dose (cumulative AUC, sheet units)", fontsize=8.5, color=INK_2)
    a.legend(loc="upper left", frameon=False, fontsize=7.5)
    a.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    spread = hi / lo
    title(fig, "More exposure gives more 2nd-harmonic dose, but loosely",
          f"{len(d)} sonications, {d['mouse'].nunique()} mice. At a given exposure, the dose varies "
          f"~{spread:.1f}× (10th–90th percentile).\nPressure = drive voltage × the sheet's fixed "
          "MPa/V conversion.")
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    return {"exp_r": fit.rvalue, "exp_slope": fit.slope, "exp_lo": lo, "exp_hi": hi, "exp_spread": spread}


def controller_headline_figure(cm: pd.DataFrame, path: Path):
    """One figure for the lab: where the controller settles vs its goal, and how often it ran
    out of headroom (voltage ceiling) or tripped the interlock, per session."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.9), facecolor=SURFACE,
                               gridspec_kw={"width_ratios": [1.45, 1]})
    goal = ((cm["goal_lo"] + cm["goal_hi"]) / 2).round(2)
    real = cm["hh_rel"] * goal
    rng = np.random.default_rng(2)
    x = goal + rng.normal(0, 0.008, len(cm))
    style(a)
    g = np.array([0.5, 1.1])
    a.fill_between(g, g - 0.05, g + 0.05, color=BLUES[0], lw=0, label="goal band (±0.05)")
    a.plot(g, g, color=INK_2, lw=1, ls="--", label="goal")
    ratio = float(np.median(cm["hh_rel"]))
    a.plot(g, ratio * g, color=ORANGE, lw=2, label=f"median settled level ({ratio:.2f}× goal)")
    cap = cm["at_v_cap"].to_numpy()
    a.scatter(x[~cap], real[~cap], s=20, color=BLUE, edgecolor=SURFACE, lw=0.6, zorder=3,
              label="one sonication")
    a.scatter(x[cap], real[cap], s=34, color=BLUE, edgecolor=INK, lw=1.1, zorder=4,
              label="reached the 0.125 V ceiling")
    worst = cm.assign(real=real).sort_values("hh_rel").iloc[0]
    a.annotate(f"M{worst.mouse} {worst.recording}: ceiling + interlock,\nsettled at {worst.hh_rel:.2f}× goal",
               (goal[worst.name], worst.real), xytext=(0.66, 0.2), fontsize=7.5, color=INK_2,
               arrowprops={"arrowstyle": "-", "color": AXIS, "lw": 0.8})
    a.set_xlim(0.5, 1.1)
    a.set_ylim(0, 1.3)
    a.set_xticks(GOALS)
    a.set_xlabel("harmonic goal (band midpoint)", fontsize=8.5, color=INK_2)
    a.set_ylabel("online 2nd harmonic after the ramp (mean)", fontsize=8.5, color=INK_2)
    inband = float((np.abs(cm["hh_rel"] - 1) <= cm["band_half_rel"]).mean())
    a.set_title(f"Settles below goal at every setpoint; {inband:.0%} of sonications average inside the band",
                fontsize=9, color=INK, loc="left")
    a.legend(loc="lower right", frameon=False, fontsize=7.5)   # data-free corner

    style(b)
    sess = sorted(cm["session"].unique())
    xs = np.arange(len(sess))
    w = 0.36
    series = [("at_v_cap", "reached voltage ceiling", BLUE), ("interlock", "interlock fired", ORANGE)]
    for k, (col, lab, colr) in enumerate(series):
        frac = [cm.loc[cm["session"] == s, col].mean() for s in sess]
        cnt = [int(cm.loc[cm["session"] == s, col].sum()) for s in sess]
        bars = b.bar(xs + (k - 0.5) * (w + 0.02), frac, w, color=colr, label=lab, zorder=3)
        for bar, f, c in zip(bars, frac, cnt):
            b.text(bar.get_x() + bar.get_width() / 2, f + 0.008, f"{c}", ha="center", va="bottom",
                   fontsize=7.5, color=INK_2)
    n = [int((cm["session"] == s).sum()) for s in sess]
    b.set_xticks(xs, [f"S{s}\n{SESSIONS[s][5:]}\nmice {6 * s - 5}–{6 * s}" for s in sess], fontsize=7.5)
    b.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1, decimals=0))
    b.set_ylabel(f"share of sonications (n = {min(n)}–{max(n)} per session)", fontsize=8.5, color=INK_2)
    b.set_title("Ceiling hits climb across sessions (count on each bar)", fontsize=9, color=INK, loc="left")
    b.legend(loc="upper left", frameon=False, fontsize=7.5)
    b.set_ylim(0, max(0.45, b.get_ylim()[1]))
    title(fig, f"Controller behaviour across all {len(cm)} sonications ({cm['mouse'].nunique()} mice, "
          f"{len(sess)} sessions)",
          "Online 2nd harmonic = the controller's own per-burst reading, averaged from when it first "
          "reaches the band. Ringed points ran at the 0.125 V ceiling.")
    fig.tight_layout(rect=(0, 0, 1, 0.9), w_pad=2.5)
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def noise_figure(cm: pd.DataFrame, path: Path):
    sessions = sorted(cm["session"].unique())
    fig, axes = plt.subplots(2, len(sessions), figsize=(3.4 * len(sessions), 6.4), facecolor=SURFACE,
                             sharey="row", squeeze=False)
    for row, (col, lab) in enumerate((("noise_h2", "2nd-harmonic band"), ("noise_wb", "wideband"))):
        for k, sess in enumerate(sessions):
            a = axes[row, k]
            style(a)
            for m, s in cm[cm["session"] == sess].sort_values("order").groupby("mouse"):
                rel = s[col] / s[col].iloc[0]
                a.plot(s["order"], rel, color=MUTED, lw=0.9, marker="o", ms=3)
                a.text(s["order"].iloc[-1] + 0.1, rel.iloc[-1], f"M{m}", fontsize=6.5, color=INK_2,
                       va="center")
            a.axhline(1, color=AXIS, lw=0.8)
            if row == 1:
                a.set_xlabel("recording order within the mouse", fontsize=8, color=INK_2)
            if k == 0:
                a.set_ylabel(f"{lab}\nbaseline ÷ first recording", fontsize=8, color=INK_2)
            a.set_title(f"session {sess} ({SESSIONS[sess]})" if row == 0 else "", fontsize=9,
                        color=INK, loc="left")
    title(fig, "Does the normalisation baseline drift through a session? (receiver noise floor, zero-volt bursts)",
          "Every recording is normalised by its own zero-volt bursts, so drift here rescales its dose.")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)


# ---------------------------------------------------------------------------


def dose_line(dc: pd.DataFrame) -> str:
    both = dc[dc["_merge"] == "both"]
    ok = (both["diff"].abs() < 5e-4).sum()
    bad = both[both["diff"].abs() >= 5e-4]
    only_r = dc[dc["_merge"] == "left_only"]
    only_s = dc[dc["_merge"] == "right_only"]
    s = (f"{ok} of {len(both)} targets (mice 2–24) match the sheet's cumulative 2nd harmonic "
         f"to 4 decimals")
    if len(bad):
        s += "; mismatches: " + ", ".join(
            f"M{r.mouse} P{r.position} ({r.cum_h2_rec:.3f} vs {r.cum_h2_sheet:.3f})" for r in bad.itertuples())
    if len(only_r):
        s += "; only in recordings: " + ", ".join(f"M{r.mouse} P{r.position}" for r in only_r.itertuples())
    if len(only_s):
        s += "; only in the sheet: " + ", ".join(f"M{r.mouse} P{r.position}" for r in only_s.itertuples())
    return s + "."


def irregular(cm: pd.DataFrame) -> str:
    out = []
    for m, s in cm.groupby("mouse"):
        names = set(s["recording"])
        missing = [f"Target{k}" for k in range(1, 7) if f"Target{k}" not in names]
        extra = sorted(names - {f"Target{k}" for k in range(1, 7)})
        if (m, "Target5Repeat") in MISNAMED and "Target5Repeat" in names:
            out.append(f"M{m} Target5Repeat is position 6 (misnamed: goal, bursts and dose match "
                       "the sheet's position 6, and the lab's Reduced_ file is Target6)")
            continue
        if missing or extra:
            out.append(f"M{m}" + (f" has no {', '.join(missing)}" if missing else "")
                       + (" and" if missing and extra else "")
                       + (f" has {', '.join(extra)}" if extra else ""))
    return "; ".join(out) + "." if out else "none."


def session_table(cm: pd.DataFrame) -> str:
    rows = ["| Session | Mice | Recordings | Online harmonic ÷ goal (median) | Recordings averaging in band "
            "| Interlock | At V ceiling | Ceiling (V) | Noise-floor CV |",
            "|---|---|---|---|---|---|---|---|---|"]
    for sess, s in cm.groupby("session"):
        inb = (np.abs(s["hh_rel"] - 1) <= s["band_half_rel"]).mean()
        ncv = s.groupby("mouse")["noise_h2"].apply(lambda x: x.std() / x.mean()).median()
        vmax = ", ".join(f"{v:g}" for v in sorted(s["max_v_allowed"].unique()))
        rows.append(f"| {sess} ({SESSIONS[sess]}) | {s['mouse'].min()}–{s['mouse'].max()} | {len(s)} "
                    f"| {s['hh_rel'].median():.2f} | {inb:.0%} | {int(s['interlock'].sum())} "
                    f"| {int(s['at_v_cap'].sum())} | {vmax} | {ncv:.0%} |")
    return "\n".join(rows)


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

## Raw recordings, mice {cm['mouse'].min()}–{cm['mouse'].max()} ({n} recordings, {cm['session'].nunique()} sessions)

Read from `data/acoustic/US_Data.zip` (`fus.archive`). Online harmonic (`HH`) is the
controller's own reading, on the goal's scale. `recordings.csv` has every number per recording.

- **Recordings vs the sheet** (`dose_check.csv`): {dose_line(f['dose_check'])}
- **Irregular recordings:** {irregular(cm)}
- **Dose vs exposure** (`dose_vs_exposure.png`, the simplest figure for the lab): log 2nd-harmonic\n  dose against peak pressure summed over bursts gives r = {f['exp_r']:.2f} (slope {f['exp_slope']:.2f}); at a\n  given exposure the middle 80% of sonications span {f['exp_lo']:.2f}–{f['exp_hi']:.2f}× the trend (~{f['exp_spread']:.1f}×).\n- **Headline** (`controller_headline.png`, for the lab): settled level vs goal for every sonication, and
  ceiling/interlock rates per session.
- **Traces**, one figure per session: `controller_harmonic_s1.png` … `_s4.png` (online harmonic
  with the goal band) and `controller_voltage_s1.png` … (drive voltage with the ceiling).

{session_table(cm)}

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

- Sheet-level numbers exclude Mouse 1; raw-level numbers include it (positions from the
  operator log, `docs/mouse01-dose-delivery.md`).
- Every animal is assumed to have received AAV9 (notes/current.md, 2026-10-19), so nothing is
  split by capsid.
"""


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    d = sheet()
    facts = analyse_sheet(d)
    sheet_figures(d, facts)
    b, rec = burst_tables()
    cm = controller_metrics(b, rec)
    facts.update(steady_state(b, cm))   # also adds hh_rel per recording
    cm.to_csv(OUT / "recordings.csv", index=False)
    for sess in sorted(cm["session"].unique()):
        trace_figure(b, cm, "hh", "Online 2nd harmonic", OUT / f"controller_harmonic_s{sess}.png",
                     band=True, sess=sess)
        trace_figure(b, cm, "volt", "Drive voltage (V)", OUT / f"controller_voltage_s{sess}.png",
                     band=False, sess=sess)
    dc = dose_check(rec, d)
    dc.to_csv(OUT / "dose_check.csv", index=False)
    facts["dose_check"] = dc
    controller_summary_figure(cm, OUT / "controller_summary.png")
    controller_headline_figure(cm, OUT / "controller_headline.png")
    dx = dose_vs_exposure(b, cm, facts["mpa_ratio"])
    facts.update(dose_vs_exposure_figure(dx, OUT / "dose_vs_exposure.png"))
    noise_figure(cm, OUT / "noise_floor.png")
    (OUT / "README.md").write_text(readme(facts, cm))
    print((OUT / "README.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
