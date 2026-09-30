"""Mouse 2: acoustic dose vs GFP coverage, once each section's orientation is recorded.

    .venv/bin/python scripts/mouse02_dose_vs_coverage.py

Unlike Mouse 1, the Mouse 2 sections are mounted at any rotation and every target
was sonicated, so nothing inside the image fixes which slot is which target. A
person clicks the front of the brain and the notch on each section:

    .venv/bin/python -m fus.orientation data/processed/histology/Mouse_02/ds4/*.ome.tif \
        data/processed/histology/Mouse_01/ds4/section_s*.ome.tif

which fills ``data/section_orientation.csv``. The histology is then measured
with those directions:

    .venv/bin/python -m fus.histology measure data/processed/histology/Mouse_02/ds4/*.ome.tif \
        --orientation data/section_orientation.csv \
        --csv results/histology/mouse02_template_slots.csv --figures results/histology/figures/mouse02_template

Assumptions (see docs/mouse01-dose-delivery.md for A1-A7):

B1  The notch is on the animal's left for Mouse 2 (N. Todd, 2026-09-29, correcting
    his earlier "bottom right"), and on the right for Mouse 1 (2026-09-18).
B2  Which side of the animal plan +x (targets 1-3) is on is *derived*, not assumed:
    the Mouse 1 notch corners plus the empty target 3 fix it. The script stops if
    the Mouse 1 sections disagree.
B3  TargetN was fired at position N (N. Todd, 2026-09-18: true for every mouse but 1).
B4  Coverage per target is the mean over the usable sections. Section depth is
    unknown, and N. Todd reports the head was rolled -- right-side targets higher,
    left-side lower -- so one section cuts the two sides at different depths.

Writes ``results/histology/mouse02_dose_vs_coverage.csv`` and
``results/histology/mouse02_dose_delivery_scatter.png``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from fus import extract
from fus import orientation, rois
from fus import histology
from fus.histology import slot_to_target

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import mouse01_dose_vs_coverage as m1  # noqa: E402  (shared plot styling)

ACOUSTIC = ROOT / "data" / "acoustic" / "20260611"
SLOTS = ROOT / "results/histology/mouse02_template_slots.csv"
OUT = ROOT / "results/histology/mouse02_dose_vs_coverage.csv"
SCATTER = ROOT / "results/histology/mouse02_dose_delivery_scatter.png"


def main() -> None:
    notes = orientation.read_sheet()
    try:
        plus_x = orientation.plus_x_side(notes)
    except ValueError as e:
        sys.exit(str(e))
    if plus_x is None:
        sys.exit("click the notch side on Mouse 1 s2, s4, s5, s6 first "
                 "(notebooks/ingest_new_histology.ipynb or `python -m fus.orientation`)")

    m2 = {sec: a for (mouse, sec), a in notes.items() if mouse == 2 and not a.exclude}
    missing = [sec for sec, a in m2.items() if not a.complete]
    if missing:
        sys.exit(f"click centre, front and notch side for: {', '.join(missing)} (or tick Exclude)")

    slots = pd.read_csv(SLOTS)
    placed = slots.groupby("section").first() if "anchor_row" in slots else pd.DataFrame()
    finetuned = rois.read_locations()
    stale = []
    for sec, a in m2.items():
        if sec not in placed.index or placed.loc[sec, "placement"] not in ("template", "finetuned"):
            stale.append(sec)
            continue
        r = placed.loc[sec]
        if r.placement == "finetuned":
            saved = finetuned.get((2, sec))
            got = slots[(slots.section == sec) & (slots.channel == "TRITC")].set_index("slot")
            if saved is None or any(abs(got.loc[roi.slot, "row"] - roi.centre[0]) > 1
                                    or abs(got.loc[roi.slot, "col"] - roi.centre[1]) > 1
                                    for roi in saved.values()):
                stale.append(sec)
            continue
        pixel_um = histology.ROI_RADIUS_MM * r.fit_scale * 1000 / r.radius_px
        t = histology.template_fit(a.centre, a.angle_deg, pixel_um)
        if abs(r.anchor_row - t.row0) > 1 or abs(r.anchor_col - t.col0) > 1:
            stale.append(sec)
    if stale:
        sys.exit(f"{SLOTS.name} was not measured with the current template or finetuned ROIs for "
                 f"{', '.join(stale)}; re-run `fus.histology measure --orientation` (see header)")

    flagged = slots[slots.angle_flag].section.unique().tolist()
    wrong = [sec for sec, a in m2.items() if a.review != "ok"]
    use = [sec for sec in m2 if sec not in flagged and sec not in wrong]
    if not use:
        sys.exit("no section is marked ✓ in the review (notebooks/ingest_new_histology.ipynb, step 5)")
    mir = {sec: orientation.is_mirrored(a, plus_x) for sec, a in m2.items()}
    d = slots[slots.section.isin(use)].copy()
    d["target"] = [slot_to_target(s, mir[sec]) for s, sec in zip(d.slot, d.section)]
    cov = d.groupby(["target", "channel"]).coverage.agg(["mean", "min", "max"]).unstack("channel")

    rows = []
    for t in range(1, 7):
        r = extract(ACOUSTIC / f"Mouse_Cntr_02_Target{t}.mat").to_dict()  # B3
        row = {"target": t, "recording": f"Target{t}", "n_sections": len(use),
               "cum_2nd_harmonic": r["cum_2nd_harmonic_xlsx"],
               "bursts_x_goal": r["n_bursts"] * r["harmonic_setpoint"],
               "interlock_fired": r["goal_was_reduced"]}
        for ch in ("TRITC", "FITC"):
            for stat in ("mean", "min", "max"):
                row[f"{ch}_{stat}"] = cov.loc[t, (stat, ch)]
        rows.append(row)
    df = pd.DataFrame(rows)
    df.round(4).to_csv(OUT, index=False)

    fig = plot_scatter(df)
    fig.savefig(SCATTER, dpi=150, facecolor=m1.SURFACE)
    fit = stats.pearsonr(df.cum_2nd_harmonic, df.TRITC_mean)
    print(f"plan +x is the animal's {plus_x}; sections used: {', '.join(use)}")
    print("mirrored: " + ", ".join(f"{k}={'yes' if v else 'no'}" for k, v in mir.items()))
    print(f"flagged by the angle check and dropped: {', '.join(flagged) or 'none'}")
    print(f"not marked ✓ in review, dropped: {', '.join(wrong) or 'none'}")
    print(df[["target", "cum_2nd_harmonic", "bursts_x_goal", "TRITC_mean", "FITC_mean"]]
          .round(3).to_string(index=False))
    print(f"Pearson r (TRITC vs measured dose, 6 targets) = {fit.statistic:.2f} (p = {fit.pvalue:.2f})")
    print(f"wrote {OUT}\nwrote {SCATTER}")


def plot_scatter(df: pd.DataFrame):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    m1._style()
    x, y = df.cum_2nd_harmonic.to_numpy(), df.TRITC_mean.to_numpy()
    fig, ax = plt.subplots(figsize=(7, 5.2), facecolor=m1.SURFACE)
    m1._clean_axes(ax)
    err = np.vstack([y - df.TRITC_min, df.TRITC_max - y])
    ax.errorbar(x, y, yerr=err, fmt="none", ecolor=m1.SERIES["TRITC"], alpha=0.5, elinewidth=1.4)
    ax.scatter(x, y, s=260, color=m1.SERIES["TRITC"], edgecolors=m1.SURFACE, linewidths=2, zorder=3)
    for r in df.itertuples():
        ax.annotate(str(r.target), (r.cum_2nd_harmonic, r.TRITC_mean), ha="center",
                    va="center_baseline", fontsize=10, weight="bold", color="#ffffff", zorder=4)
    r = stats.pearsonr(x, y).statistic
    ax.text(0.03, 0.97, f"r = {r:.2f},  n = 6 targets,  {int(df.n_sections[0])} sections",
            transform=ax.transAxes, ha="left", va="top", fontsize=11, color=m1.INK_2)
    ax.set_xlim(0, x.max() * 1.1)
    ax.set_ylim(-0.03, 1.0)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.set_xlabel("Acoustic dose (cumulative 2nd harmonic)")
    ax.set_ylabel("AAV delivery (TRITC coverage, 1.0 mm ROI)")
    ax.set_title("Mouse 2: acoustic dose vs AAV delivery", loc="left",
                 fontsize=13, weight="bold", pad=12)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    main()
