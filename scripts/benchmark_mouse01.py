"""Benchmark a placement method on Mouse 1, the one animal with an answer key.

    .venv/bin/python scripts/benchmark_mouse01.py

Mouse 1 is the reference because three things about it are known independently of
any placement method:

* target 3 was never sonicated, so it should read near zero;
* the front of the brain faces the image left in every section;
* the leave-one-out GFP fit (``results/histology/mouse01_slots.csv``) is validated
  section by section, and every Mouse 1 number in ``docs/`` comes from it.

This script places ROIs on the same four sections -- s2, s4, s5, s6; s3 is half
a section -- the way later animals are placed, and scores them against that
reference: finetuned shapes (``data/roi_locations.csv``) where a section has them,
otherwise the template from its clicks (``data/section_orientation.csv``). Mouse 2 and later animals can only be
placed by template, so the template earns trust here or nowhere.

Checks (thresholds are provisional, set before running):

B1  control      T3 coverage <= 0.05 in every section
B2  placement    median ROI-centre offset from the fit <= 0.5 mm (half the ROI
                 radius), and none > 1.0 mm
B3  delivery     each target's mean coverage within 0.10 of the fit's
B4  ranking      Spearman rho of the six target means vs the fit >= 0.9
B5  mirroring    the notch clicks reproduce the mirroring known from target 3

Writes ``results/histology/mouse01_benchmark.csv`` (one row per section x target).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from fus import histology, orientation, rois

ROOT = Path(__file__).resolve().parents[1]
SECTIONS = ["section_s2", "section_s4", "section_s5", "section_s6"]
REFERENCE = ROOT / "results/histology/mouse01_slots.csv"
OUT = ROOT / "results/histology/mouse01_benchmark.csv"
DS4 = ROOT / "data/processed/histology/Mouse_01/ds4"


def main() -> int:
    sheet = orientation.read_sheet()
    clicks = {s: sheet.get((1, s)) for s in SECTIONS}
    missing = [s for s, a in clicks.items() if a is None or not a.complete]
    if missing:
        sys.exit(f"click centre, front and notch side on Mouse 1 {', '.join(missing)} "
                 "(notebooks/ingest_new_histology.ipynb with MOUSE = 1)")
    try:
        plus_x = orientation.plus_x_side(sheet)
    except ValueError as e:
        plus_x = None
        b5_detail = str(e)
    else:
        b5_detail = f"targets 1-3 on the animal's {plus_x}"

    ref = pd.read_csv(REFERENCE)
    ref = ref[(ref.channel == "TRITC") & ref.section.isin(SECTIONS)].set_index(["section", "slot"])

    finetuned = rois.read_locations()
    rows = []
    for name in SECTIONS:
        a = clicks[name]
        mirrored = orientation.KNOWN_MIRRORED[(1, name)]
        sec = histology.load_section(DS4 / f"{name}.ome.tif")
        fit = histology.template_fit(a.centre, a.angle_deg, sec.pixel_um)
        hand = finetuned.get((1, name))
        res = histology.analyse_section(sec, template=fit, rois=hand)
        cov = {r["slot"]: r["coverage"] for r in res.rows if r["channel"] == "TRITC"}
        pts = fit.points() if hand is None else {roi.slot: roi.centre for roi in hand.values()}
        for slot in sorted(pts):
            r = ref.loc[(name, slot)]
            rows.append({
                "section": name,
                "target": histology.slot_to_target(slot, mirrored),
                "offset_mm": np.hypot(pts[slot][0] - r.row, pts[slot][1] - r.col) * sec.pixel_um / 1000,
                "coverage_fit": r.coverage,
                "coverage_template": cov[slot],
                "placement": "template" if hand is None else "finetuned",
                "click_angle_deg": a.angle_deg,
                "fit_angle_deg": r.fit_angle_deg,
            })
        print(f"{name}: {'finetuned' if hand else 'template'}; click front {a.angle_deg:5.1f}°, "
              f"fit {ref.loc[(name, 1)].fit_angle_deg:5.1f}°")
        del sec, res
    d = pd.DataFrame(rows)
    d.round(4).to_csv(OUT, index=False)

    means = d.groupby("target")[["coverage_fit", "coverage_template"]].mean()
    table = means.assign(
        diff=means.coverage_template - means.coverage_fit,
        offset_mm=d.groupby("target").offset_mm.mean(),
    ).round(2)
    print("\nPer target (mean of 4 sections, TRITC):")
    print(table.to_string())

    t3 = d[d.target == 3].coverage_template
    rho = stats.spearmanr(means.coverage_fit, means.coverage_template).statistic
    mirroring_ok = plus_x is not None and all(
        orientation.is_mirrored(clicks[s], plus_x) == orientation.KNOWN_MIRRORED[(1, s)]
        for s in SECTIONS)
    checks = [
        ("B1 control", (t3 <= 0.05).all(), f"T3 max {t3.max():.3f}"),
        ("B2 placement", d.offset_mm.median() <= 0.5 and d.offset_mm.max() <= 1.0,
         f"median {d.offset_mm.median():.2f} mm, max {d.offset_mm.max():.2f} mm"),
        ("B3 delivery", (table["diff"].abs() <= 0.10).all(),
         f"largest target difference {table['diff'].abs().max():.2f}"),
        ("B4 ranking", rho >= 0.9, f"Spearman {rho:.2f}"),
        ("B5 mirroring", mirroring_ok, b5_detail),
    ]
    print()
    for name, ok, detail in checks:
        print(f"{'PASS' if ok else 'FAIL'}  {name:13s} {detail}")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0 if all(ok for _, ok, _ in checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
