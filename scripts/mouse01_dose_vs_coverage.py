"""Mouse 1: acoustic dose vs GFP coverage, under explicit assumptions.

    PYTHONPATH=src python scripts/mouse01_dose_vs_coverage.py

Reads the per-slot histology table (``python -m fus.histology measure``) and
the seven Mouse_Cntr_01 recordings and joins them per target. Writes
``results/histology/mouse01_dose_vs_coverage.csv`` (both dose measures x both GFP
channels). Mouse 1 also appears in the all-mice figures of
``scripts/slices_dose_delivery.py``; there is no per-mouse figure.

This is exploratory. Every join below rests on an assumption that has not been
confirmed with the lab:

A1  Section orientation: s2 and s4 as scanned, s5 and s6 mirrored. Confirmed
    9/18: sections carry a notch (bottom right on Mouse 1) and target 3 was
    never sonicated, so a flipped section is identifiable.
A2  The empty site is target 3 -- confirmed by the lab.
A3  Coverage per target is the mean over s2, s4, s5, s6; bars span min-max.
A4  Recording -> position, **confirmed** by the lab (9/18) and by the
    acquisition timestamps in each file: the first three sonications all went
    to position 1 (Target1, Target2, Target3), then positions 4, 5, 6 in turn,
    then a 7th back at position 2 (Target2_Repeat). Position 3 was never
    sonicated. For every other mouse, TargetN was fired at position N.
A5  Repeat sonications at one position add: doses are summed.
A6  The no-FUS control has zero dose.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from fus import archive
from fus import extract
from fus.histology import slot_to_target

ROOT = Path(__file__).resolve().parents[1]
SLOTS = ROOT / "results" / "histology" / "mouse01_slots.csv"
OUT = ROOT / "results" / "histology" / "mouse01_dose_vs_coverage"

MIRRORED = {"section_s2": False, "section_s4": False, "section_s5": True, "section_s6": True}

MAPPINGS = {
    "confirmed": {1: ["Target1", "Target2", "Target3"], 2: ["Target2_Repeat"], 3: [],
                  4: ["Target4"], 5: ["Target5"], 6: ["Target6"]},
}



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
    for m in (m for m in archive.recordings() if m.mouse == 1):
        with archive.local_path(m) as path:
            out[m.recording] = extract(path).to_dict()
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


def main() -> None:
    sections = section_coverage()
    df = join(recordings(), coverage_by_target(sections))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.round(4).to_csv(OUT.with_suffix(".csv"), index=False)
    print(df[["mapping", "target", "recordings", "cum_2nd_harmonic", "bursts_x_goal",
              "TRITC_mean", "FITC_mean"]].round(3).to_string(index=False))
    print(f"\nwrote {OUT.with_suffix('.csv')}")


if __name__ == "__main__":
    main()
