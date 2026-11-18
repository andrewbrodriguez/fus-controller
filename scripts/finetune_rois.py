"""Finetune T1-T6 by hand on every kept section, in one napari window.

    .venv/bin/python scripts/finetune_rois.py              # every kept section of every mouse
    .venv/bin/python scripts/finetune_rois.py --mouse 5 6  # only these animals

"Kept" means the section's orientation clicks are complete and it isn't excluded in
``data/section_orientation.csv``. Mice come from ``MICE`` in ``slices_dose_delivery.py``.

The window shows one section at a time, front up and notch side left, with NeuN in grey and
the GFP stain in green. T1-T6 are six yellow shapes. They start from the saved locations in
``data/roi_locations.csv``, or from the template placed by the clicks if a section has none.
With the targets layer selected and in select mode (S):

- click a shape and drag it to **move** it
- drag a corner handle to **resize** (hold Shift to keep it round)
- drag the handle above the box to rotate
- hold Space and drag to pan

**Enter** (Save + next) saves the section and moves on, **B** goes back, **Reset section to
template** starts a section over, and **Done** saves and closes. Every move between sections
writes ``data/roi_locations.csv``, so closing the window loses at most the section on screen.
The file is tracked in git: ``git diff data/roi_locations.csv`` shows what changed. When the
window closes, this script prints how far each target moved.

Two cautions:

- Resizing changes the measurement. Coverage is the GFP+ fraction *inside the shape*, so a
  smaller shape on a bright plume reads higher.
- Placing shapes by eye on the GFP can pull a target onto its own signal. Anatomy is the
  safer guide.

New shapes change both delivery measures, so every downstream number has to be re-run.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from fus import histology, orientation, rois

sys.path.insert(0, str(Path(__file__).resolve().parent))
from slices_dose_delivery import MICE  # noqa: E402  (the one list of imaged mice)

ROOT = histology.REPO_ROOT


def kept_sections(mice) -> tuple[list[Path], list[str]]:
    """Exported sections with complete clicks and no exclusion; and kept ones not on disk."""
    sheet = orientation.read_sheet()
    paths, missing = [], []
    for m in mice:
        on_disk = {orientation.section_name(p): p for p in
                   sorted((ROOT / f"data/processed/histology/Mouse_{m:02d}/ds4").glob("*.ome.tif"))}
        for (mouse, name), a in sorted(sheet.items()):
            if mouse != m or a.exclude or not a.complete:
                continue
            if name in on_disk:
                paths.append(on_disk[name])
            else:
                missing.append(f"Mouse {m} {name}")
    return paths, missing


def report(paths, before: dict, after: dict) -> None:
    """How far each target moved and how its size changed, in mm on the tissue."""
    changed, untouched, unsaved = [], 0, []
    for path in paths:
        k = (orientation.mouse_of(path), orientation.section_name(path))
        new = after.get(k)
        if new is None:
            unsaved.append(k)
            continue
        old = before.get(k)
        mm = rois.export_pixel_um(path) / 1000
        if old is None:
            sizes = ", ".join(f"T{t} {2 * r.radius * mm:.2f}" for t, r in sorted(new.items()))
            changed.append((k, [f"first save (was the template); mm across: {sizes}"]))
            continue
        lines = []
        for t, r in sorted(new.items()):
            o = old[t]
            moved = float(np.hypot(*np.subtract(r.centre, o.centre))) * mm
            d_old, d_new = 2 * o.radius * mm, 2 * r.radius * mm
            if moved >= 0.01 or abs(d_new - d_old) >= 0.01 or r.angle_deg != o.angle_deg:
                lines.append(f"T{t} moved {moved:.2f} mm, {d_old:.2f} -> {d_new:.2f} mm across")
        if lines:
            changed.append((k, lines))
        else:
            untouched += 1

    print(f"\n{len(changed)} of {len(paths)} sections changed, {untouched} unchanged, "
          f"{len(unsaved)} never saved")
    for (m, name), lines in changed:
        print(f"  Mouse {m} {name}")
        for line in lines:
            print(f"    {line}")
    for m, name in unsaved:
        print(f"  Mouse {m} {name}: no saved shapes (measured with the template)")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="scripts/finetune_rois.py",
                                description="Move and resize T1-T6 on every kept section.")
    p.add_argument("--mouse", type=int, nargs="+", default=list(MICE),
                   help=f"mice to include (default: all imaged, {', '.join(map(str, MICE))})")
    args = p.parse_args(argv)

    paths, missing = kept_sections(args.mouse)
    for s in missing:
        print(f"{s}: kept in the sheet but not exported; skipping")
    if not paths:
        print("no kept sections to finetune")
        return 1
    print(f"{len(paths)} sections: " + ", ".join(
        f"Mouse {m} ({sum(orientation.mouse_of(q) == m for q in paths)})" for m in args.mouse))

    import napari

    before = rois.read_locations()
    rois.finetune(paths, [orientation.mouse_of(q) for q in paths])
    napari.run()
    report(paths, before, rois.read_locations())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
