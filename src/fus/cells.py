"""Pipeline B: target ROI -> segmented NeuN cells -> GFP-tagged cells.

Pipeline A (`fus.histology`) measures each target ROI by its pixels: the fraction
that is GFP+ (coverage) and the mean background-subtracted intensity. Pipeline B
measures the same ROI by its cells:

1. **ROI** -- the target's shape exactly as pipeline A uses it: the finetuned
   ellipse in ``data/roi_locations.csv`` if saved, else the template placed from
   the orientation clicks.
2. **Export** the six ROI bounding boxes (plus a margin) at full resolution
   (0.325 um/px) straight from the ``.vsi``, in one QuPath launch
   (``scripts/qupath/export_regions.groovy``), uncompressed -- they are temporary.
3. **Segment** NeuN (CY5) with StarDist 2D_versatile_fluo
   (``scripts/stardist_segment.py``), one process per ROI, each started as soon as
   its crop is written. StarDist's network runs on the Apple GPU
   (``.venv-stardist-gpu``: TensorFlow 2.18 + tensorflow-metal); its slow step, merging
   candidate outlines, is single-core, so the six ROIs run side by side. GPU and CPU
   give the same cells (5,717 of 5,717 matched on one crop).
4. **Tag** each cell whose centroid is inside the ROI: its mean raw anti-GFP
   (TRITC) divided by the section's **background** -- the median anti-GFP over the
   section's tissue -- GFP+ above `THRESHOLD` (fold over background).

`THRESHOLD` = 4.05x background is the log-Otsu split of that score on the pilot's
3 x 3 mm crop around Mouse 2 slide04_s3 T3 (``segmentation_alpha/results/
mouse02_slide04_s3_T3_3000um_gfp_threshold.json`` -> ``fold``; docs/gfp-cell-tagging.md).
Dividing by each section's own background is what lets one threshold work across
animals: Mouse 1's un-transduced tissue is ~5x brighter than Mouse 2's, and the
earlier percentile scale (0.0575) tagged 5-14% of its no-FUS control. It was set on
one section of one animal and is not yet hand-validated.

Outputs:

  results/histology/cells/mouseNN_cell_rois.csv          one row per section x target
  data/processed/histology/Mouse_NN/cells/<section>_cells.csv    one row per cell (git-ignored)
  data/processed/histology/Mouse_NN/cells/<section>_T<k>_labels.tif

Command line::

    .venv/bin/python -m fus.cells --mouse 2 --section slide04_s3
    .venv/bin/python -m fus.cells --mouse 2 --all
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

from fus import histology, orientation, rois

ROOT = histology.REPO_ROOT
THRESHOLD = 4.05          # fold over the section's background
WORKERS = 6               # StarDist processes at once (one per ROI)
MARGIN_UM = 30.0           # box margin so cells on the ROI edge are segmented whole
DS4 = 4                    # the ds4 exports are 4x downsampled from full resolution
NEUN, GFP = 3, 2           # channel order DAPI, FITC, TRITC, CY5
_GPU_ENV = ROOT / "segmentation_alpha/.venv-stardist-gpu/bin/python"
STARDIST_PYTHON = Path(os.environ.get(
    "FUS_STARDIST_PYTHON",
    _GPU_ENV if _GPU_ENV.exists() else ROOT / "segmentation_alpha/.venv-stardist/bin/python"))
SEGMENT_SCRIPT = ROOT / "scripts/stardist_segment.py"
EXPORT_SCRIPT = ROOT / "scripts/qupath/export_regions.groovy"
RESULTS = ROOT / "results/histology/cells"


def work_dir(mouse: int) -> Path:
    return ROOT / f"data/processed/histology/Mouse_{mouse:02d}/cells"


def ds4_path(mouse: int, section: str) -> Path:
    return ROOT / f"data/processed/histology/Mouse_{mouse:02d}/ds4/{section}.ome.tif"


def source(mouse: int, section: str) -> tuple[Path, int]:
    """The .vsi and series a section was exported from (section_s4, slide01_s3, ...)."""
    prefix, series = re.fullmatch(r"(.+)_s(\d+)", section).groups()
    raw = ROOT / f"data/histology/Mouse_{mouse:02d}"
    if prefix == "section":
        return raw / "Image.vsi", int(series)
    return next(raw.glob(f"*Slide_{int(prefix.removeprefix('slide')):02d}.vsi")), int(series)


def section_rois(mouse: int, section: str) -> tuple[dict[int, rois.Roi], str]:
    """``{target: Roi}`` in ds4 pixels, as pipeline A measures them, and their origin."""
    saved = rois.read_locations().get((mouse, section))
    if saved:
        return saved, "finetuned"
    sheet = orientation.read_sheet()
    a = sheet[(mouse, section)]
    mirrored = orientation.KNOWN_MIRRORED.get((mouse, section))
    if mirrored is None:
        mirrored = orientation.is_mirrored(a, orientation.plus_x_side(sheet))
    return rois.from_template(a, rois.export_pixel_um(ds4_path(mouse, section)), mirrored), "template"


def section_background(mouse: int, section: str) -> tuple[float, np.ndarray]:
    """Median anti-GFP over the section's tissue (4x export), and the tissue mask."""
    sec = histology.load_section(ds4_path(mouse, section))
    tissue = histology.tissue_mask(sec)
    return float(np.median(sec.channels["TRITC"][tissue])), tissue


def box(roi: rois.Roi, full_um: float) -> tuple[int, int, int, int]:
    """Full-resolution x, y, w, h around a ds4 ROI, with `MARGIN_UM`."""
    ext = max(roi.semi_a, roi.semi_b) * DS4 + MARGIN_UM / full_um
    r, c = roi.centre[0] * DS4, roi.centre[1] * DS4
    x, y = int(np.floor(c - ext)), int(np.floor(r - ext))
    return max(x, 0), max(y, 0), int(np.ceil(2 * ext)), int(np.ceil(2 * ext))


def export_and_segment(vsi: Path, series: int, jobs: dict, verbose: bool = True) -> None:
    """Export every ROI box in one QuPath launch; start StarDist on each as it lands.

    ``jobs`` maps target -> (xywh, crop path, labels path).
    """
    if not STARDIST_PYTHON.exists():
        raise FileNotFoundError(f"StarDist environment not found at {STARDIST_PYTHON}; "
                                "see docs/gfp-cell-tagging.md or set $FUS_STARDIST_PYTHON")
    by_crop = {str(crop.resolve()): (t, labels) for t, (_, crop, labels) in jobs.items()}
    specs = []
    for xywh, crop, _ in jobs.values():
        specs += ["-a", ",".join(map(str, xywh)) + "," + str(crop.resolve())]
    exporter = subprocess.Popen(
        [histology._qupath(), "script", "-a", str(vsi.resolve()), "-a", str(series),
         "-a", "UNCOMPRESSED", *specs, str(EXPORT_SCRIPT)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    running = []

    def start(crop: str):
        while len([p for p in running if p.poll() is None]) >= WORKERS:
            time.sleep(0.5)
        t, labels = by_crop[crop]
        if verbose:
            print(f"  T{t} exported -> segmenting", flush=True)
        log = open(Path(str(labels)).with_suffix(".log"), "w")  # not a pipe: TF is chatty
        running.append(subprocess.Popen(
            [str(STARDIST_PYTHON), str(SEGMENT_SCRIPT), "--channel", str(NEUN), crop, str(labels)],
            stdout=log, stderr=subprocess.STDOUT))
        running[-1].log = log

    for line in exporter.stdout:
        if line.startswith("WROTE "):
            start(line[6:].strip())
    if exporter.wait() != 0:
        raise RuntimeError(f"QuPath region export failed for {vsi.name} series {series}")
    for p in running:
        code = p.wait()
        p.log.close()
        if code != 0:
            raise RuntimeError(f"StarDist failed; see {p.log.name}")
    missing = [t for t, (_, _, labels) in jobs.items() if not labels.exists()]
    if missing:
        raise RuntimeError(f"no labels for targets {missing}")


def tag(crop: Path, labels_path: Path, roi_full: rois.Roi, background: float,
        threshold: float = THRESHOLD) -> pd.DataFrame:
    """One row per cell with its centroid inside ``roi_full`` (crop pixel coordinates).

    ``gfp_fold`` = the cell's mean raw anti-GFP / ``background``.
    """
    with tifffile.TiffFile(crop) as tf:
        gfp = tf.series[0].levels[0].asarray()[GFP].astype(np.float32) / background
    labels = tifffile.imread(labels_path)
    ea, eb = roi_full.axes()
    rows = []
    for lab, sl in enumerate(ndi.find_objects(labels), start=1):
        if sl is None:
            continue
        mask = labels[sl] == lab
        r, c = ndi.center_of_mass(mask)
        r, c = r + sl[0].start, c + sl[1].start
        dr, dc = r - roi_full.centre[0], c - roi_full.centre[1]
        u = (dr * ea[0] + dc * ea[1]) / roi_full.semi_a
        w = (dr * eb[0] + dc * eb[1]) / roi_full.semi_b
        if u**2 + w**2 > 1:
            continue
        rows.append({"label": lab, "row": r, "col": c, "area_px": int(mask.sum()),
                     "gfp_fold": float(gfp[sl][mask].mean())})
    cells = pd.DataFrame(rows, columns=["label", "row", "col", "area_px", "gfp_fold"])
    cells["gfp_positive"] = cells.gfp_fold > threshold
    return cells


def measure_section(mouse: int, section: str, threshold: float = THRESHOLD,
                    keep_crops: bool = False, verbose: bool = True) -> pd.DataFrame:
    """Pipeline B for every target of one section; returns and saves one row per target."""
    t0 = time.time()
    vsi, series = source(mouse, section)
    shapes, origin = section_rois(mouse, section)
    background, tissue = section_background(mouse, section)
    ds4_um = rois.export_pixel_um(ds4_path(mouse, section))
    full_um = ds4_um / DS4
    work = work_dir(mouse)
    work.mkdir(parents=True, exist_ok=True)

    jobs = {}
    for t, roi in sorted(shapes.items()):
        crop, labels = work / f"{section}_T{t}.ome.tif", work / f"{section}_T{t}_labels.tif"
        labels.unlink(missing_ok=True)
        jobs[t] = (box(roi, full_um), crop, labels)
    if verbose:
        print(f"{section}: exporting {len(jobs)} ROIs and segmenting", flush=True)
    export_and_segment(vsi, series, jobs, verbose)

    summary, all_cells = [], []
    for t, (xywh, crop, lab_path) in jobs.items():
        roi = shapes[t]
        x, y = xywh[0], xywh[1]
        roi_full = rois.Roi(t, roi.slot, (roi.centre[0] * DS4 - y, roi.centre[1] * DS4 - x),
                            roi.semi_a * DS4, roi.semi_b * DS4, roi.angle_deg)
        cells = tag(crop, lab_path, roi_full, background, threshold)
        cells["row"] += y
        cells["col"] += x
        cells["area_um2"] = cells.pop("area_px") * full_um**2
        cells.insert(0, "target", t)
        all_cells.append(cells)
        box_, inside = roi.mask(tissue.shape)
        tissue_mm2 = float((tissue[box_] & inside).sum() * (ds4_um / 1000) ** 2)
        n, pos = len(cells), int(cells.gfp_positive.sum())
        summary.append({
            "mouse": mouse, "section": section, "target": t, "slot": roi.slot, "roi_source": origin,
            "roi_tissue_mm2": round(tissue_mm2, 3), "n_cells": n,
            "cells_per_mm2": round(n / tissue_mm2, 1) if tissue_mm2 else np.nan,
            "n_gfp_pos": pos, "fraction_gfp_pos": round(pos / n, 4) if n else np.nan,
            "median_cell_fold": round(float(cells.gfp_fold.median()), 3) if n else np.nan,
            "threshold_fold": threshold, "background": round(background, 1),
        })
        if not keep_crops:
            crop.unlink(missing_ok=True)

    cells = pd.concat(all_cells, ignore_index=True)
    cells.insert(0, "section", section)
    cells.round(4).to_csv(work / f"{section}_cells.csv", index=False)
    out = pd.DataFrame(summary)
    save_summary(mouse, out)
    if verbose:
        print(f"{section}: done in {time.time() - t0:.0f} s", flush=True)
    return out


def save_summary(mouse: int, rows: pd.DataFrame) -> Path:
    """Merge one section's rows into results/histology/cells/mouseNN_cell_rois.csv."""
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f"mouse{mouse:02d}_cell_rois.csv"
    if path.exists():
        old = pd.read_csv(path)
        old = old[~old.section.isin(rows.section.unique())]
        rows = pd.concat([old, rows], ignore_index=True)
    rows.sort_values(["section", "target"]).to_csv(path, index=False)
    return path


def sections(mouse: int) -> list[str]:
    """Sections of a mouse that are clicked and not excluded."""
    return sorted(a.section for (m, _), a in orientation.read_sheet().items()
                  if m == mouse and a.complete and not a.exclude)


def main(argv: list[str] | None = None) -> int:
    import argparse

    p = argparse.ArgumentParser(prog="python -m fus.cells", description=__doc__.split("\n")[0])
    p.add_argument("--mouse", type=int, required=True)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--section", nargs="+")
    g.add_argument("--all", action="store_true", help="every clicked, non-excluded section")
    p.add_argument("--threshold", type=float, default=THRESHOLD, help="fold over background")
    p.add_argument("--keep-crops", action="store_true")
    args = p.parse_args(argv)
    for section in (sections(args.mouse) if args.all else args.section):
        print(measure_section(args.mouse, section, args.threshold, args.keep_crops)
              .drop(columns=["threshold_fold"]).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
