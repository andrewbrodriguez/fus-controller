"""GFP+ calling for StarDist NeuN cells: mean normalised GFP per cell, one threshold.

    .venv/bin/python segmentation_alpha/gfp_tagging.py [crop.ome.tif]

1. **Normalise** the anti-GFP stain (TRITC) over the *whole section*: 0 = its 1st
   percentile over tissue, 1 = its 99.9th (from the 4x export, the same tissue
   mask the coverage pipeline uses). Linear, not clipped -- a bright cell can
   exceed 1. Using the section, not the crop, keeps the scale the same whatever
   crop is analysed.
2. **Score** each cell (StarDist on NeuN; run_stardist.py makes the labels) by the
   mean normalised GFP over its pixels.
3. **Threshold**, first pass: Otsu's split of the cell scores on a **log** scale.
   Scores span two decades, and on log10 they are bimodal: a background peak near
   0.014 and a plume hump near 0.2-0.4. Log-Otsu lands in the valley between them.
   Two alternatives are reported alongside:
   * reference -- the `REFERENCE_PERCENTILE` of random cell shapes in cell-free
     tissue > `OFF_TARGET_MM` from every target ("brighter than 99% of
     un-sonicated tissue"); permissive, it takes in the plume's dim halo;
   * linear Otsu -- strict, it splits the plume population itself.
   The first pass was going to be the reference; it was switched to log-Otsu after
   seeing the histogram (2026-09-30). The notebook can apply any threshold.

**Fold over background** (2026-09-30, what pipeline B uses): each cell's mean raw
anti-GFP divided by its section's background, the median anti-GFP over the
section's tissue (4x export). Mouse 1's un-transduced tissue is ~5x brighter than
Mouse 2's, so a threshold on the percentile scale above didn't transfer between
animals; dividing by each section's own background does. Its log-Otsu split is
`threshold_fold`, the constant in `fus.cells`.

Writes results/<crop>_gfp_cells.csv, _gfp_reference.csv and _gfp_threshold.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu as otsu

from fus import histology, rois

HERE = Path(__file__).resolve().parent
CROP = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE / "mouse02_slide04_s3_T3_3000um.ome.tif"
STEM = CROP.name.split(".")[0]
LABELS = HERE / "results" / f"{STEM}_stardist_labels.tif"
GFP = 2  # TRITC (anti-GFP stain); order DAPI, FITC, TRITC, CY5

NORM_PERCENTILES = (1.0, 99.9)
OFF_TARGET_MM = 1.5
CELL_GAP_UM = 1.0        # reference patches keep this far from any cell
N_REFERENCE = 4000
REFERENCE_PERCENTILE = 99.0
SEED = 0


def crop_info() -> dict:
    return json.loads(CROP.with_name(STEM + ".json").read_text())


def section_scale(info: dict) -> tuple[float, float, np.ndarray, float]:
    """(lo, hi) of TRITC over the section's tissue, plus its tissue mask and px size."""
    ds4 = histology.REPO_ROOT / info["source"].replace("/full/", "/ds4/")
    sec = histology.load_section(ds4)
    tissue = histology.tissue_mask(sec)
    lo, hi = np.percentile(sec.channels["TRITC"][tissue], NORM_PERCENTILES)
    section_scale.background = float(np.median(sec.channels["TRITC"][tissue]))
    return float(lo), float(hi), tissue, rois.export_pixel_um(ds4)


def crop_mask(section_mask: np.ndarray, ds4_um: float, info: dict) -> np.ndarray:
    """The section tissue mask cut to the crop, at the crop's resolution."""
    f = ds4_um / info["pixel_um"]
    (r0, r1), (c0, c1) = info["rows"], info["cols"]
    return section_mask[np.ix_((np.arange(r0, r1) / f).astype(int), (np.arange(c0, c1) / f).astype(int))]


def target_centres(info: dict, ds4_um: float) -> dict[int, tuple[float, float]]:
    """ROI centre of every target of this section, in crop pixels (may lie outside)."""
    f = ds4_um / info["pixel_um"]
    shapes = rois.read_locations()[(info["mouse"], info["section"])]
    return {t: (roi.centre[0] * f - info["rows"][0], roi.centre[1] * f - info["cols"][0])
            for t, roi in shapes.items()}


def place(rng, shapes, allowed, n):
    """Random real cell shapes at positions where every pixel is ``allowed``."""
    out, tries = [], 0
    h, w = allowed.shape
    while len(out) < n and tries < 200 * n:
        tries += 1
        m = shapes[rng.integers(len(shapes))]
        r0 = int(rng.integers(0, h - m.shape[0]))
        c0 = int(rng.integers(0, w - m.shape[1]))
        if allowed[r0:r0 + m.shape[0], c0:c0 + m.shape[1]][m].all():
            out.append((m, r0, c0))
    return out


def main() -> None:
    info = crop_info()
    with tifffile.TiffFile(CROP) as tf:
        image = tf.asarray()
    px_um = info["pixel_um"]
    labels = tifffile.imread(LABELS)

    lo, hi, section_tissue, ds4_um = section_scale(info)
    background = section_scale.background
    gfp = (image[GFP].astype(np.float32) - lo) / (hi - lo)
    tissue = crop_mask(section_tissue, ds4_um, info)

    objects = ndi.find_objects(labels)
    rows = []
    for lab, sl in enumerate(objects, start=1):
        if sl is None:
            continue
        mask = labels[sl] == lab
        v = gfp[sl][mask]
        raw = float(image[GFP][sl][mask].mean())
        r, c = ndi.center_of_mass(mask)
        rows.append({"label": lab, "row": sl[0].start + r, "col": sl[1].start + c,
                     "area_um2": float(mask.sum() * px_um**2),
                     "gfp_mean": float(v.mean()), "gfp_sd": float(v.std()),
                     "gfp_fold": raw / background})
    cells = pd.DataFrame(rows)

    centres = target_centres(info, ds4_um)
    cells["dist_T_mm"] = np.min([np.hypot(cells.row - r, cells.col - c) * px_um / 1000
                                 for r, c in centres.values()], axis=0)

    # Reference: cell-shaped patches in cell-free tissue far from every target.
    yy, xx = np.ogrid[:gfp.shape[0], :gfp.shape[1]]
    near_target = np.zeros(gfp.shape, bool)
    for r, c in centres.values():
        near_target |= (yy - r) ** 2 + (xx - c) ** 2 <= (OFF_TARGET_MM * 1000 / px_um) ** 2
    cell_free = ndi.distance_transform_edt(labels == 0) * px_um > CELL_GAP_UM
    allowed = tissue & cell_free & ~near_target
    shapes = [labels[sl] == i for i, sl in enumerate(objects, start=1) if sl is not None]
    reference = pd.DataFrame()
    threshold_ref = None
    if allowed.mean() > 0.01:
        rng = np.random.default_rng(SEED)
        reference = pd.DataFrame([
            {"gfp_mean": float(gfp[r0:r0 + m.shape[0], c0:c0 + m.shape[1]][m].mean()),
             "area_um2": float(m.sum() * px_um**2)}
            for m, r0, c0 in place(rng, shapes, allowed, N_REFERENCE)])
        threshold_ref = float(np.percentile(reference.gfp_mean, REFERENCE_PERCENTILE))
    threshold_otsu = float(otsu(cells.gfp_mean.to_numpy()))
    threshold_log_otsu = float(10 ** otsu(np.log10(np.clip(cells.gfp_mean.to_numpy(), 1e-4, None))))
    threshold = threshold_log_otsu
    threshold_fold = float(10 ** otsu(np.log10(np.clip(cells.gfp_fold.to_numpy(), 1e-3, None))))
    cells["gfp_positive"] = cells.gfp_mean > threshold

    out = HERE / "results"
    cells.round(4).to_csv(out / f"{STEM}_gfp_cells.csv", index=False)
    if len(reference):
        reference.round(4).to_csv(out / f"{STEM}_gfp_reference.csv", index=False)
    summary = {
        "crop": CROP.name, "labels": LABELS.name,
        "score": "mean over the cell's pixels of normalised TRITC",
        "normalisation": {"over": "whole-section tissue (4x export)", "percentiles": NORM_PERCENTILES,
                          "lo": round(lo, 1), "hi": round(hi, 1)},
        "threshold": round(threshold, 4),
        "threshold_from": "log_otsu",
        "threshold_log_otsu": round(threshold_log_otsu, 4),
        "threshold_reference": None if threshold_ref is None else round(threshold_ref, 4),
        "reference": None if threshold_ref is None else {
            "patches": len(reference), "percentile": REFERENCE_PERCENTILE,
            "where": f"cell-free tissue > {OFF_TARGET_MM} mm from every target",
            "fraction_of_crop": round(float(allowed.mean()), 3)},
        "threshold_otsu": round(threshold_otsu, 4),
        "fold": {"background": round(background, 1),
                 "background_is": "median anti-GFP over the section's tissue (4x export)",
                 "threshold_fold": round(threshold_fold, 3),
                 "gfp_positive": int((cells.gfp_fold > threshold_fold).sum()),
                 "fraction_positive": round(float((cells.gfp_fold > threshold_fold).mean()), 3),
                 "agreement_with_percentile_scale": round(float(
                     ((cells.gfp_fold > threshold_fold) == cells.gfp_positive).mean()), 4)},
        "cells": int(len(cells)), "gfp_positive": int(cells.gfp_positive.sum()),
        "fraction_positive": round(float(cells.gfp_positive.mean()), 3),
    }
    (out / f"{STEM}_gfp_threshold.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Fraction GFP+ by distance to the nearest target (T ROI radius = 0.87 mm), per threshold.
    bands = pd.cut(cells.dist_T_mm, [0, 0.25, 0.5, 0.87, 1.25, 1.5, 1.75, 2.25])
    profile = pd.DataFrame({"cells": cells.groupby(bands, observed=True).size()})
    for name, t in (("log_otsu", threshold_log_otsu), ("reference", threshold_ref), ("otsu", threshold_otsu)):
        if t is not None:
            profile[f"positive_{name}"] = (cells.gfp_mean > t).groupby(bands, observed=True).mean().round(4)
    profile["positive_fold"] = (cells.gfp_fold > threshold_fold).groupby(bands, observed=True).mean().round(4)
    profile.index = profile.index.astype(str)
    profile.rename_axis("dist_mm").to_csv(out / f"{STEM}_gfp_profile.csv")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
