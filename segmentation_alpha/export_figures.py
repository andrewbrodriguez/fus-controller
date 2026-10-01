"""Two figures of GFP tagging on the T3 ROI edge, rendered as the napari viewer shows them.

    .venv/bin/python segmentation_alpha/export_figures.py

A 400 x 400 um window of the 3 mm crop (Mouse 2 slide04_s3), centred on the T3 ROI
edge on the side away from the other targets, so the plume's inside and outside
are both in frame. Writes to docs/figures/:

  gfp-tagging-stains.png   NeuN (grey, 0.4 opacity) + anti-GFP stain (green), additive,
                           GFP on the section's normalised 0-1 scale
  gfp-tagging-cells.png    the same window with the stains off: StarDist cells outlined
                           yellow (GFP+) or magenta (GFP-) at the first-pass threshold

napari can't screenshot a hidden window, so this reproduces its rendering: additive
blending of each layer's colormap x opacity, the same contrast limits as
gfp_tagging.ipynb, and 2 px label contours.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
STEM = "mouse02_slide04_s3_T3_3000um"
OUT = ROOT / "docs/figures"
WINDOW_UM = 400.0
EDGE_MM = 0.87          # T3 ROI radius
DIRECTION = (0, -1)     # toward the image left: no other target that side
CY5, TRITC = 3, 2
YELLOW, MAGENTA = (1.0, 0.92, 0.0), (1.0, 0.0, 1.0)


def window(info):
    px = info["pixel_um"]
    cr, cc = np.array(info["roi_centre_fullres_px"]) - [info["rows"][0], info["cols"][0]]
    r = cr + DIRECTION[0] * EDGE_MM * 1000 / px
    c = cc + DIRECTION[1] * EDGE_MM * 1000 / px
    half = int(round(WINDOW_UM / px / 2))
    return slice(int(r) - half, int(r) + half), slice(int(c) - half, int(c) + half)


def scale(x, lo, hi):
    return np.clip((x.astype(np.float32) - lo) / (hi - lo), 0, 1)


def contours(labels, width=2):
    """Pixels of each cell within ``width`` px of a different label (napari contour)."""
    size = 2 * width + 1
    inside = labels > 0
    return inside & ((ndi.maximum_filter(labels, size) != labels) | (ndi.minimum_filter(labels, size) != labels))


def save(rgb, path, px_um, legend=None):
    h, w = rgb.shape[:2]
    fig = plt.figure(figsize=(w / 200, h / 200), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(rgb, interpolation="nearest")
    ax.set_axis_off()
    bar = 100 / px_um  # 100 um scale bar, bottom right
    x1, y = w - 0.05 * w, h - 0.06 * h
    ax.plot([x1 - bar, x1], [y, y], color="white", linewidth=4, solid_capstyle="butt")
    ax.text(x1 - bar / 2, y - 0.02 * h, "100 µm", color="white", ha="center", va="bottom",
            fontsize=11, family=["Helvetica Neue", "Arial", "DejaVu Sans"])
    if legend:
        for i, (label, colour) in enumerate(legend):
            yy = 0.05 * h + i * 0.055 * h
            ax.add_patch(plt.Rectangle((0.04 * w, yy - 0.015 * h), 0.03 * h, 0.03 * h,
                                       fill=False, edgecolor=colour, linewidth=2.5))
            ax.text(0.04 * w + 0.045 * h, yy, label, color="white", va="center", fontsize=11,
                    family=["Helvetica Neue", "Arial", "DejaVu Sans"])
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    info = json.loads((HERE / f"{STEM}.json").read_text())
    summary = json.loads((HERE / f"results/{STEM}_gfp_threshold.json").read_text())
    cells = pd.read_csv(HERE / f"results/{STEM}_gfp_cells.csv")
    image = tifffile.imread(HERE / f"{STEM}.ome.tif")
    labels_all = tifffile.imread(HERE / f"results/{STEM}_stardist_labels.tif")
    px = info["pixel_um"]
    rs, cs = window(info)

    # Stains: grey NeuN x 0.4 + green GFP, additive, clipped -- napari's blending.
    cy5_lo, cy5_hi = np.percentile(image[CY5][::8, ::8], [1, 99.8])
    lo, hi = summary["normalisation"]["lo"], summary["normalisation"]["hi"]
    neun = scale(image[CY5][rs, cs], cy5_lo, cy5_hi) * 0.4
    gfp = scale(image[TRITC][rs, cs], lo, hi)
    stains = np.clip(np.dstack([neun, neun + gfp, neun]), 0, 1)
    OUT.mkdir(parents=True, exist_ok=True)
    save(stains, OUT / "gfp-tagging-stains.png", px,
         legend=[("NeuN (CY5)", "#bdbdbd"), ("GFP stain (TRITC)", "#00ff00")])

    # Cells: outlines on black, coloured by the first-pass tag.
    labels = labels_all[rs, cs]
    positive = set(cells.label[cells.gfp_positive].astype(int))
    edge = contours(labels)
    is_pos = np.isin(labels, list(positive))
    rgb = np.zeros(labels.shape + (3,), np.float32)
    rgb[edge & is_pos] = YELLOW
    rgb[edge & ~is_pos] = MAGENTA
    save(rgb, OUT / "gfp-tagging-cells.png", px,
         legend=[("GFP+", YELLOW), ("GFP−", MAGENTA)])

    ids = np.unique(labels[labels > 0])
    in_win = cells[cells.label.isin(ids)]
    print(f"window {WINDOW_UM:.0f} um at the T3 ROI edge: {len(in_win)} cells (incl. partial), "
          f"{int(in_win.gfp_positive.sum())} GFP+ at threshold {summary['threshold']:.3f}")
    print(f"wrote {OUT / 'gfp-tagging-stains.png'}\nwrote {OUT / 'gfp-tagging-cells.png'}")


if __name__ == "__main__":
    main()
