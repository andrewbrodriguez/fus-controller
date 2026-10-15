"""Every clicked section, re-rendered into one common orientation per mouse.

    .venv/bin/python scripts/aligned_sections.py            # mice 1 and 2
    .venv/bin/python scripts/aligned_sections.py --mouse 2

The ``*_slots.png`` figures in ``results/histology/figures/`` show each section as it lay
on the slide, at any rotation and sometimes face down. Here each section is resampled
from its ds4 export (pyramid level 2, 5.2 um/px) using the person's clicks
(``orientation.canonical_affine``, the same transform as ``fus.reconstruct``). The centre
click goes to the origin and the front points up. Both mice are then shown as a **dorsal
view with the animal's left on the image left**. The canonical frame puts the notch side
left, and the notch is on the animal's right in Mouse 1 and the left from Mouse 2 on
(``orientation.notch_side``), so Mouse 1 is flipped left-right. Targets 1-3 then sit on the
same side for every animal.

Image: NeuN (CY5, log, grey) with anti-GFP (TRITC, green), each stretched per section, so
brightness is not comparable between sections; use the measured coverage for that. Circles
are the finetuned target ROIs from ``data/roi_locations.csv``, labelled with target and
pixel coverage (pipeline A).

Writes ``results/histology/aligned/``:

  mouseNN/<section>.png    one aligned section each
  mouseNN_grid.png         all of that mouse's sections at the same scale and crop, plus
                           every tissue outline overlaid to show how well they line up
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse

from fus import histology, orientation, reconstruct, rois

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/histology/aligned"
SLICES = ROOT / "results/histology/slices/slices_dose_delivery.csv"
LEVEL = 2            # pyramid level: 4x the export pixel (~5.2 um)
GRID_MM = 0.01       # output pixel
HALF_MM = 6.5
ROI_COLOR = "#ffd23f"
GFP_RGB = np.array([0.25, 1.0, 0.35])


def stretch(a: np.ndarray, mask: np.ndarray, lo: float, hi: float) -> np.ndarray:
    v = a[mask] if mask.any() else a.ravel()
    p1, p2 = np.percentile(v, [lo, hi])
    return np.clip((a - p1) / (p2 - p1 + 1e-9), 0, 1)


def render(mouse: int, name: str, a: orientation.Annotation, flip: bool):
    """Aligned RGB image, tissue fraction, and ROIs as [(target, centre rc mm, a, b mm, angle)]."""
    path = reconstruct.ds4_dir(mouse) / f"{name}.ome.tif"
    export_um = rois.export_pixel_um(path)
    sec = histology.load_section(path, level=LEVEL)
    bin_ = round(sec.pixel_um / export_um)
    point, resample = reconstruct.to_frame(a, export_um, bin_=bin_, grid_mm=GRID_MM, half_mm=HALF_MM)
    tissue = histology.tissue_mask(sec, factor=1, min_area_mm2=0.5)
    cy5 = np.log1p(sec.channels["CY5"].astype(np.float32))
    gfp = sec.channels["TRITC"].astype(np.float32)
    grey = stretch(cy5, tissue, 1, 99.7)
    green = stretch(gfp, tissue, 50, 99.8) ** 0.8
    planes = [resample(x) for x in (grey, green, tissue.astype(np.float32))]
    if flip:
        planes = [p[:, ::-1] for p in planes]
    grey, green, t = planes
    rgb = np.clip(grey[..., None] * 0.75 + green[..., None] * GFP_RGB[None, None, :], 0, 1)
    rgb[t < 0.02] *= 0.25   # dim glass so the tissue edge reads

    aff = orientation.canonical_affine(a)
    m = aff[:2, :2]
    mm = export_um / 1000
    out = []
    for tgt, roi in sorted(rois.read_locations().get((mouse, name), {}).items()):
        r, c = point(roi.centre)
        ea, _ = roi.axes()
        d = m @ ea                       # a-axis direction in the frame (row, col)
        if flip:
            c, d = -c, d * [1, -1]
        ang = np.degrees(np.arctan2(-d[0], d[1]))   # matplotlib: CCW from +x, y up
        out.append((tgt, (r, c), roi.semi_a * mm, roi.semi_b * mm, ang))
    return rgb, t, out


def draw(ax, rgb, rois_, coverage, extent, crop, label=True):
    ax.imshow(rgb, extent=extent, interpolation="bilinear")
    for tgt, (r, c), sa, sb, ang in rois_:
        # imshow with extent (x: col, y: row with row increasing downward)
        ax.add_patch(Ellipse((c, r), 2 * sa, 2 * sb, angle=-ang, fill=False, ec=ROI_COLOR, lw=1.0))
        if label:
            cov = coverage.get(tgt)
            txt = f"T{tgt}" + (f"\n{cov:.0%}" if cov is not None and np.isfinite(cov) else "")
            ax.text(c, r, txt, color=ROI_COLOR, ha="center", va="center", fontsize=7,
                    fontweight="bold")
    (x0, x1), (y0, y1) = crop
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def scale_and_compass(ax, crop):
    (x0, x1), (y0, y1) = crop
    ax.plot([x1 - 0.4 - 2, x1 - 0.4], [y1 - 0.35] * 2, color="white", lw=2)
    ax.text(x1 - 1.4, y1 - 0.5, "2 mm", color="white", ha="center", va="bottom", fontsize=7)
    ax.annotate("", xy=(x0 + 0.5, y0 + 0.3), xytext=(x0 + 0.5, y0 + 1.1),
                arrowprops={"arrowstyle": "-|>", "color": "white", "lw": 1})
    ax.text(x0 + 0.5, y0 + 0.2, "A", color="white", ha="center", va="bottom", fontsize=7)
    ax.text(x0 + 0.15, y0 + 1.45, "L", color="white", ha="center", va="center", fontsize=7)
    ax.text(x0 + 0.85, y0 + 1.45, "R", color="white", ha="center", va="center", fontsize=7)


def run(mouse: int) -> list[str]:
    sheet = orientation.read_sheet()
    names = sorted(n for (m, n), a in sheet.items()
                   if m == mouse and not a.exclude and None not in (a.centre, a.front, a.notch))
    flip = orientation.notch_side(mouse) == "right"
    cov = pd.read_csv(SLICES) if SLICES.exists() else pd.DataFrame()
    ax_mm = np.linspace(-HALF_MM + GRID_MM / 2, HALF_MM - GRID_MM / 2, int(round(2 * HALF_MM / GRID_MM)))
    extent = (ax_mm[0] - GRID_MM / 2, ax_mm[-1] + GRID_MM / 2, ax_mm[-1] + GRID_MM / 2, ax_mm[0] - GRID_MM / 2)

    rendered = {}
    for name in names:
        print(f"  mouse {mouse} {name}")
        rendered[name] = render(mouse, name, sheet[(mouse, name)], flip)
    union = np.max([t for _, t, _ in rendered.values()], axis=0) > 0.5
    rows = np.where(union.any(axis=1))[0]
    cols = np.where(union.any(axis=0))[0]
    pad = 0.3
    crop = ((ax_mm[cols[0]] - pad, ax_mm[cols[-1]] + pad), (ax_mm[rows[0]] - pad, ax_mm[rows[-1]] + pad))

    def coverage(name):
        if cov.empty:
            return {}
        s = cov[(cov["mouse"] == mouse) & (cov["section"] == name)]
        return dict(zip(s["target"], s["pixel_coverage"]))

    side = "flipped left-right so the animal's left is on the left" if flip else \
        "animal's left on the left (notch side)"
    out_dir = OUT / f"mouse{mouse:02d}"
    out_dir.mkdir(parents=True, exist_ok=True)
    w = crop[0][1] - crop[0][0]
    h = crop[1][1] - crop[1][0]
    for name, (rgb, t, rr) in rendered.items():
        fig, ax = plt.subplots(figsize=(6, 6 * h / w + 0.5), facecolor="black")
        draw(ax, rgb, rr, coverage(name), extent, crop)
        scale_and_compass(ax, crop)
        ax.set_title(f"Mouse {mouse} · {name} · aligned (front up, {side})", color="white",
                     fontsize=8.5, loc="left")
        fig.tight_layout()
        fig.savefig(out_dir / f"{name}.png", dpi=150, facecolor="black")
        plt.close(fig)

    # grid + outline overlay
    n = len(rendered) + 1
    ncol = 4 if n > 6 else n
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.4 * h / w * nrow + 0.8),
                             facecolor="black", squeeze=False)
    for ax, (name, (rgb, t, rr)) in zip(axes.flat, rendered.items()):
        draw(ax, rgb, rr, coverage(name), extent, crop)
        ax.set_title(name, color="white", fontsize=8, loc="left", pad=2)
    ax = axes.flat[len(rendered)]
    ax.set_facecolor("black")
    for name, (rgb, t, rr) in rendered.items():
        ax.contour(ax_mm, ax_mm, t, [0.5], colors="white", linewidths=0.5, alpha=0.6)
        for tgt, (r, c), sa, sb, ang in rr:
            ax.plot(c, r, "+", color=ROI_COLOR, ms=4, mew=0.8)
    (x0, x1), (y0, y1) = crop
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("all outlines + ROI centres", color="white", fontsize=8, loc="left", pad=2)
    scale_and_compass(ax, crop)
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.suptitle(f"Mouse {mouse}: sections aligned by the clicks — dorsal view, front up, {side}. "
                 "Grey NeuN, green anti-GFP (stretched per section); yellow = target ROIs with pixel coverage.",
                 color="white", fontsize=9, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT / f"mouse{mouse:02d}_grid.png", dpi=130, facecolor="black")
    plt.close(fig)
    return names


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--mouse", type=int, action="append", help="default: 1 and 2")
    args = p.parse_args(argv)
    for mouse in args.mouse or [1, 2]:
        names = run(mouse)
        print(f"mouse {mouse}: {len(names)} sections -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
