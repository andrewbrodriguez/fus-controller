"""README example: one sonicated target vs the no-FUS control, Mouse 1.

    PYTHONPATH=src python scripts/mouse01_example_figure.py

Section s4 with two ROIs -- target 1 (sonicated) and target 3 (the no-FUS
control) -- plus, per target, the acoustic emission recorded during sonication
and a close-up of the GFP it produced. Writes ``docs/figures/mouse01-example.png``.

Needs ``results/histology/mouse01_slots.csv`` (``python -m fus.histology
measure``) for ROI positions and coverage. s4 is used as scanned, so slot k is
target k. That, and target 3 being the no-FUS control, are assumptions A1/A2 in
``docs/mouse01-dose-delivery.md``. Target 1 was sonicated more than once; the
trace shown is its first run, ``Target1.mat``, which both candidate
recording-to-target mappings place at target 1.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib import patheffects
from matplotlib.patches import Circle

from fus import extract

ROOT = Path(__file__).resolve().parents[1]
SECTION = ROOT / "data/processed/histology/Mouse_01/ds4/section_s4.ome.tif"
SLOTS = ROOT / "results/histology/mouse01_slots.csv"
RECORDING = ROOT / "data/acoustic/20260611/Mouse_Cntr_01_Target1.mat"
OUT = ROOT / "docs/figures/mouse01-example.png"

SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
TRACE = "#2a78d6"
OVERVIEW_BIN = 8        # 4x export -> 32x overview (10.4 um/px)
ZOOM_HALF_MM = 1.1      # close-ups are 2.2 mm square


def stretch(a, lo, hi):
    p1, p2 = np.percentile(a, [lo, hi])
    return np.clip((a - p1) / (p2 - p1 + 1e-9), 0, 1)


def composite(neun, gfp, gfp_range):
    """NeuN in dim magenta, GFP in green, on black."""
    n = stretch(neun, 1, 99.5) ** 1.2 * 0.55
    g = np.clip((gfp - gfp_range[0]) / (gfp_range[1] - gfp_range[0]), 0, 1)
    return np.clip(np.dstack([n, g + n * 0.15, n]), 0, 1)


HALO = [patheffects.withStroke(linewidth=3, foreground="black")]


def block_mean(a, f):
    r, c = (s // f * f for s in a.shape)
    return a[:r, :c].reshape(r // f, f, c // f, f).mean(axis=(1, 3))


def load():
    with tifffile.TiffFile(SECTION) as tf:
        names = [c.split('Name="', 1)[1].split('"', 1)[0]
                 for c in tf.ome_metadata.split("<Channel ")[1:]]
        um = float(tf.ome_metadata.split('PhysicalSizeX="', 1)[1].split('"', 1)[0])
        data = dict(zip(names, tf.series[0].levels[0].asarray().astype(np.float32)))
    rois = (
        pd.read_csv(SLOTS)
        .query("section == 'section_s4' and channel == 'TRITC' and slot in [1, 3]")
        .set_index("slot")
    )
    return data, um, rois


def style_axes(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(length=0, colors=MUTED, labelsize=9)


def main() -> None:
    plt.rcParams.update({
        "font.family": ["system-ui", "-apple-system", "Helvetica Neue", "Arial", "sans-serif"],
        "text.color": INK, "axes.labelcolor": INK_2,
    })
    data, um, rois = load()
    neun, gfp = data["CY5"], data["TRITC"]
    gfp_range = np.percentile(gfp[::7, ::7], [60, 99.7])

    fig = plt.figure(figsize=(14, 6.8), facecolor=SURFACE)
    grid = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.25, 0.62], wspace=0.12, hspace=0.36,
                            left=0.01, right=0.99, top=0.85, bottom=0.11)

    # -- overview -----------------------------------------------------------
    ax = fig.add_subplot(grid[:, 0])
    f = OVERVIEW_BIN
    small_gfp = block_mean(gfp, f)
    ax.imshow(composite(block_mean(neun, f), small_gfp, np.percentile(small_gfp, [60, 99.5])))
    # target 1 labelled above its circle, target 3 to its right (dark tissue there)
    for slot, text, ls in ((1, "Target 1\nsonicated", "-"), (3, "Target 3\nno FUS (control)", "--")):
        r = rois.loc[slot]
        cx, cy, rad = (r.col + 0.5) / f - 0.5, (r.row + 0.5) / f - 0.5, r.radius_px / f
        ax.add_patch(Circle((cx, cy), rad, fill=False, ec="white", lw=2, ls=ls))
        if slot == 1:
            where = dict(xy=(cx, cy - rad), xytext=(0, 6), ha="center", va="bottom")
        else:
            where = dict(xy=(cx + rad, cy), xytext=(8, 0), ha="left", va="center")
        ax.annotate(text, textcoords="offset points", color="white", fontsize=10.5,
                    weight="bold", path_effects=HALO, **where)
    bar_px = 1000 / (um * f)  # 1 mm
    h, w = block_mean(neun, f).shape
    ax.plot([w * 0.06, w * 0.06 + bar_px], [h * 0.95] * 2, color="white", lw=3)
    ax.text(w * 0.06 + bar_px / 2, h * 0.93, "1 mm", color="white", ha="center",
            va="bottom", fontsize=10, path_effects=HALO)
    ax.text(w * 0.06, h * 0.05, "← anterior", color="white", fontsize=10, va="top",
            path_effects=HALO)
    ax.set_title("Horizontal brain section (Mouse 1, s4)", loc="left", fontsize=12, color=INK)
    ax.axis("off")

    # -- per-target rows: acoustic emission, then GFP close-up ---------------
    ex = extract(RECORDING)
    norm = ex.metrics.normalised()
    t = ex.recording.burst_times
    db = 20 * np.log10(np.maximum(norm.area_2nd, 1e-12))
    nb = ex.recording.n_baseline

    half = int(ZOOM_HALF_MM * 1000 / um)
    for i, slot in enumerate((1, 3)):
        r = rois.loc[slot]

        axa = fig.add_subplot(grid[i, 1])
        if slot == 1:
            style_axes(axa)
            axa.set_ylim(-8, 55)
            axa.set_xlim(0, t[-1] + 2)
            axa.set_ylabel("2nd harmonic (dB re baseline)", fontsize=9.5)
            axa.set_xlabel("Time (s)", fontsize=9.5)
            axa.plot(t, db, color=TRACE, lw=1.6)
            axa.annotate("microbubbles arrive", (t[nb] + 1, 20), xytext=(40, 15),
                         textcoords="data", fontsize=9, color=INK_2, va="center",
                         arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1))
            axa.set_title("Acoustic emission (first of its 3 sonications)",
                          loc="left", fontsize=11, color=INK)
        else:
            axa.set_facecolor(SURFACE)
            axa.axis("off")
            axa.text(0.5, 0.5, "Not sonicated\nno acoustic emission", transform=axa.transAxes,
                     ha="center", va="center", fontsize=12, color=MUTED, linespacing=1.6,
                     bbox=dict(boxstyle="round,pad=1.2", fc="#f0efec", ec="none"))
            axa.set_title("Acoustic emission", loc="left", fontsize=11, color=INK)

        axz = fig.add_subplot(grid[i, 2])
        r0, c0 = int(r.row) - half, int(r.col) - half
        crop = np.s_[max(r0, 0):r0 + 2 * half, max(c0, 0):c0 + 2 * half]
        img = composite(block_mean(neun[crop], 2), block_mean(gfp[crop], 2), gfp_range)
        axz.imshow(img)
        axz.add_patch(Circle(((r.col - max(c0, 0)) / 2, (r.row - max(r0, 0)) / 2), r.radius_px / 2,
                             fill=False, ec="white", lw=2, ls="-" if slot == 1 else "--"))
        axz.set_title(f"AAV delivery: {r.coverage:.0%} of ROI is GFP+", loc="left",
                      fontsize=11, color=INK)
        axz.axis("off")

    fig.text(0.01, 0.975, "Same animal, same injection: sonicated tissue takes up AAV, the control does not",
             fontsize=15, weight="bold", va="top")
    fig.text(0.01, 0.925,
             "Magenta: neurons (NeuN).  Green: GFP from the delivered AAV (anti-GFP stain).  "
             "Circles: 1.5 mm ROIs placed from the six-target plan.",
             fontsize=10, color=MUTED, va="top")
    fig.text(0.99, 0.015,
             "Target labels assume the section orientation and the slide deck's no-FUS control; "
             "see docs/mouse01-dose-delivery.md.",
             fontsize=8.5, color=MUTED, ha="right")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=120, facecolor=SURFACE)
    from PIL import Image  # shrink for the README

    Image.open(OUT).convert("RGB").quantize(200).save(OUT, optimize=True)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
