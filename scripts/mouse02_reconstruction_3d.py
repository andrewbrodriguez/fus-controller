"""Rotating GIF of Mouse 2's loose 3D anti-GFP reconstruction (method: ``fus.reconstruct``).

    .venv/bin/python scripts/mouse02_reconstruction_3d.py
    .venv/bin/python scripts/mouse02_reconstruction_3d.py --spacing-mm 0.2 \
        --order slide01_s2,slide02_s2,...   # once Nick gives the real order

Interactive version: ``scripts/view_reconstruction_3d.py``.

**Assumed, not known:** the section order (default: filename order) and the spacing
(default 0.25 mm; the brain-size estimate is printed for reference). Off the section planes,
GFP depth is inferred by treating each positive region as a sphere, and tissue is
Gaussian-blended; neither is measured (see ``fus.reconstruct``). Anatomy is one opaque Allen
atlas section drawn beneath the stack, for context.

Writes to ``results/histology/reconstruction/``:

  mouse02_3d.gif             rotating 3D (GFP, dose spheres, section outlines, target columns,
                             atlas floor) beside two fixed panels: a top-down view of the
                             measured GFP (max over sections, no inference) with dose circles,
                             and cumulative delivery vs dose share per target
  mouse02_topdown.png        the two fixed panels on their own, for slides
  mouse02_sections.png       the aligned coverage maps, to check registration by eye
  mouse02_target_spread.csv  per target, how far its clicked ROI centre moves between sections
  mouse02_spacing.csv        per section, tissue area and matched atlas depth
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy import ndimage as ndi

from fus import reconstruct as rc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/histology/reconstruction"
MOUSE = 2
ISO = 0.35           # coverage level drawn as the surface

SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GFP_COLOR = "#eb6834"
TARGET_COLOR = "#2a78d6"
DOSE_COLOR = "#7a4fd6"
AXIS = "#c3c2b7"
CUMULATIVE = ROOT / "results/histology/slices/cumulative_dose_delivery.csv"
H = rc.HALF_MM


def sections_figure(stack: rc.Stack, z, path):
    cols = 4
    rows = int(np.ceil(len(stack.names) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(3 * cols, 3.2 * rows), facecolor=SURFACE)
    ext = (-H, H, H, -H)
    ax_mm = rc.frame_axis()
    for ax, name, g, ti, zi in zip(axes.flat, stack.names, stack.gfp, stack.tissue, z):
        ax.imshow(np.where(ti > 0.5, g, np.nan), cmap="magma", vmin=0, vmax=1, extent=ext)
        ax.contour(ax_mm, ax_mm, ti, [0.5], colors=MUTED, linewidths=0.6)
        for tgt, (r, c) in stack.centres.get(name, {}).items():
            ax.plot(c, r, "o", ms=9, mfc="none", mec=TARGET_COLOR, mew=1)
            ax.text(c, r, str(tgt), color=TARGET_COLOR, ha="center", va="center", fontsize=7)
        ax.set_title(f"{name}   z = {zi:.2f} mm", fontsize=8, color=INK_2)
        ax.set_xlim(-5, 5)
        ax.set_ylim(5, -5)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    for ax in axes.flat[len(stack.names):]:
        ax.axis("off")
    fig.suptitle("Mouse 2 sections in the common frame (front up, notch side left); "
                 "colour = anti-GFP+ fraction, circles = clicked target ROIs",
                 fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)
    plt.close(fig)


def _dose_rule(exponent: float) -> str:
    return ("sphere volume ∝ dose" if abs(exponent - 1 / 3) < 1e-6
            else f"sphere radius ∝ dose^{exponent:.2g}")


def _xyz(v):
    """frame (depth, row, col) mm -> plot (lateral, anterior, dorsal-up)."""
    return np.column_stack([v[:, 2], -v[:, 1], -v[:, 0]])


def _render(fig) -> np.ndarray:
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()


def static_panels(stack: rc.Stack, z, geo: rc.Geometry, args, height_in: float) -> np.ndarray:
    """Right-hand panels, drawn once: top-down measured footprint vs dose circles, and bars."""
    fig = plt.figure(figsize=(5.4, height_in), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 2.0], hspace=0.32, top=0.93, bottom=0.08,
                          left=0.12, right=0.97)
    ext = (-H, H, H, -H)
    ax_mm = rc.frame_axis()

    # top-down: measured only (max over sections), no depth inference
    a = fig.add_subplot(gs[0])
    if not args.no_atlas:
        rgba, labels = rc.atlas_plane(geo.fit)
        bg = rgba.copy()
        bg[..., 3] *= 0.35
        a.imshow(bg, extent=ext, interpolation="nearest")
    tissue = (stack.tissue > 0.5).mean(axis=0)
    a.contour(ax_mm, ax_mm, tissue, [0.5], colors=MUTED, linewidths=0.7)
    sm = np.stack([ndi.gaussian_filter(g, rc.SMOOTH_MM / rc.GRID_MM) for g in stack.gfp]).max(axis=0)
    rgb = np.array(matplotlib.colors.to_rgb(GFP_COLOR))
    foot = np.zeros((*sm.shape, 4))
    foot[..., :3] = rgb
    foot[..., 3] = np.clip(sm / 0.8, 0, 1) * 0.9
    a.imshow(foot, extent=ext, interpolation="bilinear")
    a.contour(ax_mm, ax_mm, sm, [ISO], colors=[GFP_COLOR], linewidths=1.0)
    for sp in rc.dose_spheres(stack, 0.0, args.dose_r_max, args.dose_exponent):
        _, r, c = sp["centre"]
        a.add_patch(plt.Circle((c, r), sp["radius"], fill=False, ec=DOSE_COLOR, lw=1.4))
        a.text(c, r - sp["radius"] - 0.15, f"T{sp['target']}", color=TARGET_COLOR,
               fontsize=8, ha="center", va="bottom", fontweight="bold")
    rows = np.where((tissue > 0.5).any(axis=1))[0]
    cols = np.where((tissue > 0.5).any(axis=0))[0]
    pad = 0.4
    a.set_xlim(ax_mm[cols[0]] - pad, ax_mm[cols[-1]] + pad)
    a.set_ylim(ax_mm[rows[-1]] + pad, ax_mm[rows[0]] - pad)
    a.set_aspect("equal")
    a.set_xlabel("lateral mm (notch side −)", fontsize=8, color=INK_2)
    a.set_ylabel("posterior mm  (front ↑)", fontsize=8, color=INK_2)
    a.tick_params(labelsize=7, colors=MUTED)
    for sp_ in a.spines.values():
        sp_.set_visible(False)
    a.set_title("Top-down, measured only: GFP (max over sections) vs dose circles",
                fontsize=8.5, color=INK, loc="left")

    # bars: cumulative delivery (committed numbers) vs dose share
    b = fig.add_subplot(gs[1])
    cum = pd.read_csv(CUMULATIVE)
    cum = cum[cum["mouse"] == MOUSE].sort_values("dose")
    share = cum["dose"] / cum["dose"].mean()
    x = np.arange(len(cum))
    w = 0.26
    b.bar(x - w, share, w, color=DOSE_COLOR, alpha=0.75, label="dose ÷ animal mean")
    b.bar(x, cum["gfp_area_mm2_rel"], w, color=GFP_COLOR, label="GFP+ area (A), ÷ mean")
    b.bar(x + w, cum["gfp_cells_rel"], w, color=GFP_COLOR, alpha=0.45,
          label="GFP+ cells (B), ÷ mean")
    b.axhline(1, color=AXIS, lw=0.8)
    b.set_xticks(x, [f"T{t}\n{d:.2f}" for t, d in zip(cum["target"], cum["dose"])], fontsize=7)
    b.tick_params(axis="y", labelsize=7, colors=MUTED)
    b.tick_params(axis="x", colors=INK_2, length=0)
    for k in ("top", "right"):
        b.spines[k].set_visible(False)
    for k in ("left", "bottom"):
        b.spines[k].set_color(AXIS)
    ra = np.corrcoef(cum["dose"], cum["gfp_area_mm2_rel"])[0, 1]
    rb = np.corrcoef(cum["dose"], cum["gfp_cells_rel"])[0, 1]
    b.set_title(f"Summed over all 12 sections (targets by dose)   r: A {ra:+.2f}, B {rb:+.2f}, n = 6",
                fontsize=8.5, color=INK, loc="left")
    b.set_ylabel("× animal mean", fontsize=8, color=INK_2)
    b.legend(fontsize=6.5, frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.22))
    img = _render(fig)
    plt.close(fig)
    return img


def gif(stack: rc.Stack, z, geo: rc.Geometry, spacing_note, args, path):
    sigma_z = args.sigma * args.spacing
    dz = args.spacing / 4
    vol, zz = rc.extrude(stack.gfp, z, dz, ISO, k=args.thickness, h_min=args.spacing / 4)
    tissue, tz = rc.blend(stack.tissue, z, sigma_z, dz)

    height = 7.0
    fig = plt.figure(figsize=(7.4, height), facecolor=SURFACE)
    # fixed draw order (back to front) instead of matplotlib's per-collection depth guess,
    # which can paint the floor over the data as the view turns
    ax = fig.add_subplot(projection="3d", facecolor=SURFACE, computed_zorder=False)
    handles = []
    ax_mm = rc.frame_axis()
    lim_x, lim_y = [-5.0, 5.0], [-5.0, 5.0]
    floor = -(max(zz.max(), tz.max()) + 0.15)
    if not args.no_atlas:
        # one opaque anatomy layer: an atlas horizontal section, as a floor under the stack
        rgba, labels = rc.atlas_plane(geo.fit)
        st = 3
        img = rgba[::st, ::st]
        X, Y = np.meshgrid(ax_mm[::st], -ax_mm[::st])
        ax.plot_surface(X, Y, np.full(X.shape, floor), facecolors=img, rstride=1, cstride=1,
                        linewidth=0, antialiased=False, shade=False, zorder=0)
        for text, r, c in labels:
            ax.text(c, -r, floor + 0.02, text, fontsize=6.5, color=INK, ha="center", zorder=1)
        inside = np.where(rgba[..., 3] > 0)
        lim_x = [min(-5, ax_mm[inside[1]].min()), max(5, ax_mm[inside[1]].max())]
        lim_y = [min(-5, -ax_mm[inside[0]].max()), max(5, -ax_mm[inside[0]].min())]
        handles.append(plt.Line2D([], [], marker="s", ls="", ms=8, color="#7ED04B",
                                  label=f"Allen atlas, {geo.fit.dv_c:.1f} mm deep (region colours)"))
    # tissue: a barely-there shell plus each section's outline, so it never hides the GFP
    tv, tf = rc.volume_mesh(tissue, tz, 0.5, step=3)
    ax.add_collection3d(Poly3DCollection(_xyz(tv)[tf], facecolor=MUTED, edgecolor="none",
                                         alpha=0.04, zorder=2))
    for path_ in rc.section_outlines(stack, z):
        q = _xyz(path_)
        ax.plot(q[:, 0], q[:, 1], q[:, 2], color=INK_2, lw=0.35, alpha=0.45, zorder=3)
    handles.append(plt.Line2D([], [], color=INK_2, lw=0.8, label="section outlines (tissue)"))
    if not args.no_dose:
        for sp in rc.dose_spheres(stack, float(z.mean()), args.dose_r_max, args.dose_exponent):
            v, f = rc.sphere_mesh(sp["centre"], sp["radius"], n=16)
            ax.add_collection3d(Poly3DCollection(_xyz(v)[f], facecolor=DOSE_COLOR,
                                                 edgecolor="none", alpha=0.16, zorder=4))
        handles.append(plt.Line2D([], [], marker="o", ls="", ms=8, color=DOSE_COLOR, alpha=0.5,
                                  label=f"2nd-harmonic AUC dose: {_dose_rule(args.dose_exponent)}"
                                        f" (max r {args.dose_r_max:g} mm)"))
    gv, gf = rc.volume_mesh(vol, zz, ISO, step=2)
    ax.add_collection3d(Poly3DCollection(_xyz(gv)[gf], facecolor=GFP_COLOR, edgecolor="none",
                                         alpha=0.85, zorder=5))
    handles.append(plt.Line2D([], [], marker="s", ls="", ms=8, color=GFP_COLOR,
                              label=f"anti-GFP+ ≥ {ISO:g}, depth inferred as spheres"))
    for tgt, (r, c) in stack.target_columns().items():
        ax.plot([c, c], [-r, -r], [-zz.min(), floor], color=TARGET_COLOR, lw=1.0, zorder=6)
        ax.text(c, -r, -zz.min() + 0.25, f"T{tgt}", color=TARGET_COLOR, fontsize=8,
                ha="center", fontweight="bold", zorder=7)

    lim = [lim_x, lim_y, (floor - 0.05, -zz.min() + 0.4)]
    ax.set_xlim(*lim[0])
    ax.set_ylim(*lim[1])
    ax.set_zlim(*lim[2])
    ax.set_box_aspect([(hi - lo) * (args.z_exaggeration if k == 2 else 1)
                       for k, (lo, hi) in enumerate(lim)])
    ax.set_xlabel("lateral mm (notch side −)", fontsize=8, color=INK_2)
    ax.set_ylabel("anterior mm", fontsize=8, color=INK_2)
    ax.set_zlabel("dorsal ↑ mm" + (f" (×{args.z_exaggeration:g})" if args.z_exaggeration != 1 else ""),
                  fontsize=8, color=INK_2)
    ax.tick_params(labelsize=7, colors=MUTED)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_facecolor(SURFACE)
        a.pane.set_edgecolor(SURFACE)
    ax.legend(handles=handles, loc="lower left", fontsize=7, frameon=False,
              bbox_to_anchor=(-0.05, -0.02))
    order = "filename order, first section most dorsal" if args.order is None else "order given with --order"
    fig.suptitle("Mouse 2, loose 3D reconstruction of anti-GFP delivery", fontsize=11,
                 color=INK, y=0.97)
    fig.text(0.5, 0.895,
             f"ASSUMED: {order}; sections {spacing_note}.\n"
             f"GFP depth off the section planes inferred as spheres (×{args.thickness:g}); "
             "not measured. Atlas fitted to tissue outlines, context only.",
             ha="center", fontsize=7, color=INK_2)
    right = static_panels(stack, z, geo, args, height)
    from PIL import Image
    Image.fromarray(right).save(path.with_name(path.stem.replace("_3d", "_topdown") + ".png"))
    frames = []
    for azim in np.linspace(0, 360, args.frames, endpoint=False):
        ax.view_init(elev=args.elev, azim=azim - 60)
        left = _render(fig)
        h = min(left.shape[0], right.shape[0])
        frames.append(np.hstack([left[:h], right[:h]]))
    plt.close(fig)
    save_gif(frames, path, args.fps)
    return len(gf)


def save_gif(frames: list[np.ndarray], path: Path, fps: float) -> None:
    """One palette for every frame, so the fixed panels don't shimmer as the 3D view turns."""
    from PIL import Image
    sample = np.vstack(frames[:: max(1, len(frames) // 8)])
    palette = Image.fromarray(sample).quantize(colors=255, method=Image.Quantize.MEDIANCUT)
    images = [Image.fromarray(f).quantize(palette=palette, dither=Image.Dither.NONE) for f in frames]
    images[0].save(path, save_all=True, append_images=images[1:], duration=int(1000 / fps),
                   loop=0, optimize=False, disposal=1)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--order", help="comma-separated section names, first = z 0; default filename order")
    p.add_argument("--spacing-mm", type=float, default=rc.DEFAULT_SPACING_MM,
                   help="tissue-scale mm between sections (assumed)")
    p.add_argument("--thickness", type=float, default=1.0,
                   help="GFP depth as a multiple of the sphere inferred from each positive region")
    p.add_argument("--sigma", type=float, default=1.0, help="tissue depth blur, in section spacings")
    p.add_argument("--z-exaggeration", type=float, default=1.0)
    p.add_argument("--frames", type=int, default=48)
    p.add_argument("--fps", type=float, default=12)
    p.add_argument("--elev", type=float, default=22)
    p.add_argument("--no-dose", action="store_true", help="leave out the dose spheres")
    p.add_argument("--dose-exponent", type=float, default=1 / 3,
                   help="radius ∝ dose ** this (1/3: volume ∝ dose)")
    p.add_argument("--dose-r-max", type=float, default=1.0, help="radius of the highest dose, mm")
    p.add_argument("--no-atlas", action="store_true", help="leave out the anatomy overlay")
    p.add_argument("--suffix", default="", help="appended to output names, for alternatives")
    args = p.parse_args(argv)

    stack = rc.load_stack(MOUSE, args.order.split(",") if args.order else None)
    geo = rc.geometry(stack)
    est = geo.est
    print(f"spacing from brain size: {est.describe()} in atlas mm; "
          f"{geo.spacing_mm * 1000:.0f} um on the tissue scale (x{geo.fit.scale:.2f})")
    print(f"atlas fit: scale {geo.fit.scale:.2f}, AP centre {geo.fit.ap_c:.1f} mm, "
          f"outline overlap {geo.fit.iou:.2f}")
    args.spacing = args.spacing_mm
    spacing_note = (f"{args.spacing_mm * 1000:.0f} µm apart (brain size alone suggests "
                    f"~{geo.spacing_mm * 1000:.0f})")
    z = np.arange(len(stack.names)) * args.spacing

    spread = []
    for tgt, (r, c) in stack.target_columns().items():
        pts = np.array([cs[tgt] for cs in stack.centres.values() if tgt in cs])
        spread.append({"target": tgt, "n_sections": len(pts), "mean_lateral_mm": round(c, 2),
                       "mean_anterior_mm": round(-r, 2),
                       "sd_mm": round(float(np.sqrt(pts.var(axis=0).sum())), 2)})

    OUT.mkdir(parents=True, exist_ok=True)
    s = args.suffix
    pd.DataFrame(spread).to_csv(OUT / f"mouse02_target_spread{s}.csv", index=False)
    pd.DataFrame({"section": stack.names, "tissue_area_mm2": stack.areas_mm2.round(1),
                  "atlas_depth_mm": est.depth_mm.round(2)}).to_csv(OUT / f"mouse02_spacing{s}.csv",
                                                                   index=False)
    sections_figure(stack, z, OUT / f"mouse02_sections{s}.png")
    n_faces = gif(stack, z, geo, spacing_note, args, OUT / f"mouse02_3d{s}.gif")
    print(pd.DataFrame(spread).to_string(index=False))
    print(f"{len(stack.names)} sections, {n_faces} faces -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
