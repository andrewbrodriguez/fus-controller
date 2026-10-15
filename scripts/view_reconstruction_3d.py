"""Interactive 3D viewer of the loose anti-GFP reconstruction (napari; method: ``fus.reconstruct``).

    .venv/bin/python scripts/view_reconstruction_3d.py            # Mouse 2
    .venv/bin/python scripts/view_reconstruction_3d.py --mouse 2 --order slide01_s2,...

Layers (toggle with the eye icons; opacity per layer in the controls, top left):

* **GFP (sphere depth)** -- anti-GFP+ coverage. Each positive region in a section is given
  the depth of a sphere through its middle; drawn as an isosurface, with its threshold
  slider in the layer controls.
* **tissue shell** -- see-through surface of the Gaussian-blended tissue.
* **section outlines** -- each section's tissue edge at its depth.
* **sections (raw)** -- the measured planes only, one slab per section, no inference. Off by
  default. Switch to 2D (the square button, bottom left) to step through them.
* **target columns** -- each target's mean clicked ROI centre, as a vertical line.
* **dose spheres (2f AUC)** -- one sphere per target at the stack's middle depth, scaled by
  the cumulative 2nd-harmonic AUC. By default volume ∝ dose and the highest dose has
  r = 1 mm (``reconstruct.dose_spheres``). This is the dose hypothesis drawn next to the GFP,
  not a measurement.
* **anatomy (Allen atlas)** -- one opaque atlas horizontal section in region colours, with
  labels, beneath the stack. It is fitted to the tissue outlines (scale and anterior offset;
  never to GFP), for context only: a generic brain, not this animal's anatomy.

The dock widget on the right rebuilds everything with another spacing, sphere thickness,
tissue depth blur or in-plane blur. The section order and spacing are **assumed** (see
``fus.reconstruct``). The spacing defaults to 0.25 mm, and the brain-size estimate is
printed on start. Depth increases ventrally, with the first section in the order the most
dorsal.
"""

from __future__ import annotations

import argparse

import napari
import numpy as np
from magicgui import magicgui
from napari.utils import Colormap

from fus import reconstruct as rc

ISO = 0.35
DOSE_COLOR = "#7a4fd6"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--mouse", type=int, default=2)
    p.add_argument("--order", help="comma-separated section names, first = z 0; default filename order")
    p.add_argument("--spacing-mm", type=float, default=rc.DEFAULT_SPACING_MM,
                   help="tissue-scale mm between sections (assumed)")
    p.add_argument("--no-atlas", action="store_true", help="leave out the anatomy plane")
    args = p.parse_args(argv)

    stack = rc.load_stack(args.mouse, args.order.split(",") if args.order else None)
    atlas = rc.load_atlas()
    geo = rc.geometry(stack, atlas)
    print(f"Mouse {args.mouse}: {len(stack.names)} sections in "
          f"{'filename order' if not args.order else 'the order given'}: {', '.join(stack.names)}")
    print(f"spacing: {args.spacing_mm * 1000:.0f} um (assumed). Brain size alone suggests "
          f"{geo.spacing_mm * 1000:.0f} um on the tissue scale ({geo.est.describe()} in atlas mm)")
    print(f"atlas fit: scale {geo.fit.scale:.2f}, outline overlap {geo.fit.iou:.2f}")
    plane = None if args.no_atlas else rc.atlas_plane(geo.fit, atlas)
    g = rc.GRID_MM
    y0 = rc.frame_axis()[0]

    viewer = napari.Viewer(title=f"Mouse {args.mouse} — loose 3D reconstruction (order and spacing assumed)",
                           ndisplay=3)
    layers = {}

    def build(spacing: float, thickness: float, sigma: float, smooth: float,
              dose_exponent: float = 1 / 3, dose_r_max: float = 1.0):
        z = np.arange(len(stack.names)) * spacing
        dz = spacing / 4
        gfp, zz = rc.extrude(stack.gfp, z, dz, ISO, k=thickness, h_min=spacing / 4,
                             smooth_mm=smooth)
        tissue, tz = rc.blend(stack.tissue, z, sigma * spacing, dz, smooth)
        cols = stack.target_columns()
        floor = max(zz.max(), tz.max()) + 0.15
        lines = [np.array([[zz[0], r, c], [floor, r, c]]) for r, c in cols.values()]
        tops = np.array([[zz[0], r, c] for r, c in cols.values()])
        outlines = rc.section_outlines(stack, z)
        tv, tf = rc.volume_mesh(tissue, tz, 0.5, step=2)
        spheres = rc.dose_spheres(stack, float(z.mean()), dose_r_max, dose_exponent)
        sv, sf, off = [], [], 0
        for sp in spheres:
            v, f = rc.sphere_mesh(sp["centre"], sp["radius"])
            sv.append(v)
            sf.append(f + off)
            off += len(v)
        sv, sf = np.vstack(sv), np.vstack(sf)
        dose_pts = np.array([sp["centre"] - [sp["radius"] + 0.1, 0, 0] for sp in spheres])
        dose_text = [f"T{sp['target']} {sp['dose']:.2f}" for sp in spheres]

        if plane is not None:
            pv, pf, pc = rc.plane_mesh(plane[0], floor)
            label_pts = np.array([[floor - 0.02, r, c] for _, r, c in plane[1]])
        if not layers:
            # back to front: anatomy, tissue, then the data
            if plane is not None:
                layers["anatomy"] = viewer.add_surface(
                    (pv, pf), vertex_colors=pc, name="anatomy (Allen atlas)", opacity=1.0,
                    blending="opaque", shading="none", units="mm")
                layers["anatomy labels"] = viewer.add_points(
                    label_pts, name="anatomy labels", size=0.05, face_color="#0b0b0b",
                    text={"string": [t for t, _, _ in plane[1]], "color": "#0b0b0b",
                          "size": 10, "translation": np.array([-0.05, 0, 0])}, units="mm")
            layers["tissue"] = viewer.add_surface(
                (tv, tf, np.ones(len(tv))), name="tissue shell",
                colormap=Colormap(["#b8b6ae", "#b8b6ae"], name="tissue grey"), opacity=0.10,
                blending="translucent", shading="smooth", units="mm")
            layers["outlines"] = viewer.add_shapes(
                outlines, shape_type="path", name="section outlines", edge_color="#52514e",
                edge_width=0.015, opacity=0.6, units="mm")
            layers["raw"] = viewer.add_image(
                stack.gfp, name="sections (raw)", scale=(spacing, g, g), translate=(0, y0, y0),
                colormap="magma", contrast_limits=(0, 1), rendering="mip", visible=False,
                units="mm")
            layers["gfp"] = viewer.add_image(
                gfp, name="GFP (sphere depth)", scale=(dz, g, g), translate=(zz[0], y0, y0),
                rendering="iso", iso_threshold=ISO, colormap="inferno", contrast_limits=(0, 1),
                blending="translucent", units="mm")
            layers["dose"] = viewer.add_surface(
                (sv, sf, np.ones(len(sv))), name="dose spheres (2f AUC)",
                colormap=Colormap([DOSE_COLOR, DOSE_COLOR], name="dose purple"), opacity=0.3,
                blending="translucent", shading="smooth", units="mm")
            layers["dose labels"] = viewer.add_points(
                dose_pts, name="dose labels", size=0.05, face_color=DOSE_COLOR,
                text={"string": dose_text, "color": DOSE_COLOR, "size": 10}, units="mm")
            layers["cols"] = viewer.add_shapes(
                lines, shape_type="line", name="target columns", edge_color="#2a78d6",
                edge_width=0.03, units="mm")
            layers["labels"] = viewer.add_points(
                tops, name="target labels", size=0.15, face_color="#2a78d6",
                text={"string": [f"T{t}" for t in cols], "color": "#2a78d6", "size": 12,
                      "translation": np.array([-0.15, 0, 0])}, units="mm")
        else:
            if plane is not None:
                layers["anatomy"].data = (pv, pf)
                layers["anatomy"].vertex_colors = pc
                layers["anatomy labels"].data = label_pts
            layers["tissue"].data = (tv, tf, np.ones(len(tv)))
            layers["outlines"].data = outlines
            layers["raw"].scale = (spacing, g, g)
            layers["gfp"].data = gfp
            layers["gfp"].scale, layers["gfp"].translate = (dz, g, g), (zz[0], y0, y0)
            layers["cols"].data = lines
            layers["labels"].data = tops
            layers["dose"].data = (sv, sf, np.ones(len(sv)))
            layers["dose labels"].data = dose_pts
            layers["dose labels"].text = {"string": dose_text, "color": DOSE_COLOR, "size": 10}
        viewer.status = (f"spacing {spacing * 1000:.0f} µm, sphere thickness ×{thickness:g}, "
                         f"tissue depth blur σ {sigma * spacing * 1000:.0f} µm, "
                         f"in-plane blur {smooth * 1000:.0f} µm")

    @magicgui(auto_call=False, call_button="rebuild",
              spacing_um={"min": 10, "max": 1000, "step": 5, "label": "spacing (µm)"},
              thickness={"min": 0.1, "max": 3.0, "step": 0.1, "label": "GFP thickness (× sphere)"},
              sigma={"min": 0.1, "max": 5.0, "step": 0.1, "label": "tissue depth blur (× spacing)"},
              smooth_um={"min": 0, "max": 500, "step": 10, "label": "in-plane blur (µm)"},
              dose_exponent={"min": 0.1, "max": 1.0, "step": 0.01,
                             "label": "dose sphere: radius ∝ dose^ (0.33 = volume)"},
              dose_r_max={"min": 0.1, "max": 3.0, "step": 0.1, "label": "dose sphere: max radius (mm)"})
    def controls(spacing_um: float = round(args.spacing_mm * 1000), thickness: float = 1.0,
                 sigma: float = 1.0, smooth_um: float = rc.SMOOTH_MM * 1000,
                 dose_exponent: float = 0.33, dose_r_max: float = 1.0):
        build(spacing_um / 1000, thickness, sigma, smooth_um / 1000,
              1 / 3 if abs(dose_exponent - 0.33) < 1e-9 else dose_exponent, dose_r_max)

    build(args.spacing_mm, 1.0, 1.0, rc.SMOOTH_MM)
    viewer.window.add_dock_widget(controls, area="right", name="reconstruction")
    viewer.dims.axis_labels = ("depth, ventral → (mm)", "anterior ← (mm)", "lateral, notch side ← (mm)")
    scene = getattr(viewer, "canvas", None)  # napari >= 0.9 moved these off the viewer
    (scene.overlays.scale_bar if scene else viewer.scale_bar).visible = True
    (viewer.scene.camera if hasattr(viewer, "scene") else viewer.camera).angles = (0, -35, 70)
    viewer.reset_view()
    napari.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
