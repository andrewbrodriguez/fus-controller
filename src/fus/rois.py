"""Hand-finetuned target ROIs: T1-T6 as shapes a person drags and resizes in napari.

The template (`fus.orientation`) places all six targets rigidly from three
clicks. Finetuning starts from that template and lets each target move, resize
and rotate on its own. The result is saved per target in
``data/roi_locations.csv``, in 4x-export pixels:

    centre_row, centre_col   ellipse centre
    semi_a_px, semi_b_px     half-axes
    angle_deg                direction of the ``a`` axis, degrees from +row toward +col

`histology.analyse_section(..., rois=...)` measures inside exactly these shapes.

Hand placement uses the GFP image, so a target can be pulled onto its own
signal; the review step and ``scripts/benchmark_mouse01.py`` are the checks.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from fus import histology, orientation

LOCATIONS = histology.REPO_ROOT / "data" / "roi_locations.csv"
COLUMNS = ["mouse", "section", "target", "slot", "centre_row", "centre_col",
           "semi_a_px", "semi_b_px", "angle_deg"]


@dataclass
class Roi:
    """One target's ellipse, in export pixels."""

    target: int
    slot: int
    centre: tuple[float, float]
    semi_a: float
    semi_b: float
    angle_deg: float = 0.0

    def axes(self) -> tuple[np.ndarray, np.ndarray]:
        """Unit vectors (row, col) of the a and b axes."""
        th = np.deg2rad(self.angle_deg)
        ea = np.array([np.cos(th), np.sin(th)])
        return ea, np.array([-ea[1], ea[0]])

    @property
    def radius(self) -> float:
        """Radius of the circle with the same area."""
        return float(np.sqrt(self.semi_a * self.semi_b))

    def corners(self, scale: float = 1.0) -> np.ndarray:
        """napari ellipse vertices (bounding-box corners, in order), divided by ``scale``."""
        ea, eb = self.axes()
        a, b = ea * self.semi_a, eb * self.semi_b
        c = np.asarray(self.centre)
        return np.array([c - a - b, c + a - b, c + a + b, c - a + b]) / scale

    @classmethod
    def from_corners(cls, target: int, slot: int, v: np.ndarray, scale: float = 1.0) -> Roi:
        v = np.asarray(v, float) * scale
        a, b = (v[1] - v[0]) / 2, (v[3] - v[0]) / 2
        angle = float(np.degrees(np.arctan2(a[1], a[0])))
        return cls(target, slot, tuple(v.mean(axis=0)), float(np.linalg.norm(a)),
                   float(np.linalg.norm(b)), angle)

    def mask(self, shape: tuple[int, int]) -> tuple[tuple[slice, slice], np.ndarray]:
        """Bounding-box slices into an image of ``shape``, and the ellipse inside them."""
        ext = max(self.semi_a, self.semi_b)
        r, c = self.centre
        r0, r1 = int(max(0, np.floor(r - ext))), int(min(shape[0], np.ceil(r + ext) + 1))
        c0, c1 = int(max(0, np.floor(c - ext))), int(min(shape[1], np.ceil(c + ext) + 1))
        rr, cc = np.mgrid[r0:r1, c0:c1]
        dr, dc = rr - r, cc - c
        ea, eb = self.axes()
        u = (dr * ea[0] + dc * ea[1]) / self.semi_a
        w = (dr * eb[0] + dc * eb[1]) / self.semi_b
        return (slice(r0, r1), slice(c0, c1)), u**2 + w**2 <= 1


def export_pixel_um(path) -> float:
    """Pixel size (um) of an exported section at level 0, from its OME metadata."""
    import tifffile

    with tifffile.TiffFile(path) as tf:
        return float(tf.ome_metadata.split('PhysicalSizeX="', 1)[1].split('"', 1)[0])


def from_template(a: orientation.Annotation, pixel_um: float, mirrored: bool) -> dict[int, Roi]:
    """The template's six circles as ``{target: Roi}``, in export pixels."""
    fit = histology.template_fit(a.centre, a.angle_deg, pixel_um)
    radius = histology.ROI_RADIUS_MM * fit.scale * 1000 / pixel_um
    out = {}
    for slot, (r, c) in fit.points().items():
        target = histology.slot_to_target(slot, mirrored)
        out[target] = Roi(target, slot, (r, c), radius, radius, 0.0)
    return out


# ---------------------------------------------------------------------------
# File
# ---------------------------------------------------------------------------


def read_locations(path: Path = LOCATIONS) -> dict[tuple[int, str], dict[int, Roi]]:
    """``{(mouse, section): {target: Roi}}``."""
    if not Path(path).exists():
        return {}
    d = pd.read_csv(path)
    out: dict = {}
    for r in d.itertuples():
        out.setdefault((int(r.mouse), r.section), {})[int(r.target)] = Roi(
            int(r.target), int(r.slot), (r.centre_row, r.centre_col), r.semi_a_px,
            r.semi_b_px, r.angle_deg)
    return out


def write_locations(locations: dict, path: Path = LOCATIONS) -> None:
    rows = []
    for (mouse, section), rois in sorted(locations.items()):
        for t, roi in sorted(rois.items()):
            rows.append({"mouse": mouse, "section": section, "target": t, "slot": roi.slot,
                         "centre_row": round(roi.centre[0], 1),
                         "centre_col": round(roi.centre[1], 1),
                         "semi_a_px": round(roi.semi_a, 1), "semi_b_px": round(roi.semi_b, 1),
                         "angle_deg": round(roi.angle_deg % 180, 1)})
    pd.DataFrame(rows, columns=COLUMNS).to_csv(path, index=False)


# ---------------------------------------------------------------------------
# napari
# ---------------------------------------------------------------------------


def finetune(paths, mice=None, sheet: Path = orientation.SHEET, locations: Path = LOCATIONS,
             show: bool = True):
    """Drag, resize and rotate T1-T6 on each section; **Save** writes their locations.

    Starts from saved locations when a section has them, otherwise from the
    template placed by its orientation clicks. Returns the viewer.
    """
    import napari

    paths = [Path(p) for p in paths]
    mice = mice or [orientation.mouse_of(p) for p in paths]
    clicks = orientation.read_sheet(sheet)
    saved = read_locations(locations)
    plus_x = orientation.plus_x_side(clicks)
    state = {"i": 0, "scale": 1.0, "export_um": 1.0, "rois": {}, "dirty": False}

    viewer = napari.Viewer(title="Finetune targets", show=show)
    gray, gfp = orientation._add_images(viewer)
    arrows = viewer.add_vectors(np.zeros((0, 2, 2)), name="orientation arrows", edge_width=4,
                                vector_style="arrow", opacity=0.5)
    arrows.editable = False
    shapes = viewer.add_shapes(name="targets (drag / resize)", ndim=2, edge_color="yellow",
                               face_color=[1, 1, 0, 0.08], edge_width=4)

    lay, title, note, buttons = orientation._panel(viewer)
    note("Each yellow shape is one target. With the <b>targets</b> layer selected "
         "(press <b>S</b> for select mode):<br>"
         "• click a shape, then drag it to move<br>"
         "• drag a corner handle to resize (hold <b>Shift</b> to keep it round)<br>"
         "• drag the handle above the box to rotate<br>"
         "Moving to another section, or <b>Done</b>, saves this one. Hold <b>Space</b> "
         "and drag to pan. "
         "The ROI size is the measurement: resizing changes what \"coverage\" means "
         "for that target, so keep 2 mm circles unless there's a reason not to.")
    (b_reset,) = buttons("Reset section to template")
    b_prev, b_save = buttons("◀ Previous [B]", "Save + next ▶ [Enter]")
    (b_done,) = buttons("Done — save and close")
    status = note()
    lay.addStretch()

    def key(i=None):
        i = state["i"] if i is None else i
        return (mice[i], orientation.section_name(paths[i]))

    def template_rois(k, export_um):
        a = clicks.get(k)
        if a is None or not a.complete or a.exclude:
            return None
        mirrored = orientation.KNOWN_MIRRORED.get(k)
        if mirrored is None:
            mirrored = orientation.is_mirrored(a, plus_x)
        return from_template(a, export_um, mirrored)

    def show_rois(rois):
        targets = sorted(rois)
        shapes.selected_data = set()
        shapes.data = []
        shapes.add_ellipses([rois[t].corners(state["scale"]) for t in targets])
        shapes.features = pd.DataFrame({"target": [f"T{t}" for t in targets]})
        shapes.text = {"string": "{target}", "color": "yellow", "size": 14}
        state["targets"] = targets
        shapes.mode = "select"
        viewer.layers.selection.active = shapes

    def read_back():
        """The shapes as ``{target: Roi}``; None if one was deleted or added."""
        if len(shapes.data) != len(state.get("targets", [])):
            return None
        old = state["rois"][key()]
        return {t: Roi.from_corners(t, old[t].slot, v, state["scale"])
                for t, v in zip(state["targets"], shapes.data)}

    def refresh(extra=""):
        k = key()
        done = sum(1 for i in range(len(paths)) if key(i) in saved)
        src = "saved locations" if k in saved else "template (not yet saved)"
        status.setText(f"showing: {src}<br><b>{done} / {len(paths)}</b> sections saved"
                       + (f"<br>{extra}" if extra else ""))

    def load(i):
        state["i"] = i
        k = key(i)
        sec, state["scale"] = orientation._load_view(paths[i])
        state["export_um"] = sec.pixel_um / state["scale"]
        title.setText(f"<h3>Mouse {k[0]} — {k[1]}</h3>{i + 1} of {len(paths)}")
        viewer.title = f"Finetune — Mouse {k[0]} — {k[1]}"
        orientation._show_images(gray, gfp, sec)
        a = clicks.get(k, orientation.Annotation(*k))
        orientation._set_arrows(arrows, a, state["scale"])
        rois = state["rois"].get(k) or saved.get(k) or template_rois(k, state["export_um"])
        if rois is None:
            shapes.data = []
            state["targets"] = []
            refresh("<span style='color:orange'>no orientation clicks (or excluded): "
                    "nothing to finetune</span>")
        else:
            state["rois"][k] = rois
            show_rois(rois)
            refresh()
        viewer.reset_view()

    def keep_edits() -> bool:
        if not state.get("targets"):
            return True
        rois = read_back()
        if rois is None:
            refresh("<span style='color:orange'>a shape was deleted or added; press "
                    "<b>Reset</b> to get all six back</span>")
            return False
        state["rois"][key()] = rois
        return True

    def save() -> bool:
        if not keep_edits():
            return False
        if state.get("targets"):
            saved[key()] = state["rois"][key()]
            write_locations(saved, locations)
        return True

    def go(step):
        if not save():
            return
        j = state["i"] + step
        if 0 <= j < len(paths):
            load(j)
        else:
            refresh("End of list — saved. Press <b>Done</b>.")

    def reset():
        k = key()
        rois = template_rois(k, state["export_um"])
        if rois is not None:
            state["rois"][k] = rois
            show_rois(rois)
            refresh("reset to the template (not saved yet)")

    def done():
        try:
            if not save():
                return
        except Exception:
            orientation._close_later(viewer)
            raise
        orientation._close_later(viewer)

    b_reset.clicked.connect(reset)
    b_prev.clicked.connect(lambda: go(-1))
    b_save.clicked.connect(lambda: go(+1))
    b_done.clicked.connect(done)
    orientation._bind(viewer, [shapes], {"Enter": lambda v: go(+1), "b": lambda v: go(-1)})

    load(0)
    return viewer
