"""Which way each tissue section faces, recorded by a person in napari.

The target plan is nearly symmetric, so the GFP pattern alone cannot say where
the front of the brain is (the pattern turned 180 degrees fits almost as well)
or which hemisphere is which (the pattern is mirror-symmetric, and sections are
mounted face up or face down). Three clicks per section settle both:

1. **centre** -- the middle of the brain, on the midline
2. **front** -- the front (anterior) tip of the brain
3. **notch side** -- a click out to the side of the brain the orientation notch
   is on, roughly square to the front arrow (it need not be on the notch)

Centre -> front gives the plan angle for `fit_plan`. Centre -> notch side says
which side of the midline the notched hemisphere is on; only its component
across the midline is used. Clicks are stored in ``data/section_orientation.csv``
in 4x-export pixels (level 0 of ``data/processed/histology/*/ds4/*.ome.tif``).

Two napari windows use this module (``notebooks/ingest_new_histology.ipynb``
runs both):

* `annotate` -- make the clicks
* `review` -- after measuring, show GFP with the fitted target circles and the
  clicked arrows, and record whether the placement looks right

Command line, for the first one::

    .venv/bin/python -m fus.orientation data/processed/histology/Mouse_02/ds4/*.ome.tif
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from fus import histology

SHEET = histology.REPO_ROOT / "data" / "section_orientation.csv"
COLUMNS = ["mouse", "section", "centre_row", "centre_col", "front_row", "front_col",
           "notch_row", "notch_col", "exclude", "review", "notes"]

#: Pyramid level shown in napari: 5.2 um/px, enough to see a notch.
VIEW_LEVEL = 2

#: Hemisphere the lab notched: right on Mouse 1 (N. Todd, 2026-09-18), left from
#: Mouse 2 on (2026-09-29).
NOTCH_SIDE_MOUSE01 = "right"
NOTCH_SIDE_DEFAULT = "left"

#: Mouse 1 sections whose mirroring is known independently of the notch: target
#: 3 was never sonicated, and its empty site lands in slot 3 (s2, s4) or slot 6
#: (s5, s6). Their notch clicks therefore fix which side of the animal plan +x
#: (targets 1-3) is on. See docs/mouse01-dose-delivery.md, A1-A2.
KNOWN_MIRRORED = {(1, "section_s2"): False, (1, "section_s4"): False,
                  (1, "section_s5"): True, (1, "section_s6"): True}


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def anterior_angle(centre, front) -> float:
    """`PlanFit` angle (degrees) that points plan +y from ``centre`` to ``front``.

    Plan +y lands along ``(cos t, -sin t)`` in (row, col), so t = atan2(-dc, dr).
    """
    dr, dc = np.subtract(front, centre)
    return float(np.degrees(np.arctan2(-dc, dr)) % 360)


def notch_on_plus_x(angle_deg: float, centre, notch) -> bool:
    """Whether ``notch`` is on the plan +x side (slots 1-3) of the midline."""
    dr, dc = np.subtract(notch, centre)
    th = np.deg2rad(angle_deg)
    return bool(dr * np.sin(th) + dc * np.cos(th) > 0)  # plan x, up to scale


def notch_side(mouse: int) -> str:
    return NOTCH_SIDE_MOUSE01 if mouse == 1 else NOTCH_SIDE_DEFAULT


def _other(side: str) -> str:
    return "left" if side == "right" else "right"


# ---------------------------------------------------------------------------
# The sheet
# ---------------------------------------------------------------------------


@dataclass
class Annotation:
    mouse: int
    section: str
    centre: tuple[float, float] | None = None
    front: tuple[float, float] | None = None
    notch: tuple[float, float] | None = None
    exclude: bool = False
    notes: str = ""
    review: str = ""  # "", "ok" or "wrong"

    @property
    def complete(self) -> bool:
        return self.exclude or None not in (self.centre, self.front, self.notch)

    @property
    def angle_deg(self) -> float | None:
        if self.centre is None or self.front is None:
            return None
        return anterior_angle(self.centre, self.front)

    @property
    def notch_offset_deg(self) -> float | None:
        """Angle between the front and notch arrows; 90 is ideal."""
        if None in (self.centre, self.front, self.notch):
            return None
        f, n = np.subtract(self.front, self.centre), np.subtract(self.notch, self.centre)
        cos = np.dot(f, n) / (np.linalg.norm(f) * np.linalg.norm(n) + 1e-9)
        return float(np.degrees(np.arccos(np.clip(cos, -1, 1))))

    def slots_13_side(self) -> str:
        """Animal side holding slots 1-3 of this section, from its notch click."""
        side = notch_side(self.mouse)
        return side if notch_on_plus_x(self.angle_deg, self.centre, self.notch) else _other(side)

    def to_row(self) -> dict:
        def rc(p):
            return (round(p[0], 1), round(p[1], 1)) if p is not None else ("", "")
        (cr, cc), (fr, fc), (nr, nc) = rc(self.centre), rc(self.front), rc(self.notch)
        return {"mouse": self.mouse, "section": self.section, "centre_row": cr,
                "centre_col": cc, "front_row": fr, "front_col": fc, "notch_row": nr,
                "notch_col": nc, "exclude": "yes" if self.exclude else "",
                "review": self.review, "notes": self.notes}


Sheet = dict[tuple[int, str], Annotation]


def read_sheet(path: Path = SHEET) -> Sheet:
    if not Path(path).exists():
        return {}
    d = pd.read_csv(path, dtype=str).fillna("")
    for col in COLUMNS:
        if col not in d:
            d[col] = ""

    def pt(r, a, b):
        return (float(r[a]), float(r[b])) if r[a] and r[b] else None

    out = {}
    for _, r in d.iterrows():
        a = Annotation(int(r["mouse"]), r["section"], pt(r, "centre_row", "centre_col"),
                       pt(r, "front_row", "front_col"), pt(r, "notch_row", "notch_col"),
                       r["exclude"].strip().lower() == "yes", r["notes"], r["review"].strip())
        out[(a.mouse, a.section)] = a
    return out


def write_sheet(annotations: Sheet, path: Path = SHEET) -> None:
    rows = [a.to_row() for _, a in sorted(annotations.items())]
    pd.DataFrame(rows, columns=COLUMNS).to_csv(path, index=False)


def update_sheet(key: tuple[int, str], path: Path = SHEET, **fields) -> None:
    """Change only ``fields`` of one section, re-reading the file first.

    Each napari window holds its own copy of the sheet; writing that copy back
    whole let a window that stayed open erase what another had saved since (the
    review verdicts were lost that way on 2026-09-30).
    """
    current = read_sheet(path)
    a = current.setdefault(key, Annotation(*key))
    for name, value in fields.items():
        setattr(a, name, value)
    write_sheet(current, path)


def mouse_of(path) -> int:
    m = re.search(r"Mouse_0*(\d+)", str(path))
    if not m:
        raise ValueError(f"no Mouse_NN folder in {path}; pass the mouse number")
    return int(m.group(1))


def section_name(path) -> str:
    """Section name as `histology.Section.name` gives it."""
    return Path(path).name.split(".")[0]


# ---------------------------------------------------------------------------
# Slots -> targets
# ---------------------------------------------------------------------------


def plus_x_side(sheet: Sheet) -> str | None:
    """Animal side of plan +x (targets 1-3), from the Mouse 1 notch clicks.

    Returns None until at least one of the `KNOWN_MIRRORED` sections has been
    clicked. Raises if the clicked ones disagree.
    """
    sides = {}
    for key, mirrored in KNOWN_MIRRORED.items():
        a = sheet.get(key)
        if a is None or a.exclude or not a.complete:
            continue
        side = a.slots_13_side()
        sides[key[1]] = _other(side) if mirrored else side
    if not sides:
        return None
    if len(set(sides.values())) > 1:
        raise ValueError(f"Mouse 1 notch clicks disagree on the side of plan +x: {sides}")
    return next(iter(sides.values()))


def is_mirrored(a: Annotation, plus_x: str) -> bool:
    """True when slots 1-3 of this section hold targets 4-6."""
    return a.slots_13_side() != plus_x


# ---------------------------------------------------------------------------
# napari
# ---------------------------------------------------------------------------


def _panel(viewer):
    """Dock widget skeleton: layout, title and status labels, button-row factory,
    and ``note`` for wrapped text.

    The panel sits in a scroll area and every label wraps, so the dock can be
    dragged narrow without closing it.
    """
    from qtpy.QtCore import Qt
    from qtpy.QtWidgets import (QHBoxLayout, QLabel, QPushButton, QScrollArea, QSizePolicy,
                                QVBoxLayout, QWidget)

    w = QWidget()
    lay = QVBoxLayout(w)

    def note(html=""):
        label = QLabel(html)
        label.setWordWrap(True)
        label.setMinimumWidth(1)
        label.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        lay.addWidget(label)
        return label

    def buttons(*labels):
        row = QHBoxLayout()
        out = []
        for t in labels:
            b = QPushButton(t)
            b.setMinimumWidth(1)
            b.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed)
            b.setToolTip(t)
            row.addWidget(b)
            out.append(b)
        lay.addLayout(row)
        return out

    title = note()
    scroll = QScrollArea()
    scroll.setWidget(w)
    scroll.setWidgetResizable(True)
    scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    scroll.setMinimumWidth(60)
    viewer.window.add_dock_widget(scroll, name="Sections", area="right")
    return lay, title, note, buttons


def _load_view(path):
    """The section at `VIEW_LEVEL`, and the factor back to export pixels."""
    import tifffile

    sec = histology.load_section(path, level=VIEW_LEVEL)
    with tifffile.TiffFile(path) as tf:
        levels = tf.series[0].levels
        scale = levels[0].shape[-1] / levels[VIEW_LEVEL].shape[-1]
    return sec, scale


def _stretch(a: np.ndarray, lo: float = 1, hi: float = 99.5) -> np.ndarray:
    p1, p2 = np.percentile(a[::4, ::4], [lo, hi])
    return np.clip((a - p1) / (p2 - p1 + 1e-9), 0, 1)


def _add_images(viewer):
    gray = viewer.add_image(np.zeros((2, 2)), name="NeuN (CY5, log)", colormap="gray",
                            opacity=0.4)
    gfp = viewer.add_image(np.zeros((2, 2)), name="GFP stain (TRITC)", colormap="green",
                           blending="additive", opacity=0.7)
    return gray, gfp


def _show_images(gray, gfp, sec):
    gray.data = _stretch(np.log1p(sec.channels["CY5"].astype(np.float32)))
    gfp.data = _stretch(sec.channels["TRITC"].astype(np.float32), 1, 99.9)
    gray.contrast_limits = gfp.contrast_limits = (0, 1)


def _set_arrows(layer, a: Annotation, scale: float):
    vecs, colors = [], []
    for p, colour in ((a.front, "red"), (a.notch, "yellow")):
        if p is not None and a.centre is not None:
            c = np.divide(a.centre, scale)
            vecs.append([c, np.divide(p, scale) - c])
            colors.append(colour)
    layer.data = np.array(vecs).reshape(-1, 2, 2)
    if colors:
        layer.edge_color = colors


def labels_for(a: Annotation, plus_x: str | None) -> dict[int, str]:
    """``{slot: "T#"}`` when the section's mirroring is known, else ``{slot: "S#"}``."""
    mirrored = KNOWN_MIRRORED.get((a.mouse, a.section))
    if mirrored is None and plus_x is not None and a.complete and not a.exclude:
        mirrored = is_mirrored(a, plus_x)
    return {s: f"S{s}" if mirrored is None else f"T{histology.slot_to_target(s, mirrored)}"
            for s in histology.PLAN_MM}


def _draw_template(layer, a: Annotation, view_um: float, scale: float, plus_x, coverage=None):
    """Six ROI circles from the clicks, in view pixels (``scale`` = export px per view px)."""
    if a.centre is None or a.front is None:
        layer.data = np.zeros((0, 2))
        return
    fit = histology.template_fit(np.divide(a.centre, scale), a.angle_deg, view_um)
    pts = fit.points()
    names = labels_for(a, plus_x)
    layer.data = np.array([pts[s] for s in sorted(pts)])
    layer.size = 2 * histology.ROI_RADIUS_MM * fit.scale * 1000 / view_um
    text = [names[s] + (f" {coverage[s]:.0%}" if coverage and s in coverage else "")
            for s in sorted(pts)]
    layer.text = {"string": text, "color": "yellow", "size": 14}


def canonical_affine(a: Annotation, scale: float = 1.0) -> np.ndarray | None:
    """Display transform putting the front arrow straight up and the notch-side arrow left.

    A 3x3 affine on (row, col) in view pixels (export px / ``scale``), with the
    centre click at the origin. It is a rotation, plus a mirror flip when the
    section lies face down. It is a display transform only: napari draws the
    unchanged pixels through it, so nothing is resampled, and clicks and shapes
    stay in image coordinates. None until all three clicks exist.
    """
    if None in (a.centre, a.front, a.notch):
        return None
    f = np.subtract(a.front, a.centre)
    f = f / np.linalg.norm(f)
    n = np.subtract(a.notch, a.centre)
    n = n - (n @ f) * f
    if np.linalg.norm(n) < 1e-6:
        return None
    n = n / np.linalg.norm(n)
    m = -np.column_stack([f, n]).T  # f -> (-1, 0) up, n -> (0, -1) left
    out = np.eye(3)
    out[:2, :2] = m
    out[:2, 2] = -m @ np.divide(a.centre, scale)
    return out


def _straighten(viewer, a: Annotation, scale: float) -> bool:
    """Apply `canonical_affine` to every layer (identity if clicks are missing)."""
    affine = canonical_affine(a, scale)
    for layer in viewer.layers:
        layer.affine = np.eye(3) if affine is None else affine
    viewer.reset_view()
    return affine is not None


def _draw_measured(layer, rows: pd.DataFrame, scale: float, mirrored: bool | None):
    """Shapes for the ROIs in a measured table (TRITC rows of one section)."""
    from fus.rois import Roi

    layer.selected_data = set()
    layer.data = []
    if not len(rows):
        return
    shapes, labels = [], []
    for r in rows.itertuples():
        semi_a = getattr(r, "semi_a_px", np.nan)
        if pd.isna(semi_a):
            roi = Roi(0, r.slot, (r.row, r.col), r.radius_px, r.radius_px)
        else:
            roi = Roi(0, r.slot, (r.row, r.col), semi_a, r.semi_b_px, r.roi_angle_deg)
        shapes.append(roi.corners(scale))
        name = f"S{r.slot}" if mirrored is None else f"T{histology.slot_to_target(r.slot, mirrored)}"
        labels.append(f"{name} {r.coverage:.0%}")
    layer.add_ellipses(shapes)
    layer.features = pd.DataFrame({"label": labels})
    layer.text = {"string": "{label}", "color": "yellow", "size": 14}


def _add_template_layer(viewer, name="template ROIs"):
    layer = viewer.add_points(ndim=2, name=name, face_color="transparent",
                              border_color="yellow", border_width=0.05,
                              border_width_is_relative=True)
    layer.editable = False
    return layer


def _close_later(viewer):
    """Close after the button's click handler returns: closing a window from
    inside its own handler can be dropped under Jupyter's Qt event loop."""
    from qtpy.QtCore import QTimer

    def close():
        try:
            viewer.close()
        except Exception:  # noqa: BLE001 - fall back to the raw window
            viewer.window._qt_window.close()

    QTimer.singleShot(0, close)


def _compass(angle_deg: float) -> str:
    """Plan angle -> where anterior points on screen, to the nearest 45 degrees."""
    names = {90: "left", 0: "down", 270: "right", 180: "up", 45: "down-left",
             135: "up-left", 225: "up-right", 315: "down-right"}
    return names[int(round(angle_deg / 45) * 45) % 360]


def _describe(a: Annotation) -> str:
    out = ""
    if a.angle_deg is not None:
        out += f"front points {_compass(a.angle_deg)}<br>"
    off = a.notch_offset_deg
    if off is not None:
        if abs(off - 90) > 60:
            out += (f"<span style='color:orange'>notch-side arrow is {off:.0f}° from the front "
                    "arrow — too close to the midline; click further out to the side</span><br>")
        else:
            out += f"notch-side arrow {off:.0f}° from front<br>"
    return out


def _bind(viewer, layers, keys):
    for k, fn in keys.items():
        viewer.bind_key(k, fn, overwrite=True)
        for layer in layers:
            layer.bind_key(k, fn, overwrite=True)


def annotate(paths, mice=None, sheet: Path = SHEET, show: bool = True):
    """Click centre, front and notch side on each section. Returns the viewer.

    Saves to ``sheet`` whenever the section changes, and on **Done**, which
    also closes the window. Existing clicks are loaded, so it can be resumed.
    """
    import napari
    from qtpy.QtWidgets import QCheckBox, QLineEdit

    paths = [Path(p) for p in paths]
    mice = mice or [mouse_of(p) for p in paths]
    annotations = read_sheet(sheet)
    state = {"i": 0, "scale": 1.0, "loading": False}
    order = ["centre", "front", "notch"]

    viewer = napari.Viewer(title="Section orientation", show=show)
    gray, gfp = _add_images(viewer)
    arrows = viewer.add_vectors(np.zeros((0, 2, 2)), name="arrows", edge_width=6,
                                vector_style="arrow")
    arrows.editable = False
    template = _add_template_layer(viewer)
    try:
        plus_x = plus_x_side(annotations)
    except ValueError:
        plus_x = None
    layers = {
        name: viewer.add_points(ndim=2, name=label, face_color=colour, size=40,
                                border_color="white")
        for name, label, colour in (("centre", "centre", "cyan"), ("front", "front", "red"),
                                    ("notch", "notch side", "yellow"))
    }

    lay, title, note, buttons = _panel(viewer)
    b_c, b_f, b_n = buttons("Centre [C]", "Front [F]", "Notch side [N]")
    (b_straight,) = buttons("Straighten view: front up, notch left [R]")
    exclude = QCheckBox("Exclude section")
    exclude.setToolTip("No usable GFP, torn, or partial")
    lay.addWidget(exclude)
    notes = QLineEdit()
    notes.setPlaceholderText("notes (optional)")
    lay.addWidget(notes)
    b_prev, b_next = buttons("◀ Previous [B]", "Save + next ▶ [Space]")
    (b_done,) = buttons("Done — save and close")
    status = note()
    lay.addStretch()

    def current() -> Annotation:
        return annotations[(mice[state["i"]], section_name(paths[state["i"]]))]

    def pick(name):
        layer = layers[name]
        viewer.layers.selection.active = layer
        layer.mode = "add"

    def refresh():
        a = current()
        _set_arrows(arrows, a, state["scale"])
        _draw_template(template, a, state["view_um"], state["scale"], plus_x)
        done = sum(annotations.get((m, section_name(p)), Annotation(m, "")).complete
                   for p, m in zip(paths, mice))
        marks = "   ".join(f"{n}: {'✓' if getattr(a, n) else '—'}" for n in order)
        status.setText(f"{marks}   {'<b>EXCLUDED</b>' if a.exclude else ''}<br>"
                       f"{_describe(a)}<b>{done} / {len(paths)}</b> sections complete")

    def on_click(name):
        layer = layers[name]

        def handler(event=None):
            if state["loading"]:
                return
            if len(layer.data) > 1:
                state["loading"] = True
                layer.data = layer.data[-1:]
                state["loading"] = False
            p = None if len(layer.data) == 0 else tuple(float(v) * state["scale"] for v in layer.data[-1])
            setattr(current(), name, p)
            refresh()
            nxt = order.index(name) + 1
            if p is not None and nxt < len(order) and getattr(current(), order[nxt]) is None:
                pick(order[nxt])
        return handler

    for name, layer in layers.items():
        layer.events.data.connect(on_click(name))

    def load(i):
        state["i"] = i
        path, mouse = paths[i], mice[i]
        sec, state["scale"] = _load_view(path)
        state["view_um"] = sec.pixel_um
        a = annotations.setdefault((mouse, sec.name), Annotation(mouse, sec.name))
        state["loading"] = True
        _show_images(gray, gfp, sec)
        for name, layer in layers.items():
            p = getattr(a, name)
            layer.data = np.array([np.divide(p, state["scale"])]) if p else np.zeros((0, 2))
        exclude.setChecked(a.exclude)
        notes.setText(a.notes)
        state["loading"] = False
        title.setText(f"<h3>Mouse {mouse} — {sec.name}</h3>{i + 1} of {len(paths)}")
        viewer.title = f"Orientation — Mouse {mouse} — {sec.name}"
        # Straighten sections that are already clicked; leave new ones as scanned
        # so the view doesn't turn under the cursor while clicking.
        _straighten(viewer, a, state["scale"])
        refresh()
        pick(next((n for n in order if getattr(a, n) is None), "centre"))

    def save():
        a = current()
        a.exclude = exclude.isChecked()
        a.notes = notes.text()
        update_sheet((a.mouse, a.section), sheet, centre=a.centre, front=a.front,
                     notch=a.notch, exclude=a.exclude, notes=a.notes)

    def go(step):
        save()
        j = state["i"] + step
        if 0 <= j < len(paths):
            load(j)
        else:
            refresh()
            status.setText(status.text() + "<br>End of list — saved. Press <b>Done</b>.")

    def done():
        try:
            save()
        finally:
            _close_later(viewer)

    def straighten():
        if not _straighten(viewer, current(), state["scale"]):
            status.setText(status.text() + "<br>Click centre, front and notch side first.")

    b_straight.clicked.connect(straighten)
    for b, fn in ((b_c, lambda: pick("centre")), (b_f, lambda: pick("front")),
                  (b_n, lambda: pick("notch")), (b_next, lambda: go(+1)),
                  (b_prev, lambda: go(-1)), (b_done, done)):
        b.clicked.connect(fn)
    exclude.toggled.connect(lambda v: (setattr(current(), "exclude", v), refresh()))
    _bind(viewer, layers.values(), {
        "r": lambda v: straighten(),
        "c": lambda v: pick("centre"), "f": lambda v: pick("front"),
        "n": lambda v: pick("notch"), "Space": lambda v: go(+1), "b": lambda v: go(-1)})

    load(0)
    return viewer


def review(paths, slots: pd.DataFrame, mice=None, sheet: Path = SHEET, show: bool = True):
    """Fact-check the placements: GFP, fitted circles labelled by target, clicked arrows.

    ``slots`` is the table from ``fus.histology measure`` (or `analyse_section`),
    with ``row``/``col``/``radius_px`` in 4x-export pixels. **Looks right** and
    **Wrong** write ``review`` = ok / wrong to the sheet and move on; **Done**
    saves and closes.
    """
    import napari

    paths = [Path(p) for p in paths]
    mice = mice or [mouse_of(p) for p in paths]
    annotations = read_sheet(sheet)
    try:
        plus_x = plus_x_side(annotations)
    except ValueError as e:
        plus_x = None
        print(f"warning: {e}; circles are labelled by slot")
    state = {"i": 0}

    viewer = napari.Viewer(title="Placement review", show=show)
    gray, gfp = _add_images(viewer)
    circles = viewer.add_shapes(name="measured ROIs", ndim=2, edge_color="yellow",
                                face_color=[1, 1, 0, 0.05], edge_width=4)
    circles.editable = False
    arrows = viewer.add_vectors(np.zeros((0, 2, 2)), name="clicked arrows", edge_width=6,
                                vector_style="arrow", opacity=0.8)
    arrows.editable = False

    lay, title, note, buttons = _panel(viewer)
    b_ok, b_bad = buttons("Looks right ✓ [Y]", "Wrong ✗ [X]")
    b_prev, b_next = buttons("◀ Previous [B]", "Next ▶ [Space]")
    (b_done,) = buttons("Done — save and close")
    status = note()
    lay.addStretch()

    def current() -> Annotation:
        key = (mice[state["i"]], section_name(paths[state["i"]]))
        return annotations.setdefault(key, Annotation(*key))

    def load(i):
        state["i"] = i
        sec, scale = _load_view(paths[i])
        a = current()
        _show_images(gray, gfp, sec)
        rows = slots[(slots.section == sec.name) & (slots.channel == "TRITC")].sort_values("slot")
        _set_arrows(arrows, a, scale)
        mirrored = KNOWN_MIRRORED.get((a.mouse, a.section))
        if mirrored is None and plus_x is not None and a.complete and not a.exclude:
            mirrored = is_mirrored(a, plus_x)
        _draw_measured(circles, rows, scale, mirrored)

        info = []
        if not len(rows):
            info.append("<span style='color:orange'>not measured</span>")
        else:
            info.append(f"placement: {rows.get('placement', pd.Series(['fit'])).iloc[0]}")
        if mirrored is None:
            info.append("circles labelled by slot: " + (
                "section excluded or not clicked" if plus_x else "Mouse 1 notch clicks missing"))
        else:
            info.append(f"mounted {'mirrored' if mirrored else 'as scanned'}")
        verdict = {"ok": "✓ looks right", "wrong": "✗ wrong"}.get(a.review, "not reviewed")
        status.setText("<br>".join(info) + f"<br><b>{verdict}</b>")
        title.setText(f"<h3>Mouse {a.mouse} — {sec.name}</h3>{i + 1} of {len(paths)}")
        viewer.title = f"Review — Mouse {a.mouse} — {sec.name}"
        _straighten(viewer, a, scale)

    def save_verdict():
        a = current()
        update_sheet((a.mouse, a.section), sheet, review=a.review)

    def go(step):
        save_verdict()
        j = state["i"] + step
        if 0 <= j < len(paths):
            load(j)
        else:
            status.setText(status.text() + "<br>End of list — saved. Press <b>Done</b>.")

    def verdict(v):
        current().review = v
        go(+1)

    def done():
        try:
            save_verdict()
        finally:
            _close_later(viewer)

    for b, fn in ((b_ok, lambda: verdict("ok")), (b_bad, lambda: verdict("wrong")),
                  (b_prev, lambda: go(-1)), (b_next, lambda: go(+1)), (b_done, done)):
        b.clicked.connect(fn)
    _bind(viewer, [], {"y": lambda v: verdict("ok"), "x": lambda v: verdict("wrong"),
                       "Space": lambda v: go(+1), "b": lambda v: go(-1)})

    load(0)
    return viewer


def main(argv: list[str] | None = None) -> int:
    import argparse

    import napari

    p = argparse.ArgumentParser(prog="python -m fus.orientation",
                                description="Click centre, front and notch side on each section.")
    p.add_argument("sections", nargs="+", type=Path, help="exported *.ome.tif sections")
    p.add_argument("--mouse", type=int, help="mouse number (default: from the Mouse_NN folder)")
    p.add_argument("--sheet", type=Path, default=SHEET)
    args = p.parse_args(argv)
    mice = [args.mouse or mouse_of(path) for path in args.sections]
    annotate(args.sections, mice, args.sheet)
    napari.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
