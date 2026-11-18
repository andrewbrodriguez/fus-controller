"""Rank one mouse's sections by quality, best first.

    .venv/bin/python scripts/rank_sections.py --mouse 6
    .venv/bin/python scripts/rank_sections.py --mouse 6 --size 420   # bigger thumbnails

Opens a window with one thumbnail per section, laid out in rank order (best top-left, reading
order). Every thumbnail faces the same way, resampled from the orientation clicks: front up,
T1-T3 on the right (the notch side is on the left, except for Mouse 1, whose notch was on the
animal's other side). All are on the same scale and frame, so a partial section looks partial.
NeuN is grey and the GFP stain green. The T1-T6 zones are drawn on top exactly as they are
measured: the finetuned shapes in ``data/roi_locations.csv`` where a section has them,
otherwise the template from its clicks (the caption says "template zones").

Sections without orientation clicks can't be turned to face the same way or given zones, so
they are left out and listed in the panel. Each section is contrast-stretched on its own, as
in the napari tools, so dim staining doesn't show as dim.

To rank: click a thumbnail, then **◀ Better** / **Worse ▶** (or **[** / **]**), or **Best** /
**Worst**. Dragging a thumbnail onto another also works: drop on its upper half to go before
it, lower half to go after. Type a note for the selected section (torn, folded, partial,
little GFP); it shows on the thumbnail.

The ranking is saved to ``data/section_quality.csv`` (mouse, section, rank, notes; rank 1 is
best) after every change and on **Done**. Opening the window again starts from it. The
ranking is a record, not a filter: nothing downstream reads it. To drop a section, mark it
excluded in ``data/section_orientation.csv``.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from fus import histology, orientation, rois

ROOT = histology.REPO_ROOT
QUALITY = ROOT / "data" / "section_quality.csv"
COLUMNS = ["mouse", "section", "rank", "notes"]
LEVEL = 3  # ~10 um/px: plenty for a thumbnail, loads in well under a second
MARGIN_MM = 0.3


# ---------------------------------------------------------------------------
# Ranking file
# ---------------------------------------------------------------------------


def sections_of(mouse: int) -> list[Path]:
    return sorted((ROOT / f"data/processed/histology/Mouse_{mouse:02d}/ds4").glob("*.ome.tif"))


def read_ranking(mouse: int, path: Path = QUALITY) -> tuple[list[str], dict[str, str]]:
    """Saved order (best first) and notes for one mouse."""
    if not Path(path).exists():
        return [], {}
    d = pd.read_csv(path, dtype={"section": str, "notes": str}).fillna({"notes": ""})
    d = d[d["mouse"] == mouse].sort_values("rank")
    return d["section"].tolist(), dict(zip(d["section"], d["notes"]))


def write_ranking(mouse: int, order: list[str], notes: dict[str, str],
                  path: Path = QUALITY) -> None:
    """Replace one mouse's rows, keeping every other mouse's."""
    rows = pd.DataFrame([{"mouse": mouse, "section": s, "rank": i + 1,
                          "notes": notes.get(s, "").strip()} for i, s in enumerate(order)],
                        columns=COLUMNS)
    if Path(path).exists():
        old = pd.read_csv(path, dtype={"section": str, "notes": str})
        rows = pd.concat([old[old["mouse"] != mouse], rows], ignore_index=True)
    rows.sort_values(["mouse", "rank"]).to_csv(path, index=False, columns=COLUMNS)


def initial_order(names: list[str], saved: list[str]) -> list[str]:
    """Saved order for the sections shown, then the unranked ones in filename order."""
    ranked = [s for s in saved if s in names]
    return ranked + [s for s in names if s not in ranked]


def status_of(a: orientation.Annotation | None) -> str:
    if a is None or not (a.exclude or a.complete):
        return "not clicked"
    if a.exclude:
        return "excluded"
    if a.review == "wrong":
        return "marked wrong"
    return ""


# ---------------------------------------------------------------------------
# Straightened thumbnails
# ---------------------------------------------------------------------------


@dataclass
class Zone:
    target: int
    centre: tuple[float, float]  # (row, col) in thumbnail pixels
    semi_a: float
    semi_b: float
    angle_deg: float  # a axis, degrees from +col toward +row (QPainter's rotate)


@dataclass
class Thumb:
    name: str
    status: str
    finetuned: bool
    rgb: np.ndarray  # (rows, cols, 3) uint8
    zones: list[Zone] = field(default_factory=list)


@dataclass
class _Loaded:
    name: str
    path: Path
    a: orientation.Annotation
    sec: histology.Section
    scale: float  # export px per view px
    affine: np.ndarray  # view px -> straightened px, centre click at the origin

    @property
    def mm_per_px(self) -> float:
        return self.sec.pixel_um / 1000


def _load(path: Path, a: orientation.Annotation) -> _Loaded | None:
    import tifffile

    with tifffile.TiffFile(path) as tf:
        levels = tf.series[0].levels
        level = min(LEVEL, len(levels) - 1)
        scale = levels[0].shape[-1] / levels[level].shape[-1]
    affine = orientation.canonical_affine(a, scale)
    if affine is None:
        return None
    if orientation.notch_side(a.mouse) != orientation.NOTCH_SIDE_DEFAULT:
        # Mouse 1 was notched on the animal's other side; flip it so that every mouse
        # shows the same anatomy on the same side (T1-T3 on the right).
        affine = np.diag([1.0, -1.0, 1.0]) @ affine
    sec = histology.load_section(path, level=level)
    return _Loaded(orientation.section_name(path), path, a, sec, scale, affine)


def _tissue_extent_mm(s: _Loaded) -> np.ndarray:
    """(row_lo, row_hi, col_lo, col_hi) of the tissue once straightened, in mm."""
    mask = histology.tissue_mask(s.sec, factor=1)
    rc = np.argwhere(mask[::2, ::2]) * 2.0
    if len(rc) < 100:  # no tissue found: fall back to the image corners
        h, w = mask.shape
        rc = np.array([[0, 0], [h, 0], [0, w], [h, w]], float)
    world = (rc @ s.affine[:2, :2].T + s.affine[:2, 2]) * s.mm_per_px
    lo, hi = np.percentile(world, [0.5, 99.5], axis=0)
    return np.array([lo[0], hi[0], lo[1], hi[1]])


def common_frame(loaded: list[_Loaded]) -> np.ndarray:
    """One frame (row_lo, row_hi, col_lo, col_hi, in mm) that holds every section."""
    ext = np.array([_tissue_extent_mm(s) for s in loaded])
    return np.array([ext[:, 0].min() - MARGIN_MM, ext[:, 1].max() + MARGIN_MM,
                     ext[:, 2].min() - MARGIN_MM, ext[:, 3].max() + MARGIN_MM])


def _resample(img: np.ndarray, s: _Loaded, frame: np.ndarray, px_mm: float, shape) -> np.ndarray:
    """``img`` (view px) drawn into the frame: output pixel -> mm -> straightened px -> view px."""
    m_inv = np.linalg.inv(s.affine[:2, :2])
    step = 1 / (px_mm * s.mm_per_px)  # straightened px per output px
    origin = np.array([frame[0], frame[2]]) / s.mm_per_px
    return ndi.affine_transform(img, m_inv * step, offset=m_inv @ (origin - s.affine[:2, 2]),
                                output_shape=shape, order=1, cval=0.0)


def _zones(s: _Loaded, mouse: int, saved: dict, plus_x, frame, px_mm) -> tuple[list[Zone], bool]:
    k = (mouse, s.name)
    shapes = saved.get(k)
    finetuned = shapes is not None
    if not finetuned:
        mirrored = orientation.KNOWN_MIRRORED.get(k)
        if mirrored is None:
            mirrored = orientation.is_mirrored(s.a, plus_x)
        shapes = rois.from_template(s.a, rois.export_pixel_um(s.path), mirrored)
    m, t = s.affine[:2, :2], s.affine[:2, 2]
    to_px = s.mm_per_px * px_mm  # straightened px -> thumbnail px
    out = []
    for target, r in sorted(shapes.items()):
        c = (m @ (np.asarray(r.centre) / s.scale) + t) * s.mm_per_px
        ea, _ = r.axes()
        va = m @ (ea * r.semi_a / s.scale) * to_px
        out.append(Zone(target, ((c[0] - frame[0]) * px_mm, (c[1] - frame[2]) * px_mm),
                        float(np.hypot(*va)), r.semi_b / s.scale * to_px,
                        float(np.degrees(np.arctan2(va[0], va[1])))))
    return out, finetuned


def make_thumbs(mouse: int, size: int = 320, progress=print) -> tuple[list[Thumb], list[str]]:
    """Straightened thumbnails ``size`` px tall, and the sections left out (no clicks)."""
    sheet = orientation.read_sheet()
    plus_x = orientation.plus_x_side(sheet)
    saved = rois.read_locations()
    loaded, skipped = [], []
    paths = sections_of(mouse)
    for i, p in enumerate(paths, 1):
        name = orientation.section_name(p)
        progress(f"loading {i}/{len(paths)} {name}")
        a = sheet.get((mouse, name))
        s = _load(p, a) if a is not None else None
        if s is None:
            skipped.append(name)
        else:
            loaded.append(s)
    if not loaded:
        return [], skipped

    frame = common_frame(loaded)
    px_mm = size / (frame[1] - frame[0])
    shape = (size, int(round((frame[3] - frame[2]) * px_mm)))
    thumbs = []
    for s in loaded:
        gray = _resample(orientation._stretch(np.log1p(s.sec.channels["CY5"].astype(np.float32))),
                         s, frame, px_mm, shape)
        gfp = _resample(orientation._stretch(s.sec.channels["TRITC"].astype(np.float32), 1, 99.9),
                        s, frame, px_mm, shape)
        rgb = np.repeat(0.5 * gray[..., None], 3, axis=2)
        rgb[..., 1] += 0.75 * gfp
        zones, finetuned = _zones(s, mouse, saved, plus_x, frame, px_mm)
        thumbs.append(Thumb(s.name, status_of(s.a), finetuned,
                            (np.clip(rgb, 0, 1) * 255).astype(np.uint8), zones))
    return thumbs, skipped


# ---------------------------------------------------------------------------
# Window
# ---------------------------------------------------------------------------


def paint(t: Thumb, rank: int, note: str, selected: bool = False):
    """The thumbnail as a QPixmap, with zones, a caption, the note, and a border if selected."""
    from qtpy.QtCore import QPointF, QRectF, Qt
    from qtpy.QtGui import QColor, QFont, QImage, QPainter, QPen, QPixmap

    h, w, _ = t.rgb.shape
    buf = np.ascontiguousarray(t.rgb).tobytes()
    img = QImage(buf, w, h, 3 * w, QImage.Format.Format_RGB888).copy()
    p = QPainter(img)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    font = QFont()
    font.setPixelSize(max(11, h // 26))
    font.setBold(True)
    p.setFont(font)
    yellow = QColor(255, 220, 0)
    p.setPen(QPen(yellow, 2))
    p.setBrush(Qt.BrushStyle.NoBrush)
    for z in t.zones:
        p.save()
        p.translate(z.centre[1], z.centre[0])
        p.rotate(z.angle_deg)
        p.drawEllipse(QPointF(0, 0), z.semi_a, z.semi_b)
        p.restore()
        p.drawText(QRectF(z.centre[1] - 30, z.centre[0] - 12, 60, 24),
                   Qt.AlignmentFlag.AlignCenter, f"T{z.target}")

    band = font.pixelSize() + 10
    p.fillRect(QRectF(0, 0, w, band), QColor(0, 0, 0, 170))
    p.setPen(QColor("white"))
    p.drawText(QRectF(6, 0, w - 12, band), Qt.AlignmentFlag.AlignVCenter, f"#{rank}  {t.name}")
    flags = [f for f in (t.status, "" if t.finetuned else "template zones") if f]
    if flags:
        p.setPen(QColor(255, 150, 60))
        p.drawText(QRectF(6, 0, w - 12, band),
                   Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignRight, ", ".join(flags))
    if note.strip():
        p.fillRect(QRectF(0, h - band, w, band), QColor(0, 0, 0, 170))
        p.setPen(QColor("white"))
        text = p.fontMetrics().elidedText(note.strip(), Qt.TextElideMode.ElideRight, w - 12)
        p.drawText(QRectF(6, h - band, w - 12, band), Qt.AlignmentFlag.AlignVCenter, text)
    if selected:
        p.setPen(QPen(QColor(0, 200, 255), 6))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(QRectF(3, 3, w - 6, h - 6))
    p.end()
    return QPixmap.fromImage(img)


def window(mouse: int, thumbs: list[Thumb], skipped: list[str], path: Path = QUALITY):
    """The ranking window (not shown). Returns it; ``window.ranker`` exposes the actions."""
    from qtpy.QtCore import QSize, Qt
    from qtpy.QtGui import QIcon, QKeySequence
    from qtpy.QtWidgets import (QAbstractItemView, QHBoxLayout, QLabel, QLineEdit, QListView,
                                QListWidget, QListWidgetItem, QPushButton, QShortcut,
                                QVBoxLayout, QWidget)

    by_name = {t.name: t for t in thumbs}
    saved, notes = read_ranking(mouse, path)
    order = initial_order(list(by_name), saved)
    role = Qt.ItemDataRole.UserRole

    win = QWidget()
    win.setWindowTitle(f"Rank sections — Mouse {mouse}")
    outer = QHBoxLayout(win)

    class Grid(QListWidget):
        def dropEvent(self, event):
            super().dropEvent(event)
            changed()

    grid = Grid()
    grid.setViewMode(QListView.ViewMode.ListMode)  # ListMode drops reorder the items
    grid.setFlow(QListView.Flow.LeftToRight)
    grid.setWrapping(True)
    grid.setResizeMode(QListView.ResizeMode.Adjust)
    grid.setUniformItemSizes(True)
    grid.setSpacing(6)
    h, w, _ = thumbs[0].rgb.shape
    grid.setIconSize(QSize(w, h))
    grid.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
    grid.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
    for name in order:
        item = QListWidgetItem()
        item.setData(role, name)
        grid.addItem(item)
    outer.addWidget(grid, stretch=1)

    side = QWidget()
    side.setFixedWidth(260)
    col = QVBoxLayout(side)
    head = QLabel(f"<h3>Mouse {mouse}: {len(thumbs)} sections</h3>"
                  "Best first, in reading order. Click a section, then move it.")
    head.setWordWrap(True)
    col.addWidget(head)

    def button(text, fn, row=None):
        b = QPushButton(text)
        b.clicked.connect(fn)
        (row or col).addWidget(b)
        return b

    r1, r2 = QHBoxLayout(), QHBoxLayout()
    col.addLayout(r1)
    col.addLayout(r2)
    selected_label = QLabel()
    selected_label.setWordWrap(True)
    col.addWidget(selected_label)
    col.addWidget(QLabel("Note on the selected section:"))
    edit = QLineEdit()
    edit.setPlaceholderText("e.g. torn at T2, folded, little GFP")
    col.addWidget(edit)
    if skipped:
        gone = QLabel("<span style='color:gray'>Not shown (no orientation clicks, so they "
                      "can't be turned to match): " + ", ".join(skipped) + "</span>")
        gone.setWordWrap(True)
        col.addWidget(gone)
    col.addStretch()
    status = QLabel()
    status.setWordWrap(True)
    col.addWidget(status)
    outer.addWidget(side)

    def current_order() -> list[str]:
        return [grid.item(i).data(role) for i in range(grid.count())]

    def selected() -> str | None:
        item = grid.currentItem()
        return None if item is None else item.data(role)

    def redraw():
        sel = selected()
        for i in range(grid.count()):
            item = grid.item(i)
            name = item.data(role)
            pm = paint(by_name[name], i + 1, notes.get(name, ""), name == sel)
            icon = QIcon()
            for mode in (QIcon.Mode.Normal, QIcon.Mode.Selected, QIcon.Mode.Active):
                icon.addPixmap(pm, mode)
            item.setIcon(icon)
        if sel:
            extra = status_of_name(sel)
            selected_label.setText(f"Selected: <b>#{current_order().index(sel) + 1} {sel}</b>"
                                   + (f" ({extra})" if extra else ""))

    def status_of_name(name):
        t = by_name[name]
        return ", ".join(f for f in (t.status, "" if t.finetuned else "template zones") if f)

    def save():
        write_ranking(mouse, current_order(), notes, path)
        rel = Path(path).relative_to(ROOT) if Path(path).is_relative_to(ROOT) else path
        status.setText(f"saved to {rel}")

    def changed():
        redraw()
        save()

    def move_to(row: int):
        cur = grid.currentRow()
        if cur < 0:
            return
        row = max(0, min(grid.count() - 1, row))
        if row != cur:
            item = grid.takeItem(cur)
            grid.insertItem(row, item)
            grid.setCurrentRow(row)
            changed()

    def on_select(*_):
        edit.blockSignals(True)
        edit.setText(notes.get(selected() or "", ""))
        edit.blockSignals(False)
        redraw()

    def on_edit(text):
        if selected():
            notes[selected()] = text

    def on_note_done():
        if selected():
            redraw()
            save()

    def done():
        save()
        win.close()

    button("◀ Better  [", lambda: move_to(grid.currentRow() - 1), r1)
    button("Worse ▶  ]", lambda: move_to(grid.currentRow() + 1), r1)
    button("⇤ Best", lambda: move_to(0), r2)
    button("Worst ⇥", lambda: move_to(grid.count() - 1), r2)
    for key, fn in (("[", lambda: move_to(grid.currentRow() - 1)),
                    ("]", lambda: move_to(grid.currentRow() + 1))):
        sc = QShortcut(QKeySequence(key), grid)
        sc.setContext(Qt.ShortcutContext.WidgetShortcut)
        sc.activated.connect(fn)
    b_done = QPushButton("Done — save and close")
    b_done.clicked.connect(done)
    col.insertWidget(col.count() - 1, b_done)
    grid.currentItemChanged.connect(on_select)
    edit.textEdited.connect(on_edit)
    edit.editingFinished.connect(on_note_done)

    redraw()
    grid.setCurrentRow(0)
    status.setText("starting from the saved ranking" if saved else
                   "not ranked yet: filename order. Saved once you move a section, add a "
                   "note or press Done")
    win.resize(min(1600, 4 * (w + 12) + 300), min(1000, 3 * (h + 12) + 40))
    win.ranker = {"order": current_order, "select": grid.setCurrentRow, "move_to": move_to,
                  "notes": notes, "edit": edit, "on_note_done": on_note_done, "done": done,
                  "grid": grid}
    return win


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="scripts/rank_sections.py",
                                description="Rank one mouse's sections by quality, best first.")
    p.add_argument("--mouse", type=int, required=True)
    p.add_argument("--size", type=int, default=320, help="thumbnail height in px (default 320)")
    args = p.parse_args(argv)

    from qtpy.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    thumbs, skipped = make_thumbs(args.mouse, args.size)
    for name in skipped:
        print(f"{name}: no orientation clicks, not shown")
    if not thumbs:
        print(f"Mouse {args.mouse}: no clicked sections to rank")
        return 1
    win = window(args.mouse, thumbs, skipped)
    win.show()
    app.exec()

    order, notes = read_ranking(args.mouse)
    print(f"\nMouse {args.mouse}, best first (saved to {QUALITY.relative_to(ROOT)}):")
    sheet = orientation.read_sheet()
    for i, name in enumerate(order, 1):
        extra = status_of(sheet.get((args.mouse, name)))
        print(f"  {i:2d}. {name}" + (f"  [{extra}]" if extra else "")
              + (f"  {notes[name]}" if notes.get(name) else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
