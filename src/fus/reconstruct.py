"""Loose 3D reconstruction of anti-GFP delivery from one animal's horizontal sections.

N. Todd (2026-10-02): delivery is loosely a sphere around each focal spot, which is the
reason for summing delivery over slices. This module stacks the slices so the shape can
be looked at. It is a picture, not a measurement: nothing here feeds a number in ``docs/``.

1. **Section maps.** Each section's pipeline-A map (anti-GFP residual > k SD,
   ``histology.gfp_residual``) is binned to ~21 um and resampled into a common frame using
   the person's clicks: centre at the origin, front up, notch side left
   (``orientation.canonical_affine``). Rigid only: no correction for shrinkage or tilt
   between sections.
2. **Spacing** (`estimate_spacing`). Neither the depth order nor the spacing was recorded.
   Spacing is inferred from brain size: each section's tissue area, divided by the
   tissue shrinkage squared, is matched to the horizontal-section area of the Allen CCF
   2017 atlas (100 um) on the dorsal side of its widest level. The sections hold both
   hippocampi and the cerebellum, which puts them there. The spread of the matched depths
   gives the span, and span / (n - 1) gives the spacing. This assumes even spacing and
   that area differences come from depth rather than from torn or missing tissue. Tears
   shrink area and widen the span, so treat the estimate as an upper bound. For Mouse 2 it
   is ~0.1 mm (0.07-0.15 over Mouse 1's shrinkage range 0.82-0.91). **The default is
   not this estimate.** It is ``DEFAULT_SPACING_MM`` = 0.25 mm (A. Rodriguez, 2026-10-15:
   the sections are probably farther apart; 12 x 0.25 mm spans about the ~3 mm focal
   column). The estimate is printed for reference.
   **Order is not inferred.** The default is filename order. Ordering by GFP would bake
   the answer into the picture, and ordering by area is ambiguous above and below the
   widest level.
3. **Filling depth.** Tissue (`blend`): each section spreads into depth as a Gaussian of
   SD ``sigma_z``. Each plane is the weighted mean of the sections, divided by
   ``max(sum of weights, 1)``, so it fills between neighbours and fades beyond the ends.
   GFP (`extrude`): following "delivery is loosely a sphere", each GFP-positive region
   (smoothed coverage >= ``level``) is given the depth a sphere would have. A connected
   region whose deepest interior point is R from its edge gets half-thickness
   ``k * sqrt(2 R e - e^2)`` at a point e from the edge: R at the middle, 0 at the rim,
   so a round region becomes a sphere cut through its equator. Sections combine by
   **max**, so an empty neighbouring section doesn't erode a positive one. Below
   ``level`` a section gets only a thin slab (``h_min``). Everything off the section
   planes is inferred, not measured. The sphere assumes each section cut its plume
   through the middle, which overstates depth for sections that only grazed one.
4. **Anatomy** (`fit_atlas`, `atlas_plane`). One opaque Allen CCF horizontal section,
   coloured by region and labelled, drawn beneath the stack for context. It is taken at
   the sections' mean matched depth. Its in-plane placement is fitted to the sections'
   mean tissue mask, with scale and anterior offset free and the midline fixed by the
   clicks. The fit uses tissue, never GFP. Mouse 2: scale 0.90, outline overlap (IoU)
   0.79. It is a generic brain, not this animal's anatomy.
5. **Dose spheres** (`dose_spheres`). One sphere per target, centred on its column (mean
   clicked ROI centre) at the middle depth of the stack. Its radius is scaled by the
   cumulative 2nd-harmonic AUC (``results/histology/slices/doses.csv``):
   ``r = r_max * (dose / max dose) ** exponent``. The default exponent 1/3 makes **volume
   proportional to dose**, and r_max = 1.0 mm is half the ~2 mm focal spot (N. Todd,
   2026-09-18). This is the dose hypothesis drawn in the same space as the GFP, not a
   measurement.
"""

from __future__ import annotations

import gzip
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

from fus import histology, orientation, rois

REPO = histology.REPO_ROOT
ATLAS_URL = ("http://download.alleninstitute.org/informatics-archive/current-release/"
             "mouse_ccf/annotation/ccf_2017/annotation_100.nrrd")
ATLAS = REPO / "data/processed/atlas/annotation_100.nrrd"

BIN = 16             # export px per map px (1.3 um -> ~21 um)
GRID_MM = 0.04       # common-frame pixel
HALF_MM = 6.0        # common frame spans +-HALF_MM
SMOOTH_MM = 0.15     # in-plane Gaussian before blending
SHRINKAGE = histology.TEMPLATE_SCALE        # tissue mm per in-vivo mm (Mouse 1 mean)
SHRINKAGE_RANGE = (0.82, 0.91)              # Mouse 1's fitted range
DEFAULT_SPACING_MM = 0.25   # assumed (see module doc, 2); not the brain-size estimate


def ds4_dir(mouse: int) -> Path:
    return REPO / f"data/processed/histology/Mouse_{mouse:02d}/ds4"


def _cache(mouse: int) -> Path:
    return REPO / f"data/processed/histology/Mouse_{mouse:02d}/recon_maps.npz"


def section_maps(mouse: int, names: list[str]) -> dict[str, dict]:
    """Binned GFP+ fraction and tissue fraction per section, in export-px / BIN coordinates.

    Cached in ``data/processed``: each section takes ~6 s from the ds4 export.
    """
    path = _cache(mouse)
    cached = dict(np.load(path)) if path.exists() else {}
    out = {}
    for name in names:
        if f"{name}/gfp" not in cached:
            print(f"  measuring {name}")
            sec = histology.load_section(ds4_dir(mouse) / f"{name}.ome.tif")
            tissue = histology.tissue_mask(sec)
            pos = histology.gfp_residual(sec, "TRITC", tissue).positive() & tissue
            cached[f"{name}/gfp"] = histology.block_mean(pos.astype(np.float32), BIN)
            cached[f"{name}/tissue"] = histology.block_mean(tissue.astype(np.float32), BIN)
            cached[f"{name}/um"] = np.float64(sec.pixel_um)
            path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(path, **cached)
        out[name] = {k: cached[f"{name}/{k}"] for k in ("gfp", "tissue", "um")}
    return out


def to_frame(a: orientation.Annotation, um: float):
    """(export px -> frame mm) for points, and a resampler for binned maps into the frame."""
    aff = orientation.canonical_affine(a)  # export px -> centred, front up, notch left (px)
    m, t = aff[:2, :2], aff[:2, 2]
    mm_per_px = um / 1000

    def point(p):
        return (m @ np.asarray(p, float) + t) * mm_per_px

    n = int(round(2 * HALF_MM / GRID_MM))
    u0 = -HALF_MM + GRID_MM / 2
    # frame px (i, j) -> frame mm u -> export px m^T (u / mm_per_px - t) -> map px / BIN
    matrix = m.T * GRID_MM / mm_per_px / BIN
    offset = (m.T @ (np.full(2, u0) / mm_per_px - t)) / BIN

    def resample(img):
        return ndi.affine_transform(img, matrix, offset, output_shape=(n, n), order=1, cval=0.0)

    return point, resample


@dataclass
class Stack:
    """One animal's sections in the common frame, in the assumed order (first = z 0)."""

    mouse: int
    names: list[str]
    gfp: np.ndarray        # (n, rows, cols) GFP+ fraction
    tissue: np.ndarray     # (n, rows, cols) tissue fraction
    centres: dict = field(default_factory=dict)  # name -> {target: (row_mm, col_mm)}

    @property
    def areas_mm2(self) -> np.ndarray:
        return (self.tissue > 0.5).sum(axis=(1, 2)) * GRID_MM**2

    def target_columns(self) -> dict[int, np.ndarray]:
        """Mean (row_mm, col_mm) of each target's clicked ROI centre over sections."""
        pts: dict = {}
        for per in self.centres.values():
            for t, rc in per.items():
                pts.setdefault(t, []).append(rc)
        return {t: np.mean(v, axis=0) for t, v in sorted(pts.items())}


def available(mouse: int) -> list[str]:
    """Clicked sections that are neither excluded nor marked wrong, in filename order."""
    sheet = orientation.read_sheet()
    return sorted(n for (m, n), a in sheet.items()
                  if m == mouse and not a.exclude and a.review != "wrong" and a.complete)


def load_stack(mouse: int, order: list[str] | None = None) -> Stack:
    names = order or available(mouse)
    missing = set(names) - set(available(mouse))
    if missing:
        raise ValueError(f"not available (unclicked, excluded or marked wrong): {sorted(missing)}")
    sheet = orientation.read_sheet()
    maps = section_maps(mouse, names)
    locations = rois.read_locations()
    gfp, tissue, centres = [], [], {}
    for name in names:
        point, resample = to_frame(sheet[(mouse, name)], float(maps[name]["um"]))
        gfp.append(resample(maps[name]["gfp"]))
        tissue.append(resample(maps[name]["tissue"]))
        centres[name] = {t: point(r.centre) for t, r in locations.get((mouse, name), {}).items()}
    return Stack(mouse, names, np.stack(gfp), np.stack(tissue), centres)


# ---------------------------------------------------------------------------
# Spacing from brain size
# ---------------------------------------------------------------------------


def atlas_area_curve() -> tuple[np.ndarray, np.ndarray]:
    """Allen CCF 2017 horizontal-section area (mm^2) by depth (mm from the dorsal surface)."""
    if not ATLAS.exists():
        ATLAS.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(ATLAS_URL, ATLAS)
    raw = ATLAS.read_bytes()
    header, body = raw.split(b"\n\n", 1)
    sizes = [int(s) for s in header.split(b"sizes:")[1].split(b"\n")[0].split()]  # AP DV ML
    vol = np.frombuffer(gzip.decompress(body), "<u4").reshape(sizes[::-1])  # ML, DV, AP
    area = (vol > 0).sum(axis=(0, 2)) * 0.01  # 100 um voxels
    return np.arange(sizes[1]) * 0.1, area


@dataclass
class SpacingEstimate:
    spacing_mm: float
    span_mm: float
    depth_mm: np.ndarray            # matched atlas depth per section (same order as areas)
    low_mm: float                   # spacing at the high shrinkage end
    high_mm: float                  # spacing at the low shrinkage end

    def describe(self) -> str:
        return (f"{self.spacing_mm * 1000:.0f} um apart (span {self.span_mm:.2f} mm; "
                f"{self.low_mm * 1000:.0f}-{self.high_mm * 1000:.0f} um over shrinkage "
                f"{SHRINKAGE_RANGE[0]}-{SHRINKAGE_RANGE[1]})")


def estimate_spacing(areas_mm2: np.ndarray, shrinkage: float = SHRINKAGE) -> SpacingEstimate:
    """Even spacing that best explains the sections' tissue areas, from the atlas (see module doc).

    The span drops the smallest and largest match and rescales (order statistics of n evenly
    spaced points), so one torn section doesn't set it.
    """
    depth, area = atlas_area_curve()
    peak = int(area.argmax())
    up_d, up_a = depth[: peak + 1], area[: peak + 1]
    n = len(areas_mm2)

    def fit(s):
        d = np.interp(np.asarray(areas_mm2) / s**2, up_a, up_d)
        ds = np.sort(d)
        span = (ds[-2] - ds[1]) * (n - 1) / (n - 3) if n > 3 else ds[-1] - ds[0]
        return d, span

    d, span = fit(shrinkage)
    spans = [fit(s)[1] for s in SHRINKAGE_RANGE]
    return SpacingEstimate(span / (n - 1), span, d, min(spans) / (n - 1), max(spans) / (n - 1))


# ---------------------------------------------------------------------------
# Gaussian blending in depth
# ---------------------------------------------------------------------------


def blend(planes: np.ndarray, z: np.ndarray, sigma_z: float, dz: float,
          smooth_mm: float = SMOOTH_MM, pad: float = 2.0) -> tuple[np.ndarray, np.ndarray]:
    """Volume (nz, rows, cols) on a regular z grid, and that grid (mm).

    Each section is smoothed in-plane by ``smooth_mm`` and spread in depth as a Gaussian of
    SD ``sigma_z``. Plane = sum(w_i * section_i) / max(sum w_i, 1): a weighted mean where
    sections overlap, fading to zero beyond the stack (``pad`` x ``sigma_z`` past each end).
    """
    s = smooth_mm / GRID_MM
    sm = np.stack([ndi.gaussian_filter(p, s) for p in planes]).astype(np.float32)
    zz = np.arange(z.min() - pad * sigma_z, z.max() + pad * sigma_z + dz / 2, dz)
    w = np.exp(-0.5 * ((zz[:, None] - z[None, :]) / sigma_z) ** 2).astype(np.float32)
    vol = np.tensordot(w, sm, axes=(1, 0)) / np.maximum(w.sum(axis=1), 1.0)[:, None, None]
    return vol, zz


def extrude(planes: np.ndarray, z: np.ndarray, dz: float, level: float, k: float = 1.0,
            h_min: float = 0.05, smooth_mm: float = SMOOTH_MM) -> tuple[np.ndarray, np.ndarray]:
    """GFP volume with sphere-like depth inferred from each positive region (see module doc).

    Returns (nz, rows, cols) and the z grid (mm). Value at (z, x) is
    ``max_i c_i(x) * 2 ** -((z - z_i) / h_i(x)) ** 4``: flat to ~0.8 h, half at h, then gone.
    """
    s = smooth_mm / GRID_MM
    sm = np.stack([ndi.gaussian_filter(p, s) for p in planes]).astype(np.float32)
    h = np.full(sm.shape, h_min, np.float32)
    for i, c in enumerate(sm):
        labels, n = ndi.label(c >= level)
        if not n:
            continue
        e = ndi.distance_transform_edt(labels > 0) * GRID_MM
        r = ndi.maximum(e, labels, index=np.arange(1, n + 1))
        big_r = np.concatenate([[0.0], r])[labels]
        h[i] = np.maximum(k * np.sqrt(np.clip(2 * big_r * e - e**2, 0, None)), h_min)
    reach = float(h.max()) * 1.3
    zz = np.arange(z.min() - reach, z.max() + reach + dz / 2, dz)
    vol = np.zeros((len(zz), *sm.shape[1:]), np.float32)
    for zi, zk in enumerate(zz):
        near = np.abs(z - zk) <= reach
        for i in np.where(near)[0]:
            np.maximum(vol[zi], sm[i] * np.exp2(-(((zk - z[i]) / h[i]) ** 4)), out=vol[zi])
    return vol, zz


def frame_axis() -> np.ndarray:
    n = int(round(2 * HALF_MM / GRID_MM))
    return -HALF_MM + (np.arange(n) + 0.5) * GRID_MM


# ---------------------------------------------------------------------------
# Atlas overlay (anatomy for context; fitted to tissue outlines, never to GFP)
# ---------------------------------------------------------------------------

ONTOLOGY_URL = "http://api.brain-map.org/api/v2/structure_graph_download/1.json"
ONTOLOGY = REPO / "data/processed/atlas/structure_graph_1.json"
ATLAS_MM = 0.1
#: Regions labelled on the anatomy plane: (Allen acronym, label).
LABELLED = [("Isocortex", "cortex"), ("HIP", "hippocampus"), ("CP", "striatum"),
            ("TH", "thalamus"), ("CB", "cerebellum"), ("MB", "midbrain"),
            ("OLF", "olfactory"), ("fiber tracts", "white matter"), ("VS", "ventricles")]


def load_atlas() -> np.ndarray:
    """Allen CCF 2017 annotation at 100 um, indexed (AP, DV, ML): anterior, dorsal, left = 0."""
    atlas_area_curve()  # downloads if missing
    header, body = ATLAS.read_bytes().split(b"\n\n", 1)
    sizes = [int(s) for s in header.split(b"sizes:")[1].split(b"\n")[0].split()]
    return np.frombuffer(gzip.decompress(body), "<u4").reshape(sizes[::-1]).transpose(2, 1, 0)


def _ontology() -> dict:
    import json
    if not ONTOLOGY.exists():
        urllib.request.urlretrieve(ONTOLOGY_URL, ONTOLOGY)
    return json.loads(ONTOLOGY.read_text())["msg"][0]


def structure_info(acronyms) -> dict[str, dict]:
    """acronym -> {name, color '#rrggbb', ids: the structure and all its descendants}."""
    out = {}

    def ids(node):
        yield node["id"]
        for c in node["children"]:
            yield from ids(c)

    def walk(node):
        if node["acronym"] in acronyms:
            out[node["acronym"]] = {"name": node["name"], "color": "#" + node["color_hex_triplet"],
                                    "ids": np.fromiter(ids(node), np.uint32)}
        for c in node["children"]:
            walk(c)

    walk(_ontology())
    return out


@dataclass
class AtlasFit:
    """Atlas mm -> frame mm: ``frame = (atlas - centre) * scale``, per axis.

    ``scale`` and ``ap_c`` come from overlaying the atlas's horizontal outline (over the
    sections' matched depth range) on the sections' mean tissue mask. The midline ``ml_c`` is
    fixed by the centre clicks. ``dv_c`` is the mean matched depth (`estimate_spacing`),
    placed at the middle of the stack, with depth increasing ventrally. The first section in
    the assumed order is therefore the most dorsal, which is itself an assumption.
    """

    scale: float
    ap_c: float
    ml_c: float
    dv_c: float
    iou: float

    def resample(self, plane: np.ndarray, order: int = 0) -> np.ndarray:
        """An atlas horizontal plane (AP, ML) resampled into the frame grid."""
        ax = frame_axis()
        m = np.diag([GRID_MM / self.scale / ATLAS_MM] * 2)
        off = np.array([(ax[0] / self.scale + self.ap_c) / ATLAS_MM - 0.5,
                        (ax[0] / self.scale + self.ml_c) / ATLAS_MM - 0.5])
        return ndi.affine_transform(plane.astype(np.float64), m, off,
                                    output_shape=(len(ax), len(ax)), order=order, cval=0.0)


def fit_atlas(stack: Stack, est: SpacingEstimate, atlas: np.ndarray | None = None) -> AtlasFit:
    brain = (load_atlas() if atlas is None else atlas) > 0
    ml = np.where(brain.any(axis=(0, 1)))[0]
    ml_c = (ml[0] + ml[-1] + 1) / 2 * ATLAS_MM
    lo, hi = (int(round(d / ATLAS_MM)) for d in (est.depth_mm.min(), est.depth_mm.max()))
    slab = brain[:, lo: hi + 1, :].mean(axis=1)
    tissue = (stack.tissue > 0.5).mean(axis=0) > 0.5
    best = (-1.0, 0.0, 0.0)
    for s in np.arange(0.76, 1.041, 0.02):
        for ap_c in np.arange(4.0, 12.01, 0.1):
            a = AtlasFit(s, ap_c, ml_c, 0.0, 0.0).resample(slab, order=1) > 0.5
            iou = (a & tissue).sum() / max((a | tissue).sum(), 1)
            if iou > best[0]:
                best = (float(iou), float(s), float(ap_c))
    iou, s, ap_c = best
    return AtlasFit(s, ap_c, ml_c, float(est.depth_mm.mean()), iou)


def atlas_plane(fit: AtlasFit, atlas: np.ndarray | None = None, depth_mm: float | None = None):
    """Opaque anatomy image: RGBA (rows, cols, 4) in the frame grid, and region labels.

    The Allen horizontal section at ``depth_mm`` (default: the sections' mean matched depth),
    coloured with the atlas's own region colours, transparent outside the brain. Labels are
    ``[(text, row_mm, col_mm)]``, one per `LABELLED` region present, at the point deepest
    inside its largest piece.
    """
    atlas = load_atlas() if atlas is None else atlas
    dv = int(round((fit.dv_c if depth_mm is None else depth_mm) / ATLAS_MM))
    ids = fit.resample(atlas[:, dv, :], order=0).astype(np.uint32)
    colours = {}

    def walk(node):
        colours[node["id"]] = node["color_hex_triplet"]
        for c in node["children"]:
            walk(c)

    walk(_ontology())
    rgba = np.zeros((*ids.shape, 4), np.float32)
    for u in np.unique(ids):
        if u == 0:
            continue
        hexc = colours.get(int(u), "AAAAAA")
        rgba[ids == u] = [int(hexc[j:j + 2], 16) / 255 for j in (0, 2, 4)] + [1.0]
    ax = frame_axis()
    labels = []
    info = structure_info([acr for acr, _ in LABELLED])
    for acr, text in LABELLED:
        if acr not in info:
            continue
        m = np.isin(ids, info[acr]["ids"])
        if m.sum() * GRID_MM**2 < 0.3:
            continue
        lab, n = ndi.label(m)
        big = lab == (np.bincount(lab.ravel())[1:].argmax() + 1)
        r, c = np.unravel_index(ndi.distance_transform_edt(big).argmax(), big.shape)
        labels.append((text, float(ax[r]), float(ax[c])))
    return rgba, labels


def plane_mesh(rgba: np.ndarray, depth: float, stride: int = 2):
    """`atlas_plane` image as a flat mesh at ``depth``: verts (depth, row, col) mm, faces, RGBA.

    Only cells whose four corners are inside the brain get faces, so the plane has the
    brain's outline.
    """
    ax = frame_axis()[::stride]
    img = rgba[::stride, ::stride]
    nr, nc = img.shape[:2]
    rr, cc = np.meshgrid(ax[:nr], ax[:nc], indexing="ij")
    verts = np.column_stack([np.full(rr.size, depth), rr.ravel(), cc.ravel()])
    idx = np.arange(rr.size).reshape(nr, nc)
    inside = img[..., 3] > 0
    cell = inside[:-1, :-1] & inside[1:, :-1] & inside[:-1, 1:] & inside[1:, 1:]
    a, b = idx[:-1, :-1][cell], idx[1:, :-1][cell]
    c, d = idx[:-1, 1:][cell], idx[1:, 1:][cell]
    faces = np.concatenate([np.column_stack([a, b, d]), np.column_stack([a, d, c])])
    return verts, faces, img.reshape(-1, 4)


DOSES = REPO / "results/histology/slices/doses.csv"


def dose_spheres(stack: Stack, z_mid: float, r_max: float = 1.0,
                 exponent: float = 1 / 3) -> list[dict]:
    """[{target, dose, centre (depth, row, col) mm, radius mm}] for the stack's animal."""
    import pandas as pd
    d = pd.read_csv(DOSES)
    d = d[d["mouse"] == stack.mouse].set_index("target")["dose"]
    top = float(d.max())
    out = []
    for tgt, (r, c) in stack.target_columns().items():
        if tgt not in d.index or d[tgt] <= 0:
            continue
        out.append({"target": tgt, "dose": float(d[tgt]), "centre": np.array([z_mid, r, c]),
                    "radius": r_max * (float(d[tgt]) / top) ** exponent})
    return out


def sphere_mesh(centre, radius: float, n: int = 24):
    """UV sphere: verts (N, 3) around ``centre``, faces."""
    th = np.linspace(0, np.pi, n)
    ph = np.linspace(0, 2 * np.pi, 2 * n, endpoint=False)
    t, p = np.meshgrid(th, ph, indexing="ij")
    verts = np.column_stack([np.cos(t).ravel(), (np.sin(t) * np.cos(p)).ravel(),
                             (np.sin(t) * np.sin(p)).ravel()]) * radius + np.asarray(centre)
    idx = np.arange(verts.shape[0]).reshape(n, 2 * n)
    nxt = np.roll(idx, -1, axis=1)
    a, b, c, d = idx[:-1], idx[1:], nxt[:-1], nxt[1:]
    faces = np.concatenate([np.column_stack([a.ravel(), b.ravel(), d.ravel()]),
                            np.column_stack([a.ravel(), d.ravel(), c.ravel()])])
    return verts, faces


def volume_mesh(vol: np.ndarray, zz: np.ndarray, level: float, step: int = 2):
    """Mesh of a blended volume (`blend`) at ``level``, verts in frame mm (depth, row, col)."""
    from skimage import measure
    padded = np.pad(vol, ((1, 1), (0, 0), (0, 0)))
    verts, faces, _, _ = measure.marching_cubes(padded, level, step_size=step)
    ax = frame_axis()
    out = np.column_stack([np.interp(verts[:, 0] - 1, np.arange(len(zz)), zz),
                           np.interp(verts[:, 1], np.arange(len(ax)), ax),
                           np.interp(verts[:, 2], np.arange(len(ax)), ax)])
    return out, faces


def section_outlines(stack: Stack, z: np.ndarray, min_len: int = 40) -> list[np.ndarray]:
    """Each section's tissue outline as (N, 3) frame-mm paths (depth, row, col) at its z."""
    from skimage import measure
    ax = frame_axis()
    idx = np.arange(len(ax))
    paths = []
    for t, zk in zip(stack.tissue, z):
        for c in measure.find_contours(ndi.gaussian_filter(t, 1), 0.5):
            if len(c) > min_len:
                paths.append(np.column_stack([np.full(len(c), zk), np.interp(c[:, 0], idx, ax),
                                              np.interp(c[:, 1], idx, ax)]))
    return paths


@dataclass
class Geometry:
    """Spacing and atlas placement for one stack, in frame mm."""

    est: SpacingEstimate
    fit: AtlasFit

    @property
    def spacing_mm(self) -> float:
        """Inferred spacing in frame (tissue) mm: atlas-mm spacing x the fitted scale."""
        return self.est.spacing_mm * self.fit.scale


def geometry(stack: Stack, atlas: np.ndarray | None = None) -> Geometry:
    est = estimate_spacing(stack.areas_mm2)
    return Geometry(est, fit_atlas(stack, est, atlas))
