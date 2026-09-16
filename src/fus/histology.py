"""GFP area coverage at FUS target sites, measured on whole-slide sections.

The ground-truth half of the project (see ``docs/ground-truth-spec.md``): turn
an imaged brain into one delivery number per sonication target, so it can be
joined to the acoustic features from `fus.extract`.

Geometry
--------
The target plan in every recording (``data.posx.xyz``) is a 6 x 3 table in mm.
The lab's planning screenshot (``Mouse_01_AAV_Images.pptx``) draws it on an
**axial** MRI slice, and the circle spacing matches the table on both axes, so

* ``x`` is medio-lateral (targets 1-3 on one side, 4-6 mirrored on the other),
* ``y`` is antero-posterior (+y anterior),
* ``z`` is depth, and is 0 for all six targets.

The Mouse_01 sections are **horizontal** -- both hippocampi and the cerebellum
appear in every one -- so a section at the target depth crosses all six
targets at once. In these scans anterior points to the image left.

Two consequences shape the method:

1. **The plan can be placed from the GFP pattern itself.** The six targets form
   a rigid pattern, so `fit_plan` finds the rotation, translation, and
   shrinkage that best puts the pattern on the GFP plumes. To stop a target's
   own signal from pulling its ROI onto it, each target's ROI comes from a fit
   to the *other five* (leave-one-out). An empty target is then measured where
   the rest of the pattern says it should be.
2. **The image cannot say which hemisphere is which.** The plan is
   mirror-symmetric (target k <-> k +/- 3), and sections can be mounted either
   face up. Results are therefore reported per *slot* -- the plan position
   assuming targets 1-3 sit in the image's lower hemisphere when anterior faces
   left -- and `slot_to_target` maps slots to targets once the orientation of a
   section is known from outside the image.

Pipeline
--------
::

    Image.vsi
        |  export_sections()     QuPath convert-ome: one OME-TIFF per section
        v
    load_section()               4 channels at 1.3 um/px (4x downsample)
        |
        |- tissue_mask()         log(CY5) + log(TRITC), Otsu
        |- gfp_residual()        signal minus a smooth local background
        |- fit_plan()            six-target pattern on the GFP plumes
        '- measure_slots()       coverage + intensity inside each ROI

Usage
-----
>>> from fus import histology
>>> sec = histology.load_section("data/processed/histology/Mouse_01/ds4/section_s4.ome.tif")
>>> result = histology.analyse_section(sec)
>>> result.table()
"""

from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass, field
from glob import glob
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from scipy import optimize
from skimage import filters, measure, morphology

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

#: Target plan in mm, ``{target: (x, y)}``, from ``data.posx.xyz``. Identical in
#: every Mouse_Cntr_01 recording; z is 0 for all six.
PLAN_MM: dict[int, tuple[float, float]] = {
    1: (1.5, 5.0),
    2: (3.5, 3.0),
    3: (1.5, 2.5),
    4: (-1.5, 5.0),
    5: (-3.5, 3.0),
    6: (-1.5, 2.5),
}

#: ROI radius, in plan mm (scaled by the fitted shrinkage). Fixed before any
#: coverage was computed. The closest pair of targets (2 and 3) is 2.06 mm
#: apart, so 0.75 mm keeps neighbouring ROIs from overlapping. Provisional
#: until the focal spot size at 837 kHz is confirmed.
ROI_RADIUS_MM = 0.75

#: GFP+ means more than this many noise SDs above local background. The
#: spec's sensitivity analysis reports k-1 and k+1 alongside.
THRESHOLD_K = 5.0

#: Width of the local background estimate. Wider than a plume (~1 mm), so a
#: plume is not absorbed into its own background.
BACKGROUND_SIGMA_MM = 1.0

#: Channels as named in the Olympus VS scan.
GFP_CHANNELS = ("FITC", "TRITC")  # native GFP, anti-GFP stain

#: Anatomy says anterior faces left in the Mouse_01 scans, which is theta = 90
#: degrees in `PlanFit`. A fit further than this from it gets flagged.
EXPECTED_ANGLE_DEG = 90.0
ANGLE_WARN_DEG = 30.0

#: Mirror pairs: the target on the other side of the midline.
MIRROR = {1: 4, 2: 5, 3: 6, 4: 1, 5: 2, 6: 3}

REPO_ROOT = Path(__file__).resolve().parents[2]


# ---------------------------------------------------------------------------
# Export from the slide scanner format
# ---------------------------------------------------------------------------


def _qupath() -> str:
    """Path to the QuPath launcher, from ``$QUPATH`` or /Applications."""
    if os.environ.get("QUPATH"):
        return os.environ["QUPATH"]
    hits = sorted(glob("/Applications/QuPath*.app/Contents/MacOS/QuPath*"))
    if not hits:
        raise FileNotFoundError("QuPath not found; set $QUPATH to its launcher")
    return hits[-1]


def list_series(vsi: str | Path) -> list[dict]:
    """Every image series in a whole-slide file, via QuPath's Bio-Formats.

    A ``.vsi`` holds a label image, an overview, and one series per scanned
    section. Each dict has ``index``, ``name``, ``width``, ``height``,
    ``pixel_um``, ``channels``, ``downsamples``.
    """
    script = REPO_ROOT / "scripts" / "qupath" / "list_series.groovy"
    out = subprocess.run(
        [_qupath(), "script", "-a", str(Path(vsi).resolve()), str(script)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    for line in out.splitlines():
        if "SERIES_JSON " in line:
            return json.loads(line.split("SERIES_JSON ", 1)[1])
    raise RuntimeError(f"no series listing in QuPath output for {vsi}")


def export_sections(
    vsi: str | Path, out_dir: str | Path, downsample: float = 4
) -> list[Path]:
    """Write each fluorescence section of a ``.vsi`` as a pyramidal OME-TIFF.

    Sections are the series with more than 3 channels (the label and overview
    images are RGB). Output is ``section_s<series>.ome.tif``. At the default
    4x a section is ~9000 px square at 1.3 um/px, ~0.4 GB compressed.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for s in list_series(vsi):
        if len(s["channels"]) <= 3:
            continue
        dest = out_dir / f"section_s{s['index']}.ome.tif"
        subprocess.run(
            [
                _qupath(),
                "convert-ome",
                f"--series={s['index']}",
                "-d",
                str(downsample),
                "-c",
                "ZLIB",
                "--overwrite",
                str(Path(vsi).resolve()),
                str(dest.resolve()),
            ],
            capture_output=True,
            check=True,
        )
        written.append(dest)
    return written


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


@dataclass
class Section:
    """One tissue section, all channels, at the export resolution.

    ``channels`` maps the scanner's channel names (DAPI, FITC, TRITC, CY5) to
    ``(rows, cols)`` uint16 arrays.
    """

    path: Path
    channels: dict[str, np.ndarray]
    pixel_um: float

    @property
    def name(self) -> str:
        return self.path.name.split(".")[0]

    @property
    def shape(self) -> tuple[int, int]:
        return next(iter(self.channels.values())).shape

    def binned(self, channel: str, factor: int) -> np.ndarray:
        """Block-mean of a channel, as float32, for coarse-scale work."""
        return block_mean(self.channels[channel].astype(np.float32), factor)


def load_section(path: str | Path, level: int = 0) -> Section:
    """Read one pyramid level of an exported OME-TIFF section.

    Level 0 is the export resolution (1.3 um/px at 4x); each further level
    halves it. Level 3 (~10 um/px) loads in well under a second and is enough
    to see the whole section.
    """
    path = Path(path)
    with tifffile.TiffFile(path) as tf:
        pixels = tf.ome_metadata.split("<Pixels ", 1)[1]
        pixel_um = float(pixels.split('PhysicalSizeX="', 1)[1].split('"', 1)[0])
        names = [
            chunk.split('Name="', 1)[1].split('"', 1)[0]
            for chunk in tf.ome_metadata.split("<Channel ")[1:]
        ]
        levels = tf.series[0].levels
        data = levels[level].asarray()
        pixel_um *= levels[0].shape[-1] / data.shape[-1]
    if data.ndim != 3 or data.shape[0] != len(names):
        raise ValueError(f"{path}: expected (C, Y, X), got {data.shape}")
    return Section(path, dict(zip(names, data)), pixel_um)


def block_mean(img: np.ndarray, factor: int) -> np.ndarray:
    """Downsample by averaging ``factor`` x ``factor`` blocks (edges cropped)."""
    r, c = (s // factor * factor for s in img.shape)
    return img[:r, :c].reshape(r // factor, factor, c // factor, factor).mean(axis=(1, 3))


def upsample(img: np.ndarray, factor: int, shape: tuple[int, int]) -> np.ndarray:
    """Inverse of `block_mean`: bilinear, then padded or cropped to ``shape``."""
    up = ndi.zoom(img, factor, order=1, grid_mode=True, mode="nearest")
    out = np.empty(shape, dtype=up.dtype)
    r, c = min(shape[0], up.shape[0]), min(shape[1], up.shape[1])
    out[:r, :c] = up[:r, :c]
    out[r:, :] = out[r - 1 : r, :]
    out[:, c:] = out[:, c - 1 : c]
    return out


# ---------------------------------------------------------------------------
# Tissue and signal
# ---------------------------------------------------------------------------


def tissue_mask(sec: Section, factor: int = 4, min_area_mm2: float = 0.5) -> np.ndarray:
    """Boolean mask of tissue at full export resolution.

    Glass reads ~10-40 counts in TRITC and CY5 and tissue reads hundreds, so on
    a log scale Otsu separates them cleanly -- including the dim deep nuclei,
    which a linear threshold on NeuN loses. DAPI is not used: out-of-focus
    haze makes the glass beside some sections as bright as tissue.

    Holes are filled only if smaller than 0.05 mm^2, so ventricles and tears
    stay excluded. Fragments smaller than ``min_area_mm2`` are dropped.
    """
    um = sec.pixel_um * factor
    signal = sum(np.log1p(sec.binned(ch, factor)) for ch in ("CY5", "TRITC"))
    signal = ndi.gaussian_filter(signal, 1.5)
    mask = signal > filters.threshold_otsu(signal)
    # scipy rather than skimage.morphology: skimage renamed these parameters in
    # 0.26 and drops binary_opening in 0.28. The border values reproduce
    # skimage's opening exactly (erode with border on, dilate with it off).
    disk = morphology.disk(1)
    mask = ndi.binary_dilation(ndi.binary_erosion(mask, disk, border_value=1), disk, border_value=0)
    mask = ~_drop_small(~mask, int(0.05e6 / um**2))
    mask = _drop_small(mask, int(min_area_mm2 * 1e6 / um**2))
    return upsample(mask.astype(np.float32), factor, sec.shape) > 0.5


def _drop_small(mask: np.ndarray, min_px: int) -> np.ndarray:
    """Remove 4-connected components smaller than ``min_px`` pixels."""
    labels, _ = ndi.label(mask)
    sizes = np.bincount(labels.ravel())
    keep = sizes >= min_px
    keep[0] = False
    return keep[labels]


@dataclass
class Residual:
    """One GFP channel with its local background removed.

    ``values`` is signal minus background, full resolution. ``noise_sd`` is the
    robust SD of ``values`` over tissue that is not signal; GFP+ is
    ``values > k * noise_sd``.
    """

    channel: str
    values: np.ndarray
    background: np.ndarray  # coarse, for figures
    noise_sd: float

    def positive(self, k: float = THRESHOLD_K) -> np.ndarray:
        return self.values > k * self.noise_sd


def gfp_residual(
    sec: Section,
    channel: str,
    tissue: np.ndarray,
    factor: int = 8,
    sigma_mm: float = BACKGROUND_SIGMA_MM,
    k: float = THRESHOLD_K,
    iterations: int = 3,
) -> Residual:
    """Subtract a smooth, plume-excluding background from a GFP channel.

    Autofluorescence differs between regions of the same section (s5 has two
    distinct background levels), so a single threshold per section marks whole
    bright regions as GFP+. Instead the background is a normalised Gaussian
    average of tissue pixels at ``factor`` x coarser scale, recomputed
    ``iterations`` times with pixels already above threshold left out, so
    plumes do not raise their own background.

    Tile-seam shading (~0.67 mm period) is narrower than ``sigma_mm`` and is
    *not* removed here.
    """
    img = sec.binned(channel, factor)
    t = block_mean(tissue.astype(np.float32), factor) > 0.5
    sigma = sigma_mm * 1000 / (sec.pixel_um * factor)
    weight = t.astype(np.float32)
    for _ in range(iterations):
        num = ndi.gaussian_filter(img * weight, sigma)
        den = ndi.gaussian_filter(weight, sigma)
        bg = num / np.maximum(den, 1e-3)
        resid = img - bg
        v = resid[weight > 0]
        sd = 1.4826 * np.median(np.abs(v - np.median(v)))
        weight = (t & (resid < k * sd)).astype(np.float32)

    # Pixel noise at full resolution is much larger than in the binned image,
    # so re-estimate it there: robust SD over tissue, then again with the
    # pixels above threshold removed.
    full = sec.channels[channel].astype(np.float32)
    full -= upsample(bg.astype(np.float32), factor, sec.shape)
    sample = full[tissue][::7]
    noise_sd = _mad_sd(sample)
    noise_sd = _mad_sd(sample[sample < k * noise_sd])
    return Residual(channel, full, bg, noise_sd)


def _mad_sd(v: np.ndarray) -> float:
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


# ---------------------------------------------------------------------------
# Placing the target plan
# ---------------------------------------------------------------------------


@dataclass
class PlanFit:
    """Similarity transform from plan mm to image pixels.

    A plan point ``(x, y)`` lands at::

        col = col0 + s * ppm * ( x cos(theta) - y sin(theta) )
        row = row0 + s * ppm * ( x sin(theta) + y cos(theta) )

    relative to the plan centroid, with ``ppm`` pixels per mm. At theta = 90
    degrees anterior (+y) faces left and +x (targets 1-3) faces down.
    ``scale`` < 1 is tissue shrinkage.
    """

    row0: float
    col0: float
    angle_deg: float
    scale: float
    pixel_um: float
    plan: dict[int, tuple[float, float]] = field(default_factory=lambda: dict(PLAN_MM))
    score: float = float("nan")

    def points(self) -> dict[int, tuple[float, float]]:
        """``{slot: (row, col)}`` in pixels."""
        return _plan_points(
            self.plan, self.row0, self.col0, self.angle_deg, self.scale, self.pixel_um
        )

    def rescaled(self, factor: float) -> PlanFit:
        """The same fit expressed on an image ``factor`` x finer."""
        return PlanFit(
            (self.row0 + 0.5) * factor - 0.5,
            (self.col0 + 0.5) * factor - 0.5,
            self.angle_deg,
            self.scale,
            self.pixel_um / factor,
            self.plan,
            self.score,
        )


def _centred(plan: dict[int, tuple[float, float]]) -> dict[int, np.ndarray]:
    xy = np.array(list(plan.values()), dtype=float)
    return {k: np.asarray(v, float) - xy.mean(axis=0) for k, v in plan.items()}


def _mm_to_px(x, y, angle_deg, scale, pixel_um):
    """Rotate and scale a plan vector (mm) into an image offset (rows, cols)."""
    th = np.deg2rad(angle_deg)
    ppm = scale * 1000.0 / pixel_um
    return ppm * (x * np.sin(th) + y * np.cos(th)), ppm * (x * np.cos(th) - y * np.sin(th))


def _plan_points(plan, row0, col0, angle_deg, scale, pixel_um):
    out = {}
    for k, (x, y) in _centred(plan).items():
        dr, dc = _mm_to_px(x, y, angle_deg, scale, pixel_um)
        out[k] = (row0 + dr, col0 + dc)
    return out


def _grid_search(density, plan, pixel_um, angles, scales):
    """Best integer placement for each (angle, scale), by summing shifted maps."""
    best = (-np.inf, None)
    rows, cols = density.shape
    for a in angles:
        for s in scales:
            offs = _plan_points(plan, 0.0, 0.0, a, s, pixel_um)
            offs = {k: (int(round(r)), int(round(c))) for k, (r, c) in offs.items()}
            pad = max(max(abs(r), abs(c)) for r, c in offs.values()) + 1
            padded = np.pad(density, pad)
            total = np.zeros_like(density)
            for dr, dc in offs.values():
                total += padded[pad + dr : pad + dr + rows, pad + dc : pad + dc + cols]
            idx = np.unravel_index(np.argmax(total), total.shape)
            if total[idx] > best[0]:
                best = (total[idx], (float(idx[0]), float(idx[1]), a, s))
    return best


def _refine(density, plan, pixel_um, start, scale_bounds=(0.7, 1.3)):
    """Continuous (row0, col0, angle, scale) maximising density at the points."""

    def neg(p):
        r0, c0, a, s = p
        if not scale_bounds[0] <= s <= scale_bounds[1]:
            return 1e9
        pts = np.array(list(_plan_points(plan, r0, c0, a, s, pixel_um).values())).T
        return -ndi.map_coordinates(density, pts, order=1, mode="constant").sum()

    res = optimize.minimize(
        neg,
        np.asarray(start, float),
        method="Nelder-Mead",
        options={"xatol": 0.05, "fatol": 1e-6, "maxiter": 4000,
                 "initial_simplex": _simplex(start)},
    )
    return res.x, -res.fun


def _simplex(start):
    r0, c0, a, s = start
    base = np.array([r0, c0, a, s], float)
    steps = np.diag([6.0, 6.0, 4.0, 0.05])
    return np.vstack([base, base + steps])


def plan_density(
    positives: list[np.ndarray], tissue: np.ndarray, factor: int, sigma_px: float
) -> np.ndarray:
    """Smoothed local GFP+ fraction at ``factor`` x coarser scale.

    The union of the channels is used, so the placement does not depend on
    which channel is later measured.
    """
    union = np.logical_or.reduce(positives) & tissue
    return ndi.gaussian_filter(block_mean(union.astype(np.float32), factor), sigma_px)


def fit_plan(
    density: np.ndarray,
    pixel_um: float,
    plan: dict[int, tuple[float, float]] | None = None,
    leave_out: int | None = None,
    start: PlanFit | None = None,
) -> PlanFit:
    """Place the six-target plan on a GFP density map.

    Without ``start``, searches all rotations (5 degree steps) and scales
    0.7-1.3 at integer offsets, then refines. With ``start``, only refines
    from it -- used for the leave-one-out fits, which start from the full fit
    so that dropping a target cannot make the pattern jump to a different
    arrangement of plumes.

    ``leave_out`` drops one target from the objective; its point is still
    placed by the transform.
    """
    plan = dict(plan or PLAN_MM)
    used = {k: v for k, v in plan.items() if k != leave_out}
    if start is None:
        # Coarse search on a 2x binned map, blurred to ~0.5 mm so the score
        # surface is smooth enough for 5-degree steps.
        search_um = pixel_um * 2
        coarse = ndi.gaussian_filter(block_mean(density, 2), 430 / search_um)
        _, (r0, c0, a, s) = _grid_search(
            coarse, used, search_um,
            angles=np.arange(0, 360, 5), scales=np.arange(0.7, 1.31, 0.1),
        )
        r0, c0 = (r0 + 0.5) * 2 - 0.5, (c0 + 0.5) * 2 - 0.5
        # _grid_search offsets are relative to the centroid of `used`; convert
        # to the centroid of the full plan so all fits share one origin.
        r0, c0 = _recentre(r0, c0, a, s, pixel_um, used, plan)
    else:
        r0, c0, a, s = start.row0, start.col0, start.angle_deg, start.scale

    def to_used(r0, c0, a, s):
        return _recentre(r0, c0, a, s, pixel_um, plan, used)

    x0 = (*to_used(r0, c0, a, s), a, s)
    (ru, cu, a, s), score = _refine(density, used, pixel_um, x0)
    r0, c0 = _recentre(ru, cu, a, s, pixel_um, used, plan)
    return PlanFit(r0, c0, float(a % 360), float(s), pixel_um, plan, float(score))


def _recentre(r0, c0, a, s, pixel_um, src, dst):
    """Move the origin from the centroid of plan ``src`` to that of ``dst``."""
    csrc = np.mean(list(src.values()), axis=0)
    cdst = np.mean(list(dst.values()), axis=0)
    dr, dc = _mm_to_px(*(cdst - csrc), a, s, pixel_um)
    return r0 + dr, c0 + dc


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------


def disc(shape, centre, radius) -> tuple[tuple[slice, slice], np.ndarray]:
    """Bounding-box slices and boolean disc within them."""
    r, c = centre
    r0, r1 = int(max(0, np.floor(r - radius))), int(min(shape[0], np.ceil(r + radius) + 1))
    c0, c1 = int(max(0, np.floor(c - radius))), int(min(shape[1], np.ceil(c + radius) + 1))
    rr, cc = np.ogrid[r0:r1, c0:c1]
    return (slice(r0, r1), slice(c0, c1)), (rr - r) ** 2 + (cc - c) ** 2 <= radius**2


@dataclass
class SectionResult:
    section: str
    pixel_um: float
    full_fit: PlanFit
    slot_fits: dict[int, PlanFit]
    rows: list[dict]
    tissue: np.ndarray
    residuals: dict[str, Residual]

    def table(self):
        import pandas as pd

        return pd.DataFrame(self.rows)


def measure_slots(
    sec: Section,
    tissue: np.ndarray,
    residuals: dict[str, Residual],
    slot_fits: dict[int, PlanFit],
    full_fit: PlanFit,
    radius_mm: float = ROI_RADIUS_MM,
    k: float = THRESHOLD_K,
) -> list[dict]:
    """Coverage and intensity for every slot and channel.

    ``coverage`` is GFP+ tissue pixels / tissue pixels inside the ROI.
    ``coverage_k-1`` / ``coverage_k+1`` are the spec's threshold sensitivity.
    ``mean_residual`` is background-subtracted intensity per tissue pixel.
    ``tissue_fraction`` is how much of the ROI is tissue -- low means the ROI
    fell on a tear or off the section edge. ``loo_shift_mm`` is how far the
    slot moved when its own signal was left out of the fit.
    """
    rows = []
    full_pts = full_fit.points()
    for slot, fit in sorted(slot_fits.items()):
        centre = fit.points()[slot]
        radius = radius_mm * fit.scale * 1000 / sec.pixel_um
        box, d = disc(sec.shape, centre, radius)
        t = tissue[box] & d
        n_tissue = int(t.sum())
        shift = np.hypot(*(np.subtract(centre, full_pts[slot]))) * sec.pixel_um / 1000
        base = {
            "section": sec.name,
            "slot": slot,
            "row": round(centre[0], 1),
            "col": round(centre[1], 1),
            "radius_px": round(radius, 1),
            "tissue_fraction": round(n_tissue / max(int(d.sum()), 1), 3),
            "loo_shift_mm": round(float(shift), 3),
            "fit_angle_deg": round(fit.angle_deg, 1),
            "fit_scale": round(fit.scale, 3),
        }
        for ch, res in residuals.items():
            v = res.values[box][t]
            row = dict(base, channel=ch, noise_sd=round(res.noise_sd, 1), threshold_k=k)
            for kk, name in [(k, "coverage"), (k - 1, "coverage_k-1"), (k + 1, "coverage_k+1")]:
                row[name] = float(np.mean(v > kk * res.noise_sd)) if v.size else np.nan
            row["mean_residual"] = float(v.mean()) if v.size else np.nan
            rows.append(row)
    return rows


def analyse_section(
    sec: Section,
    plan: dict[int, tuple[float, float]] | None = None,
    radius_mm: float = ROI_RADIUS_MM,
    k: float = THRESHOLD_K,
    fit_factor: int = 8,
) -> SectionResult:
    """Tissue, background, plan placement, and measurement for one section."""
    plan = dict(plan or PLAN_MM)
    tissue = tissue_mask(sec)
    residuals = {ch: gfp_residual(sec, ch, tissue, k=k) for ch in GFP_CHANNELS}

    coarse_um = sec.pixel_um * fit_factor
    density = plan_density(
        [r.positive(k) for r in residuals.values()], tissue, fit_factor,
        sigma_px=250 / coarse_um,
    )
    full = fit_plan(density, coarse_um, plan)
    loo = {s: fit_plan(density, coarse_um, plan, leave_out=s, start=full) for s in plan}

    full_px = full.rescaled(fit_factor)
    loo_px = {s: f.rescaled(fit_factor) for s, f in loo.items()}
    rows = measure_slots(sec, tissue, residuals, loo_px, full_px, radius_mm, k)
    for r in rows:
        r["angle_flag"] = abs(((full.angle_deg - EXPECTED_ANGLE_DEG + 180) % 360) - 180) > ANGLE_WARN_DEG
    return SectionResult(sec.name, sec.pixel_um, full_px, loo_px, rows, tissue, residuals)


def slot_to_target(slot: int, mirrored: bool) -> int:
    """Target number for a slot, given the section's mounting.

    ``mirrored=False``: targets 1-3 lie in the image's lower hemisphere when
    anterior faces left. This has to come from outside the image -- a side
    mark, or the mounting convention -- because the plan is symmetric.
    """
    return MIRROR[slot] if mirrored else slot


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------


def plot_section(sec: Section, result: SectionResult, factor: int = 8):
    """Composite with tissue outline and ROIs, plus the GFP+ masks."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    def stretch(a, lo=1, hi=99.8):
        p1, p2 = np.percentile(a, [lo, hi])
        return np.clip((a - p1) / (p2 - p1 + 1e-9), 0, 1)

    cy5 = stretch(np.log1p(sec.binned("CY5", factor)))
    tissue = block_mean(result.tissue.astype(np.float32), factor) > 0.5
    masks = {
        ch: block_mean((r.positive() & result.tissue).astype(np.float32), factor)
        for ch, r in result.residuals.items()
    }
    table = result.table()

    fig, axes = plt.subplots(1, 3, figsize=(21, 7.5))
    panels = [
        ("NeuN (log) + GFP+ union", np.dstack([cy5 * 0.6, np.clip(masks["FITC"] + masks["TRITC"], 0, 1), cy5 * 0.0])),
        ("FITC (native) GFP+", masks["FITC"]),
        ("TRITC (stain) GFP+", masks["TRITC"]),
    ]
    for ax, (title, img) in zip(axes, panels):
        ax.imshow(img, cmap="gray", vmin=0, vmax=1 if img.ndim == 2 else None)
        ax.contour(tissue, [0.5], colors="c", linewidths=0.6)
        ch = title.split()[0]
        for slot, fit in result.slot_fits.items():
            r, c = fit.points()[slot]
            rad = ROI_RADIUS_MM * fit.scale * 1000 / sec.pixel_um
            ax.add_patch(Circle(((c + 0.5) / factor - 0.5, (r + 0.5) / factor - 0.5),
                                rad / factor, fill=False, ec="yellow", lw=1.4))
            label = f"S{slot}"
            if ch in masks:
                cov = table.query("slot == @slot and channel == @ch")["coverage"].iloc[0]
                label += f"\n{cov:.1%}"
            ax.text((c + 0.5) / factor - 0.5, (r + 0.5) / factor - 0.5, label,
                    color="yellow", ha="center", va="center", fontsize=9, weight="bold")
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
    f = result.full_fit
    fig.suptitle(
        f"{sec.name}   fit: angle {f.angle_deg:.0f} deg, scale {f.scale:.2f}   "
        f"ROI r = {ROI_RADIUS_MM} mm, GFP+ = residual > {THRESHOLD_K:g} SD   "
        "(slots assume targets 1-3 in the lower hemisphere)"
    )
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    import argparse

    import pandas as pd

    p = argparse.ArgumentParser(prog="python -m fus.histology", description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("export", help="write each section of a .vsi as OME-TIFF")
    e.add_argument("vsi", type=Path)
    e.add_argument("out_dir", type=Path)
    e.add_argument("--downsample", type=float, default=4)

    m = sub.add_parser("measure", help="GFP coverage per slot for exported sections")
    m.add_argument("sections", nargs="+", type=Path)
    m.add_argument("--csv", type=Path, required=True)
    m.add_argument("--figures", type=Path, metavar="DIR")
    m.add_argument("--radius-mm", type=float, default=ROI_RADIUS_MM)
    m.add_argument("-k", type=float, default=THRESHOLD_K)
    args = p.parse_args(argv)

    if args.cmd == "export":
        for path in export_sections(args.vsi, args.out_dir, args.downsample):
            print(f"wrote {path}")
        return 0

    tables = []
    for path in args.sections:
        sec = load_section(path)
        res = analyse_section(sec, radius_mm=args.radius_mm, k=args.k)
        tables.append(res.table())
        f = res.full_fit
        print(f"{sec.name}: angle {f.angle_deg:.0f} deg, scale {f.scale:.2f}")
        if args.figures:
            import matplotlib

            matplotlib.use("Agg")
            args.figures.mkdir(parents=True, exist_ok=True)
            fig = plot_section(sec, res)
            fig.savefig(args.figures / f"{sec.name}_slots.png", dpi=110)
        del sec, res
    out = pd.concat(tables, ignore_index=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.csv, index=False)
    print(f"wrote {len(out)} rows to {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
