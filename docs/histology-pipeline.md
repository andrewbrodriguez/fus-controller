# Histology pipeline — GFP coverage per target

**Status:** first pass, run on Mouse_01 only (2026-09-16). Every parameter marked
*provisional* still needs a decision from Nick or Bernie.

This is the ground-truth half of the project, implementing Phase 2 of
[`ground-truth-spec.md`](ground-truth-spec.md). It turns a whole-slide `.vsi` scan into
GFP area coverage for each of the six sonication targets, for every section. Code is in
[`src/fus/histology.py`](../src/fus/histology.py).

It measures area coverage (option A in
[decision 001](decisions/001-delivery-endpoint.md)). Everything up to the last step is
the same whichever endpoint is chosen.

---

## Run it

```bash
# 1. One OME-TIFF per section, 4x downsampled (1.3 µm/px). Needs QuPath 0.7;
#    set $QUPATH if it is not in /Applications. ~20 s per section.
PYTHONPATH=src python -m fus.histology export \
    data/histology/Mouse_01/Image.vsi data/processed/histology/Mouse_01/ds4

# 2. Coverage per slot, plus an overlay figure per section. ~11 s per section.
PYTHONPATH=src python -m fus.histology measure \
    data/processed/histology/Mouse_01/ds4/section_s*.ome.tif \
    --csv results/histology/mouse01_slots.csv --figures results/histology/figures
```

Geometry tests: `PYTHONPATH=src python -m pytest tests/`

---

## Geometry — two corrections to the spec

**`posx.xyz` is an axial layout.** The planning image in `Mouse_01_AAV_Images.pptx` is
the axial panel of a three-view MRI viewer, and the six circles fit `xyz` at the same
scale on both axes (271 k vs 288 k EMU/mm on the slide; every circle within 0.07 mm). So:

| Column | Meaning |
|---|---|
| `x` | medio-lateral — targets 1–3 on one side, 4–6 mirrored on the other |
| `y` | antero-posterior, +y anterior |
| `z` | depth — **0 for all six targets** |

Spec §3 reads `y` as dorsal/ventral, which is wrong. Targets 1/4 are the anterior pair,
2/5 the lateral pair, and 3/6 the posterior-medial pair, 0.5 mm behind 2/5.

**The Mouse_01 sections are horizontal, not coronal.** Every section contains both
hippocampi and the cerebellum. Anterior faces the image left. The targets share one
depth, so **a section at that depth crosses all six**, and no section has to be matched
to a target by AP position. (Spec §4 and the README banner caption say "coronal".)

**The image cannot tell the hemispheres apart.** The plan is mirror-symmetric (target
k ↔ k ± 3), and the sections were evidently mounted both ways up (see below). Results
are therefore reported per **slot**: the plan position *assuming targets 1–3 lie in the
image's lower hemisphere when anterior faces left*. `slot_to_target(slot, mirrored)`
converts a slot to a target once a section's orientation is known.

---

## Method

| Step | How | Parameter |
|---|---|---|
| Resolution | QuPath `convert-ome`, pyramid level 4× | 1.3 µm/px |
| Tissue mask | Otsu on log(CY5) + log(TRITC) at 5.2 µm/px. DAPI is not used because out-of-focus haze makes the glass beside some sections as bright as tissue | holes < 0.05 mm² filled; fragments < 0.5 mm² dropped |
| Background | Normalised Gaussian average of tissue, recomputed 3× with pixels above threshold excluded. Handles regional autofluorescence (s5 has two background levels) | σ = 1 mm *(provisional)* |
| GFP+ | residual > k × robust pixel-noise SD, per channel | k = 5, k ± 1 also reported *(provisional)* |
| Plan placement | Similarity fit (rotation, shift, shrinkage) of the six-point plan to the GFP+ density map, using both channels combined. Full search, then a **leave-one-out** refinement per target | scale 0.7–1.3 |
| ROI | Disc at each target's leave-one-out position | r = 0.75 mm × fitted scale *(provisional)* |
| Measure | coverage = GFP+ tissue px / tissue px in ROI; mean background-subtracted intensity | both channels |

**Why leave-one-out.** Placing ROIs from the GFP itself would put each ROI on its own
plume and inflate its coverage. Instead, target k's ROI comes from a fit to the other
five, so its own signal plays no part in where it lands. An empty target is measured
where the rest of the pattern says it should be.

**Why not the no-FUS target as the threshold reference** (spec §8)? That site is also
the only check that the pipeline reads zero when nothing was delivered, and it can't do
both jobs. Whole-section robust statistics are used instead. Revisit once there is an
untreated animal or a second negative region.

---

## Mouse_01 results

![Section s4: GFP+ masks and the six ROIs](figures/histology-s4-slots.png)

Section s4. Left: NeuN (log) with the GFP+ union. Middle and right: GFP+ masks per
channel. Slot 3 reads 0.0 % (FITC) and 1.2 % (TRITC). Slot 6 is 83 % in FITC and 3.6 %
in TRITC — see known problem 1.

### Per section

| Section | Fit angle | Fit scale | Empty slot | Notes |
|---|---|---|---|---|
| s2 | 107° | 0.89 | 3 | |
| s3 | 217° | 0.70 | — | **Fit failed** (half a section; angle flagged, scale at bound). Excluded |
| s4 | 91° | 0.82 | 3 | |
| s5 | 79° | 0.91 | 6 | |
| s6 | 71° | 0.84 | 6 | |

Scale ≈ 0.8–0.9 is tissue shrinkage relative to the plan, which is plausible for fixed,
sectioned tissue.

### Orientation

The empty slot switches between 3 and 6, so the sections were not all mounted the same
way. To test this without the empty pair deciding its own label, each combination of
mirrorings was scored by how well sections agree on **targets 1, 2, 4, 5 only**. The best
is s2 and s4 as-is, s5 and s6 mirrored: mean between-section SD 0.10, against 0.14 for
the runner-up. Under that orientation the empty site is the **same target in all four
sections**. That agreement is independent evidence.

Which target it is — 3 or 6 — depends on the global flip, which the image can't give.
The table below follows the deck ("No FUS control" at position 3).

### Coverage by target (s2, s4, s5, s6)

| Target | Deck label | TRITC mean (range) | FITC mean (range) |
|---|---|---|---|
| 1 | 420 bursts, 0.85 | 0.82 (0.66–0.90) | 0.41 (0.23–0.53) |
| 2 | 45 bursts, 1.05 | 0.65 (0.59–0.77) | 0.36 (0.27–0.46) |
| 3 | No FUS control | **0.02** (0.01–0.05) | **0.00** |
| 4 | 90 bursts, 0.95 | 0.52 (0.48–0.60) | 0.26 (0.15–0.45) |
| 5 | 60 bursts, 0.85 | 0.31 (0.00–0.42) | 0.30 (0.27–0.32) |
| 6 | 240 bursts, 0.55 | 0.19 (0.00–0.57) | 0.66 (0.17–0.84) |

**Not a dose-response.** This is one animal. The deck labels are position labels, and
their correspondence to recordings is still unresolved (spec §3; the spreadsheet's
`Brain Region` column disagrees with the deck for positions 4–6).

---

How these numbers were joined to acoustic dose and plotted:
[mouse01-dose-delivery.md](mouse01-dose-delivery.md).

---

## Known problems

1. **FITC and TRITC disagree at target 6** (and at target 5 in s5). The site is bright in
   native GFP and near background in the antibody channel. At full resolution in s4, the
   99th-percentile brightness is 21.9 k (FITC) against 1.7 k (TRITC); at target 1 TRITC
   reaches 22.8 k. FITC also marks a detached cerebellum fragment in s2 that TRITC does
   not, so FITC picks up at least some non-GFP signal. **Needs a look at full
   resolution:** s4 is series 4 of `Image.vsi`, and the site is centred near
   x = 11 070, y = 14 400 px.
   The combined map used for placement includes FITC, so in s5 this signal probably pulls
   slot 2 off target.
2. **Tile seams** are visible at full resolution and are not corrected.
3. **The tissue mask accepts bright non-tissue** — a strip beside s3 and an object in the
   corner of s5. It only matters if an ROI lands there.
4. **Placement depends on the other plumes.** If another target's plume sits off its
   planned spot, the whole pattern shifts. `loo_shift_mm` > ~0.5 flags this (s4 slot 5:
   0.76; s6 slot 2: 0.67).

## Output columns

| Column | Meaning |
|---|---|
| `section`, `slot`, `channel` | One row per section × slot × channel |
| `coverage`, `coverage_k-1`, `coverage_k+1` | GFP+ fraction of tissue in the ROI at k and k ± 1 |
| `mean_residual` | Background-subtracted intensity per tissue pixel |
| `tissue_fraction` | Share of the ROI that is tissue — low means a tear or the section edge |
| `loo_shift_mm` | How far the ROI moved when its own target was left out of the fit |
| `fit_angle_deg`, `fit_scale` | That slot's leave-one-out transform |
| `angle_flag` | Full fit more than 30° from anterior-left — treat the section as failed |
| `row`, `col`, `radius_px`, `noise_sd`, `threshold_k` | Provenance, in 4× export pixels |

## Open questions

1. **Which image side is the animal's right?** Is there a side mark, or a mounting
   convention? Were sections mounted both ways up (s2/s4 vs s5/s6)?
2. **FITC vs TRITC** — which is primary, and what is the FITC-only signal at target 6?
3. **ROI radius** — what is the focal spot size at 837 kHz?
4. **Section order and depth** — series 2–6 need not be in depth order, and the spacing
   is unknown.
