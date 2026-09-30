# Where the project is — 2026-09-30

Week ~5 of 13. Mouse 2 histology arrived 2026-09-29; on 2026-09-30 the ROI radius moved to
1.0 mm and Mouse 2 was exported and measured. This file is the moment-in-time picture; the
durable method write-ups are in `docs/`.

---

## One-paragraph summary

Both halves of the pipeline are built and verified on the single imaged animal: the acoustic
side reproduces the lab's own numbers exactly, and the tissue side turns a whole-slide scan
into a GFP coverage number per sonication target. The bookkeeping that blocked everything —
which recording was fired at which target — was resolved with Nick on 2026-09-18. Mouse 2
is the first clean animal (TargetN = position N, all six sonicated, recordings already
local). Its sections are measured per slot but **not yet mapped to targets**. That needs a
person to click each section's front and notch in `python -m fus.orientation`. Two of 24 brains are
imaged, and one of four acoustic sessions is synced.

## What exists and is trusted

| Piece | State |
|---|---|
| `src/fus/extract.py` | Python port of `nt_ExtractHarmonicData.m`. Reproduces the sheet's `Cumulative 2nd Harmonic`, `Mean Voltage`, `Cumulative Voltage` exactly for all 19 recordings checked (mice 1–3) |
| `src/fus/histology.py` | `.vsi` → per-target GFP coverage. ~20 s/section to export, ~11 s/section to measure. ROI r = 1.0 mm. `export --prefix` for multi-slide animals; `measure --orientation` takes each section's clicked front |
| `src/fus/orientation.py` | napari: `annotate` (click centre, front, notch side) and `review` (circles labelled by target on GFP, ✓/✗). Saves to `data/section_orientation.csv` |
| `notebooks/ingest_new_histology.ipynb` | The per-animal workflow: export → orient → finetune (optional) → measure → review → coverage per target |
| `src/fus/rois.py` + `data/roi_locations.csv` | Hand-finetuned T1–T6 shapes (drag/resize/rotate in napari), measured exactly as drawn |
| `tests/test_histology.py` | 13 synthetic geometry tests, passing |
| `data/section_orientation.csv` | Clicked centre, front, notch side, and review verdict per section (4x-export px). All 16 sections clicked 2026-09-30; **not yet reviewed** |
| `scripts/mouse02_dose_vs_coverage.py` | Mouse 2 join. Stops until the orientation sheet is filled and `measure` re-run with it |
| `notebooks/01_mouse01_target_regions.ipynb` | Opens the sections, draws labelled targets, optional napari viewer |
| `scripts/mouse01_dose_vs_coverage.py` | Joins dose to delivery, writes the CSV and both plots |
| `scripts/mouse01_example_figure.py` | The README example figure |
| Docs | `docs/histology-pipeline.md` (method), `docs/mouse01-dose-delivery.md` (how the plots were made, with assumptions), `docs/acoustic-extraction.md`, `docs/ground-truth-spec.md`, `docs/decisions/001-delivery-endpoint.md` |

**Mouse 1 result (TRITC coverage, 1.0 mm ROI, mean of 4 sections):** T1 0.66, T2 0.48,
T3 0.02 (no FUS), T4 0.48, T5 0.25, T6 0.18. Dose vs delivery gives r = 0.85 across six
targets, but that is carried by the control and by T1 (three sonications); among the four
single-sonication targets r = 0.12. Treat it as proof the pipeline works, not as a result.
(At the old 0.75 mm radius: 0.82 / 0.65 / 0.02 / 0.52 / 0.31 / 0.19, r = 0.83 — same
ranking.)

**Mouse 2 doses** (from the recordings; match the sheet): T1 0.62 (135 @ 0.65), T2 1.02
(105 @ 0.85), T3 0.77 (75 @ 0.95), T4 1.07 (135 @ 0.75), T5 0.37 (75 @ 0.65), T6 1.60
(255 @ 0.85, interlock fired).

## Settled with Nick (2026-09-18)

- **Recording → position.** `TargetN` was fired at position N for every animal except Mouse 1,
  where the first three sonications all landed at position 1 (420 bursts total), then 4, 5, 6,
  then a 7th back at position 2. Position 3 was never sonicated.
- **Orientation.** Sections are notched: right on Mouse 1, left from Mouse 2 on. Nick first
  wrote "bottom right" for Mouse 2 (9/29 11:59) and corrected it to left the same day,
  saying left "makes more sense with the delivery and harmonic doses". That reasoning is
  circular for us, so check the notch in the images.
- **Delivery measure.** The anti-GFP stain (TRITC), not native GFP.
- **Focal spot.** ~2 × 2 mm in x/y, ~3 mm in z. Precise dimensions still to come.
- **Feature priority.** Cumulative 2nd-harmonic AUC first; other bands are a secondary study.
- **Translation framing.** 4–5 minutes of sonication per target is impractical in a patient and
  long exposure carries risk — that's the argument for predicting delivery rather than
  applying more energy.

## New from Nick (2026-09-18 → 09-29), not yet in the repo

- **Mouse 2 summary PDF** (9/29, the updated version assumes the left notch). His
  per-target reading is the natural cross-check for ours. Save to
  `docs/lab/Mouse_02_summary.pdf` and run `pdftext` on it to make a `.pdf.md` sidecar.
- **Brain template + target mask** (9/18): 6 targets as 2 × 2 × 2 mm spheres, 0.1 mm voxels.
  Could replace the hand-coded `PLAN_MM` geometry, or check it. Save to `docs/lab/`.
- **Head roll in Mouse 2:** right-side targets sat higher, left-side lower. One horizontal
  section therefore cuts the two sides at different depths (assumption A8).

## Open — blocking or slowing

1. **Mouse 2 placement review.** Placement is now a rigid template from the clicks (centre,
   front, notch side), not a GFP fit. See `docs/histology-pipeline.md#mouse_02`. All
   sections are clicked, and Mouse 2 is re-measured with the template. Verdicts given on
   the old fit placements were moved into `notes` (10 wrong, 2 ok) and cleared. Next:
   open `notebooks/ingest_new_histology.ipynb`, check the centres in the orientation window
   (the template is drawn live; 5 Mouse 2 centres are still the automatic ones), re-run the
   measure cell, mark ✓/✗ in review, then run `scripts/mouse02_dose_vs_coverage.py`.
   **First, benchmark on Mouse 1:** run the notebook with `MOUSE = 1`, place the centres
   by hand, then `scripts/benchmark_mouse01.py`. With the automatic centres it fails on
   placement (median offset 0.61 mm), delivery (max difference 0.14) and ranking (ρ 0.83).
   The Mouse 2 template numbers are provisional until it passes.
2. **Tissue throughput.** Two brains imaged. Everything downstream is gated on this.
3. **Capsid assignment.** The study is AAV9 ×12 and AAV.CPP16 ×12, but no file records which
   mouse is which. Any delivery model needs it as a factor. Ask Nick.
4. **Remaining acoustic sessions.** Three of four days are still on Dropbox, unsynced.
5. **The FITC-only signal at target 6.** Bright in native GFP, near background in the stain.
   Unexplained. Someone should look at s4 (series 4 of `Image.vsi`) near x = 11 070,
   y = 14 400 px at full resolution.
6. **Section order, spacing, and which sections fall inside the 3 mm focal column.**
7. **Imaging exposure settings** — fixed across sessions or not; between-animal intensity
   comparisons depend on it.

## Findings that exist only here

Computed from `data/Mouse_Controller_Data.xlsx` (all 24 mice) but **not yet scripted or
written into `docs/`** — re-derive before citing:

- Every mouse has the same reference target, 120 bursts at goal 0.75. Across animals the
  cumulative 2nd harmonic ranges 0.70–1.73 (21% CV) and mean voltage 0.067–0.116 V (14% CV).
  That is the animal-to-animal spread at identical settings, and it gives a per-animal
  normaliser.
- After removing the N × goal trend, only ~21% of the remaining variance is animal-level;
  ~79% is target-level. Argues for within-animal contrasts over a per-animal offset alone.
- `Mean MPa` = 4.05 × `Mean Voltage` exactly, and `Mean MI` = MPa / √0.837. Both are fixed
  conversions with no skull correction — do not use them as independent predictors.
- Dose→emission slope by region runs 88 (region 2) to 183 (region 6).
- Stored burst count = programmed + 15 for mice 1–6, + 12 for mice 7–24, which identifies the
  first session.

## Traps for whoever picks this up

- **The spreadsheet's Mouse 1 rows 3–6 are rotated by one run.** `Brain Region` is right;
  the values in those rows belong to different sonications. Read doses from the recordings.
- **Slot ≠ target.** The pipeline reports plan slots; converting to targets needs the
  section's mirroring (`MIRRORED` in `scripts/mouse01_dose_vs_coverage.py`). For new animals
  it comes from `data/section_orientation.csv`.
- **Template runs write `mouseNN_template_slots.csv`**; `mouse01_slots.csv` is the GFP-fit run every Mouse 1 number comes from. Don't overwrite it.
- **Several Mouse 2 sections have little GFP anywhere** (slide01_s2, slide01_s4,
  slide04_s2, slide03_s4), probably outside the focal column. Consider `exclude=yes`
  once the order and depth are known.
- **`docs/figures/histology-sections-overview.png` still shows 0.75 mm circles.** No script
  generates it; either write one or drop the figure. The notebook's saved outputs are
  also from the 0.75 mm run.
- **`data/histology/Mouse_02.zip` (17.8 GB)** is still there after extraction; delete it
  once you're happy with the extracted copy.
- **Section s3 is half a section**; its plan fit fails and the code flags it. It is excluded
  from every number quoted.
- **Use `.venv/bin/python`.** System `python3` is 3.9 and cannot run this; the system 3.12 has
  older scikit-image. `histology.tissue_mask` deliberately avoids APIs that differ between
  scikit-image 0.24 and 0.26.
- **Figures are generated**, not hand-made. Re-run the scripts rather than editing images.
- **Agents here are text-only** (see `AGENTS.md`): no reading PDFs or figures directly, so any
  regenerated figure needs a human to eyeball it.

## Next steps, in the order I would do them

1. Record Mouse 2 (and Mouse 1 notch) orientation, re-measure, run the Mouse 2 join, and
   compare target by target with Nick's summary PDF.
2. Script the spreadsheet EDA above into `scripts/` or a notebook so the numbers are
   reproducible, and write it up as `docs/controller-data-eda.md`.
3. Chase the capsid list (Mouse 1 and 2 may be in different cohorts, which matters as soon
   as they're pooled) and the remaining acoustic sessions.
4. When new tissue arrives: export (`--prefix` per slide), record orientation, measure,
   join.
5. Later, once results hold up: port the feature extraction to MATLAB for the lab's workflow.

## Deliverable and timing

Draft journal manuscript by end of reading period (week 13). On current data that report is a
validated pipeline plus a preliminary relationship, not a trained predictive model. If only
three or four animals are imaged by December, say so in the report's design section rather
than presenting it as a shortfall.
