# Where the project is — updated 2026-11-11

Week ~9 of 13. This file is the moment-in-time picture; the durable method write-ups are in
`docs/`.

---

## One-paragraph summary

Both halves of the pipeline work. The acoustic side reproduces the lab's own numbers, and
the tissue side turns a whole-slide scan into GFP coverage per target. **Two brains are
measured.** All 16 usable sections of Mouse 1 and Mouse 2 were oriented by hand (centre,
front and notch side, plus finetuned target circles in napari) and joined to acoustic dose.
Delivery is measured two ways in the same ROIs: **pixel coverage (pipeline A)** and the
**fraction of NeuN cells tagged GFP+ (pipeline B**, StarDist; GFP+ above 4.05× the slice's
background). They agree at r = 0.98 across 96 ROIs. Within each slice, higher-dose targets
tend to carry more GFP in both animals: positive rank correlation in 11/11 slices with clear
signal under A (median ρ +0.54) and 10/11 under B (+0.61). Neither measure has been checked
against a hand count, and the review verdicts still need re-tagging. Tissue remains the
limit: 2 of 24 brains imaged. **Four more are sectioned and stained and are being imaged Mon
10/26.** A call with the lab on normalising delivery across
animals was held 10/28. Meeting notes and email drafts are kept locally (git-ignored).

## Since 11/04

- **Lab update 11/09** (correspondence kept locally in `docs/lab/correspondence.md`). **Meeting Wed 11/18, 11:15.**
  - **All imaging is done**, except the missing slide for Mouse 5 well 11 slide 1.
  - **Not yet downloaded:** more Mouse 5 scans and updated IF decks were uploaded. Locally,
    `Mouse_05/` holds only well 10 and the Mouse 3–5 decks are still the 11/01 versions.
    Fetch both from Dropbox.
  - **The updated decks are arranged roughly dorsal to ventral.**
    That is a section order from the lab, not from the GFP. It could replace the assumed
    filename order, but the decks are images only, so a person has to match them to the scans.
  - **Virus:** mice 1–12 got AAV9 and 13–24 AAV.CPP16. Every imaged mouse (1–6) is AAV9.
  - **Analyses the lab suggests:**
    - by region: the thalamus (T3/T6) looks much stronger;
    - laterality: right anterior striatum (T1) looks lower than left (T4), which fits T1
      reading low in Mice 3–5;
    - same cumulative dose reached by different parameters: few bursts at high intensity
      vs many bursts at low intensity.
- **Region analysis, first pass (11/13)** (`scripts/region_delivery.py` →
  `results/histology/region/README.md`, Mice 1–5, 28 single-sonication targets). Delivery is
  compared after taking out dose (logit delivery vs log dose, one offset per mouse);
  p-values come from shuffling targets within each mouse.
  - **Thalamus is strongest raw, but not at matched dose.** T3/T6 got the most dose in these
    four mice by chance (position 3: 1.44 in mice 2–5 vs 1.29 over mice 2–24; dose doesn't
    differ by position overall, p 0.76). Dose-adjusted, thalamus vs the rest is +0.21, p 0.62.
  - **T1 < T4 holds at matched dose**: lower in 4/4 mice on both measures (cells p 0.018,
    pixels p 0.034), with equal mean doses (0.80 vs 0.81). The other right–left pairs go both
    ways, so it isn't the whole right side.
  - **Not independent:** the lab suggested both after viewing these animals. Mouse 6 is the
    first clean test. T1 is also always fired first, so side and firing order are confounded.
- **Mouse 6 scans arrived (11/11).** `Mouse_Cntr_06.zip` was extracted to
  `data/histology/Mouse_06/` (all 32 files CRC-match the zip; the zip was then deleted). **New format again: one `.vsi`
  per slide with several sections in it** (10.1 and 11.1: 4 each; 10.2 and 11.2: 2 each; 12
  total), no Overview scans. Ingest like Mouse 2 (`--prefix w<well>sl<slide>`); channels not yet
  checked. Not yet exported, clicked or measured.
- **New lab files, in `docs/lab/`** (both from the lab; see its README):
  - `AAV_By_TreatmentLocation.pptx`: bursts, MPa and dose per target for Mice 1–6. Mouse 6's
    doses match our recordings (T1 1.04, T2 0.47, T3 2.02, T4 0.92, T5 0.43, T6 3.33).
  - `IF_Images_Table.xlsx`: per-section imaging QC for Mice 2–6. This likely explains the
    unexplained exclusions in Mice 3–5 (e.g. M3 11.2.3 "too torn", M5 10.1.3 "something wrong
    with TRITC"), **but its section counts exceed our scans**: it lists 10 / 9 / 7 imaged
    sections for Mice 3 / 4 / 5 against our 7 / 5 / 3 scans. Mouse 5 well 11 slide 2 (sections
    1–2) is marked imaged but isn't on disk; more Mouse 5 scans were uploaded 11/09, not yet
    downloaded. Mapping table sections to our scan numbers is not yet done.
- `Mouse_06_IF_Images.pptx` in `data/histology/` was replaced by a 12-slide version.

## Since 11/01

- **Repo simplified (11/04): aggregate figures only.** Removed ~60 per-mouse, per-section and
  per-session images (`results/histology/figures/`, `aligned/`, the Mouse 1 and Mouse 2
  plots, per-slice heatmaps and panels, per-session controller traces) and
  `scripts/aligned_sections.py`. The scripts no longer draw them; the numbers stay in the
  CSVs. The README headline is now the five-mouse pooled figure
  (`results/histology/slices/cumulative_pooled.png`). Old images are in git history.

## Since 10/28 (updated 2026-11-01)

- **Mice 3, 4, 5 imaged and measured.** Raw scans are in `data/histology/Mouse_0{3,4,5}/`, and
  exports in `data/processed/histology/.../ds4`. The new format is one section per `.vsi`, named
  `Image_Cntr_<mouse>.<well>.<slide>_<scan>` (export prefix `w<well>sl<slide>_<scan>`).
  There are **three channels, CY5/TRITC/FITC, with no DAPI**: `export_sections` and `fus.cells`
  now pick channels by name, not by position.
- **Section counts:** Mouse 3 has 7 (2 excluded), Mouse 4 has 5 (1 excluded) and Mouse 5 has
  3 (1 excluded). The excluded sections have no note saying why. Andrew clicked, finetuned
  and reviewed all of them (✓). Mouse 5 has well 10 only. Every QuPath reader finds just these
  three scans. More Mouse 5 data was uploaded 11/09 (see above).
  Mouse 6's zip was empty (re-sent 11/11; see above). Mouse 3's `_Image_Cntr_3.11.2_03_/` also holds a
  stray `_Image_Cntr 3.10.2_01_` stack and a `.vsi` that won't open.
- **Pipelines A and B were run on all 11 new sections; `scripts/slices_dose_delivery.py` now covers
  Mice 1–5** (`results/histology/slices/README.md`):
  - Within-slice ρ with dose is positive in 17/17 clear slices under A and 16/17 under B
    (median +0.60 and +0.61), and in 25/27 and 24/27 counting every slice.
  - Pooled over 30 targets: ΣGFP+/ΣNeuN r = +0.66, area ÷ animal mean +0.61, cells ÷ animal
    mean +0.56.
  - Per mouse, ΣGFP+/ΣNeuN r runs from +0.64 to +0.85.
  - Mouse 4 is the weakest: area r = +0.45 (its T3, dose 2.59, is high, but T2 is high in two
    T2-only sections).
  - **The cross-animal metric is still not decided.** The pooled figure shows all three side by
    side.
- **Coverage varies a lot between sections, most likely because of cut depth.** Sections with more
  tissue carry more GFP. Some sections show only T2, in Mouse 4 and in Mouse 2 alike;
  that suggests the lateral plumes reach higher, into cortex. This is a hypothesis, not tested.
- **Meeting 10/28:** the later mice got a different virus (resolved 11/09: mice 13–24 got
  AAV.CPP16).

## Since 10/20

- **Lab update 10/23**:
  - **Four brains** are being imaged 10/26.
  - **Sections:** 30 µm thick, ~360 µm apart (to be confirmed with Josh). Mouse 2 is probably
    two offset series of 6, some sections adjacent.
  - **Metrics:** the lab's preferred outcome measures are GFP+ cells summed over sections, then
    area above a per-mouse threshold (non-GFP ROI mean + 1–2 SD), then background-normalised
    mean intensity.
  - **Ceiling hits:** these are across the four surgery days, not time within a session.
- **Checked 10/25** (not yet scripted; listed under "Findings that exist only here"):
  - Raw totals don't compare across mice.
  - ΣGFP+ ÷ ΣNeuN needs no per-animal rescale and correlates best: r = +0.77 pooled.
  - Ceiling hits are mostly at the lateral targets.
  - NeuN profiles hint at three adjacent pairs in Mouse 2.

## Since 9/30

- **Meeting 10/2.** Delivery is now **summed over all slices of an animal**
  instead of normalised per slice. Open questions: how to choose which
  horizontal cut to analyse, and whether cuts are aligned or consistent across mice.
- **Cumulative analysis** (`scripts/slices_dose_delivery.py`, `results/histology/slices/README.md`).
  Each target's total GFP+ area (A) or GFP+ cells (B) is summed over every non-excluded
  slice, low-signal ones included, then divided by its animal's mean target. Pooled over 12
  targets: A r = +0.73 (p 0.0065), B r = +0.68 (p 0.016). Without the no-FUS control: A
  +0.62, B +0.57. The p-values are optimistic, because the 12 points come from 2 animals.
- **Acoustic features vs delivery** (`src/fus/features.py`, `scripts/acoustic_features_vs_delivery.py`,
  `results/acoustic_features/README.md`). This compared 15 pre-registered features with
  2f AUC alone, scored by leave-one-target-out predictions on 11 labelled targets. No
  multi-feature model beats AUC convincingly. Log-ridge on A has Q² 0.30 vs AUC's 0.06,
  but the gain has permutation p 0.23, and leaving out a whole mouse fails. Under B, AUC
  alone is best (Q² 0.15, permutation p 0.04). Treat the feature screen as hypotheses
  for the next animals.
- **Housekeeping (10/15).** Deleted `data/histology/Mouse_02.zip`, after checking that all
  32 files are extracted at matching sizes, and the full-resolution exports
  `Mouse_02/full/` (`slide01_s2`, `slide04_s3`). `view_slice.ipynb` re-exports on demand.
  The ROI size is quoted two ways, and both are right: r = 1.0 *template* mm
  (`ROI_RADIUS_MM`) × 0.87 shrinkage = 0.87 mm radius, **1.74 mm diameter on the tissue**.
  95 of 96 ROIs kept the default size; finetuning moved centres, not radii.
- **3D reconstruction of Mouse 2 (10/15)** (`src/fus/reconstruct.py`; GIF:
  `scripts/mouse02_reconstruction_3d.py`; interactive: `scripts/view_reconstruction_3d.py`;
  spacing CSVs in `results/histology/reconstruction/`; the GIF and PNGs are no longer committed). The 12 pipeline-A maps are aligned by
  the clicks, stacked, and Gaussian-blended in depth. Alignment checks out: per-ROI coverage from
  the resampled maps matches the table (r = 1.00, max diff 0.01), and each target's ROI centre
  varies by 0.06–0.21 mm SD across sections. **Spacing is inferred, ~93 µm** (73–139 µm over
  shrinkage 0.82–0.91). Each section's tissue area is matched to Allen CCF horizontal-section
  area: 47–69 mm² on tissue maps to 1.7–2.7 mm below the dorsal surface, all on the rising
  dorsal side, consistent with both hippocampi and the cerebellum being present. So the 12
  sections cover only about 1 mm of the ~3 mm focal column. Tears shrink area, so the
  span is an upper bound. In tissue units the estimate is ~84 µm, but **the
  default spacing is 250 µm** (Andrew thinks the sections are farther apart; 12 × 0.25 mm ≈
  the 3 mm focal column). GFP depth off the planes is inferred by treating each positive
  region as a sphere through its middle (`reconstruct.extrude`), which overstates depth
  for sections that only grazed a plume. Anatomy is one opaque Allen horizontal section
  (≈2.3 mm deep, fitted to the tissue outlines: scale 0.90, IoU 0.79), drawn under the stack.
  A **dose-sphere layer** puts one sphere per target at mid-stack, sized by the 2f AUC
  (volume ∝ dose, max r = 1 mm; `reconstruct.dose_spheres`), next to the GFP for comparison.
  The GIF is now a dashboard. Beside the rotating 3D view are a top-down panel of the
  *measured* GFP (max over sections, no depth inference) with dose circles, and bars of
  cumulative delivery against dose share per target (Mouse 2: r = +0.79 A, +0.64 B,
  n = 6). **The order is still assumed** (filename order), and in that order
  coverage jumps between slices (T6: .01 .87 .37 .02 .94 …). Don't choose the order from the GFP.

## What exists and is trusted

| Piece | State |
|---|---|
| `src/fus/extract.py` | Python port of `nt_ExtractHarmonicData.m`. Run on 19 recordings from mice 1–3; all 15 that can be matched to a spreadsheet row unambiguously match it exactly (`docs/acoustic-extraction.md`) |
| `src/fus/histology.py` | `.vsi` → per-target GFP coverage. ROI r = 1.0 template mm (1.74 mm diameter on tissue). Placement by GFP fit (Mouse 1 reference), clicked template, or finetuned shapes (`analyse_section(..., template=, rois=)`). `export --prefix/--series` |
| `src/fus/orientation.py` | napari: `annotate` (click centre, front, notch side; the view straightens to front-up, notch-left) and `review` (✓/✗). Each window writes only its own fields, so two open windows can't overwrite each other |
| `src/fus/rois.py` + `data/roi_locations.csv` | Hand-finetuned T1–T6 shapes per section, measured exactly as drawn. **All 16 sections finetuned** |
| `data/section_orientation.csv` | Clicks for all 17 sections (Mouse 1 s3 excluded). **No review verdicts**: the ones given were lost to the overwrite bug, since fixed |
| `notebooks/ingest_new_histology.ipynb` | Per-animal workflow: export → orient → finetune → measure → review → coverage per target |
| `notebooks/view_slice.ipynb` | Any section at full resolution (0.325 µm/px) in napari with ROIs. Exports on first open: ~7 min, ~7 GB per section |
| `scripts/slices_dose_delivery.py` | Dose vs coverage per slice and pooled; within-slice z-scores. Output and README in `results/histology/slices/` |
| `scripts/benchmark_mouse01.py` | Scores a placement method against Mouse 1's answer key (control, known front, validated fit) |
| `src/fus/cells.py` | **Pipeline B**: per-ROI cell count and fraction GFP+ (cell mean anti-GFP > 4.05× the slice's background). **Run on all 16 slices** (2026-09-30); ~1–1.5 min per slice with one QuPath launch and parallel StarDist on the GPU env (`segmentation_alpha/.venv-stardist-gpu`). One slice in `notebooks/cell_pipeline_one_slice.ipynb` |
| `segmentation_alpha/` | Cell-level pilot (crops, model comparison, threshold derivation, figures) and the model environments. Its own README |
| `tests/` | 33 tests, passing |

## Results so far

**Mouse 1 (GFP fit, 1.0 mm ROI, mean of 4 sections):** T1 0.66, T2 0.48, T3 0.02 (no FUS),
T4 0.48, T5 0.25, T6 0.18. Across 6 targets r = 0.85, but that leans on the control and on
T1 (three sonications); among the four single-sonication targets r = 0.12.

**Benchmark of hand placement on Mouse 1** (`results/histology/mouse01_benchmark.csv`):
finetuned shapes pass 4 of 5 checks. Placement is within a median 0.33 mm of the fit, and
the ranking against the fit is ρ = 0.94. Delivery fails on one target only: T4 reads 0.14
above the fit, against a tolerance of 0.10. With the automatic centres, before finetuning,
it failed 3 of 5.

**Both mice, per slice** (`results/histology/slices/`): within-slice Spearman of coverage vs
dose is positive in 11/11 slices with clear signal (median +0.54): Mouse 1 4/4, Mouse 2 7/7.
Five Mouse 2 slices are low-signal (mean coverage < 0.2, a cutoff chosen after looking),
probably cut outside the focal column. Mouse 2 T2 is the main exception: mid dose, low
coverage. It's a right-side lateral target, and the head was rolled with the right side
higher (lab, 9/29).

**Both pipelines, all 16 slices** (`results/histology/slices/README.md`): within-slice ρ with
dose is positive in 11/11 clear slices for pixels (A, median +0.54) and 10/11 for cells (B, median
+0.61). The one B exception is Mouse 2 `slide04_s4`: +0.09 under A, −0.09 under B. Across the 96
ROIs, B against A gives r = 0.98. Mouse 1's control T3 reads 0.7–1.3% GFP+ under B. The first
production threshold (0.0575 on a percentile scale) gave 5–14% there, because Mouse 1's slice
backgrounds (435–1,927 counts) are several times Mouse 2's (275); B now divides by each
slice's own background.

**Cell-level pilot** (`docs/gfp-cell-tagging.md`):
- **Model choice:** on 500 µm of NeuN, StarDist finds 428 cells in 2 s and Cellpose-SAM 308
  in 22 s; they agree on 271. StarDist was chosen for speed.
- **3 mm crop around Mouse 2 slide04_s3 T3:** 16,723 cells, 28% GFP+ at the log-Otsu
  threshold. That's 97–100% within 0.5 mm of the target centre and under 1% beyond 1.25 mm.

**Mouse 2 doses** (from the recordings; match the sheet): T1 0.62, T2 1.02, T3 0.77, T4 1.07,
T5 0.37, T6 1.60 (interlock fired).

## Settled with the lab

- **Recording → position (9/18).** `TargetN` = position N for every animal except Mouse 1 (first
  three sonications all at position 1, then 4, 5, 6, then a 7th at 2; position 3 never
  sonicated).
- **Notch (9/18, 9/29).** Right hemisphere on Mouse 1, left from Mouse 2 on. Our clicks agree:
  Mouse 1's notch clicks plus its empty target 3 put targets 1–3 on the animal's right,
  consistently across all four sections.
- **Delivery measure (9/18).** The anti-GFP stain (TRITC), not native GFP.
- **Focal spot (9/18).** ~2 × 2 mm in x/y, ~3 mm in z; precise dimensions to come.
- **Feature priority (9/18).** Cumulative 2nd-harmonic AUC first.
- **Answered 10/19** (Andrew, plus the Mouse 2 summary PDF, now `docs/lab/Mouse_02_Summary.pdf`):
  - **Capsid:** Mouse 2 is "AAV9-GFP @ 7.6e12 vg/kg" (summary PDF). **Assume AAV9 for every
    animal**; the per-mouse capsid doesn't matter for this work (decision: Andrew). Note the
    tension: up to ~28% GFP+ neurons near a target (pipeline B, `docs/gfp-cell-tagging.md`)
    is high for AAV9 at this dose in the literature, which makes the hand count (open 2)
    more important.
  - **CY5 = NeuN** on Mouse 2 (labelled on the summary PDF). Pipeline B counts neurons.
  - **No contrast MRI and no confocal.** Plan without either.
  - **Exposures differ** between slides and animals and we live with it: compare
    thresholded measures (A coverage, B per-slice fold), never raw intensity across slices.
  - **Head position:** nothing assumed per animal beyond the lab's Mouse 2 note (right side
    higher).

- **Sonication EDA (10/15)** (`scripts/sonication_eda.py` → `results/acoustic_eda/README.md`).
  This scripts most of "Findings that exist only here" below. Spreadsheet, mice 2–24: the
  design (6 burst counts × 6 goals) explains R² 0.77 of log dose; 22% of the residual dose
  and 41% of the residual drive voltage is animal-level; firing order = position for
  every target. Raw recordings, mice 1–6: the controller settles at ~0.93× its goal at
  every goal level; the interlock fired in 7/37 recordings; 2/37 hit the 0.125 V ceiling.

- **All 24 mice' recordings on disk (10/19).** `data/acoustic/US_Data.zip` (48.6 GB, from
  Dropbox) is read in place by `fus.archive` and is the **only copy**: the old unpacked
  `20260611/` folder was CRC-identical to it and was deleted, and every script (and
  `python -m fus Mouse_Cntr_01_Target1`) reads from the zip. Regenerated doses and features
  are identical. 27 GB free afterwards.
  Every one of the 145 recordings was reduced into the burst cache, and **all 138 targets for
  mice 2–24 match the sheet's cumulative 2nd harmonic to 4 decimals**. One file is misnamed:
  Mouse 17's `Target5Repeat` is position 6 (goal, bursts and dose match the sheet's P6; the
  lab's Reduced_ file says Target6). The controller findings hold in every session: it
  settles at 0.90–0.93× the goal; 31–39% of recordings average in band; the interlock fired in
  36/145 (25%). New: **voltage-ceiling hits rise across sessions, 2 → 1 → 4 → 8 of ~36**, so in
  session 4, 22% of targets may not have reached their goal at 0.125 V.
  (`results/acoustic_eda/README.md`, per-session table.)
- **Target template and mask (10/20)** (`docs/lab/gABI_template_100um.nii`,
  `Mask_AAVTargets.nii`; analysis in `docs/lab/README.md`). These are the lab's 9/18 files: the
  Allen template plus six 2 mm target spheres, matching the `posx.xyz` plan exactly. New
  information: **T1/T4 sit in anterior striatum, T2/T5 in lateral striatum/cortex, and T3/T6 in
  thalamus (VPM)**. Focal centres are 2.9–3.7 mm below the dorsal surface, while our area
  matching puts Mouse 2's sections at 1.7–2.7 mm, i.e. above most focal centres. Our assumed
  pattern offset from the brain centre (1.6 plan mm) agrees with the mask within ~0.3 mm.
- **Aligned section images (10/15; script and images removed 2026-11-04, in git history).** Every clicked section is re-rendered from its ds4 export as
  a dorsal view: front up, animal's left on the left (Mouse 1 flipped, since its notch is on
  the right), with ROIs labelled by target and coverage. There is one PNG per section and a
  grid per mouse. Checks: T1–T3 land right and T1/T4 front in both mice, at matching positions;
  green brightness inside each ROI ranks with its coverage (ρ 0.94 M2, 0.75 M1); outline
  IoU median 0.73 (M2) and 0.67 (M1).

## Open — blocking or slowing

1. **Re-tag the review verdicts** in the ingest notebook (step 6), so bad slices are excluded
   by your judgement rather than the post hoc low-signal cutoff. Then re-run
   `scripts/slices_dose_delivery.py`; `scripts/mouse02_dose_vs_coverage.py` also needs ✓s.
2. **Hand count for the cell pilot.** A blinded neuron and GFP+ count in a plume core, a
   plume edge and background. This is decision 001's pre-registered test.
3. **Cross-animal normalisation**, discussed 10/28.
   Decide the metric and threshold rule **before** measuring the four new brains.
   **Still open on sectioning:** which Mouse 2 sections pair up, the cutting order, the
   starting depth, and Mouse 1's protocol. Asked 10/20 and 11/01.
4. **Tissue throughput.** 2 of 24 brains imaged; 4 more due after 10/26.
5. ~~Remaining acoustic sessions~~ **Done 10/19**: all four sessions are in `data/acoustic/US_Data.zip`.
6. **FITC-only signal at Mouse 1 T6.** Unexplained; `notebooks/view_slice.ipynb` can now show
   it at full resolution (`section_s4`).

## Findings that exist only here

Computed from `data/Mouse_Controller_Data.xlsx` (all 24 mice) but **not yet scripted or
written into `docs/`**. Re-derive before citing:

- Every mouse has the same reference target, 120 bursts at goal 0.75. Across animals the
  cumulative 2nd harmonic ranges 0.70–1.73 (21% CV) and mean voltage 0.067–0.116 V (14% CV).
- After removing the N × goal trend, ~21% of the remaining variance is animal-level and ~79%
  target-level. That argues for within-animal contrasts.
- `Mean MPa` = 4.05 × `Mean Voltage` and `Mean MI` = MPa / √0.837. These are fixed
  conversions, so don't use them as independent predictors.
- The dose→emission slope by region runs from 88 (region 2) to 183 (region 6).
- Stored burst count = programmed + 15 for mice 1–6 and + 12 for mice 7–24.
- **(10/25) Cross-animal delivery**, from `results/histology/slices/cumulative_dose_delivery.csv`,
  12 targets. Pooled r with dose for ΣGFP+ cells ÷ ΣNeuN cells is +0.77 (+0.69 without the
  control; ρ +0.84). That beats dividing by each animal's mean target (+0.68). Raw totals
  give +0.09 (cells) and +0.31 (area), because Mouse 1 has 4 slices and Mouse 2 has 12, and
  NeuN per ROI per slice is ~1,600 (M1) vs ~5,100 (M2).
- **(10/25) Ceiling hits by position** (`results/acoustic_eda/recordings.csv`, mice 2–24): P1–P6
  had 0, 6, 1, 0, 5 and 3 hits of 23. Lateral P2/P5 account for 11/46 against 4/92 elsewhere
  (Fisher p 0.001). Lateral hits per session were 1, 1, 3 and 6.
- **(10/25) Candidate adjacent Mouse 2 pairs** (NeuN-count profile over the six ROIs, r ≥ 0.98):
  slide01_s2↔slide04_s2, slide01_s3↔slide04_s3, slide02_s3↔slide02_s4. Weak evidence; ask
  Josh, and don't set the order from it.
- DAPI vs NeuN in the 500 µm crop: no global shift (~1 px), but many individual cells don't
  line up. That fits each colour being in focus at a different depth.

## Traps for whoever picks this up

- **The spreadsheet's Mouse 1 rows 3–6 are rotated by one run.** Read doses from the recordings.
- **Slot ≠ target.** Mirroring comes from the notch clicks (`orientation.is_mirrored`), or from
  `KNOWN_MIRRORED` for Mouse 1.
- **`mouse01_slots.csv` is the GFP-fit reference** every Mouse 1 number in `docs/` comes from.
  Click-based runs write `mouseNN_template_slots.csv`; don't overwrite the reference.
- **Don't run `histology.tissue_mask` on a crop.** It separates tissue from glass by Otsu, and on
  an all-tissue crop it splits bright from dim tissue instead. Cut the section-level mask.
- **The segmentation models have their own environments** in `segmentation_alpha/`:
  `.venv-stardist-gpu` (TF 2.18 + tensorflow-metal; what `fus.cells` uses), `.venv-stardist`
  (CPU) and `.venv` (Cellpose). The project `.venv` has none of them. tensorflow-metal fails
  with TF 2.21, hence the pin. StarDist's slow step is its single-core outline merging, not
  the network, which is why `fus.cells` runs one process per ROI.
- **Pipeline B scores are fold over the slice's background.** Don't compare them with the
  pilot's earlier 0.0575 percentile scale; that scale tagged 5–14% of Mouse 1's no-FUS control.
- **Full-resolution exports are 7 GB each**, in `data/processed/histology/Mouse_NN/full/`
  (none there now). Delete them when done. Pipeline B doesn't need them: it
  exports only the ROI boxes and deletes them afterwards.
- **`docs/figures/histology-sections-overview.png` still shows 0.75 mm circles**, and no
  script makes it.
- **Use `.venv/bin/python`**: system `python3` is 3.9.
- **Agents here are text-only**: any regenerated figure needs a human to look at it.

## Next steps, in the order I would do them

**Before the 11/18 meeting** (added 11/11):

- Download the new Mouse 5 scans and the updated Mouse 3–5 decks from Dropbox.
- Ingest Mouse 6 (12 sections, one `.vsi` per slide) and the new Mouse 5 scans in
  `notebooks/ingest_new_histology.ipynb`.
- Match `IF_Images_Table.xlsx` rows to our scans: why the counts differ, and the QC reason for
  each exclusion.
- Rough cuts of the lab's three suggested comparisons: ~~region~~ and ~~T1 vs T4~~ done 11/13
  (`results/histology/region/`); same dose by different parameters still to do.

Older list:

1. Look over the pipeline B figures and the 10/6 ones (`cumulative_pooled.png`,
   `results/acoustic_features/model_comparison.png`). No human has looked at them yet.
2. **Call 10/28** (held). Settle the cross-animal metric and
   thresholds, then script them into `scripts/slices_dose_delivery.py`, and ingest the four
   new brains under those rules. Re-tag review verdicts first.
3. Hand-count the cell pilot (core, edge, background) and score StarDist plus the GFP+ calls.
   That settles decision 001.
4. Look at Mouse 1 T4 in the benchmark (the one failing check) and decide whether 0.10 is
   the right tolerance.
5. Script the spreadsheet EDA into `docs/controller-data-eda.md`.
6. When new tissue arrives (imaging 10/26): `notebooks/ingest_new_histology.ipynb`.
7. Later: port the feature extraction to MATLAB for the lab.

## Deliverable and timing

Draft journal manuscript by end of reading period (week 13). On current data that report is a
validated pipeline plus a consistent within-animal dose–delivery pattern in two animals, not a
trained predictive model. If only three or four animals are imaged by December, say so in
the design section rather than presenting it as a shortfall.
