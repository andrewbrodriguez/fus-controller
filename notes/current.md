# Where the project is — updated 2026-11-16

Week ~10 of 13. This file is the moment-in-time picture; the durable method write-ups are in
`docs/`.

---

## One-paragraph summary

Both halves of the pipeline work, and **all six imaged brains are measured**. The acoustic
side reproduces the lab's own numbers on all 138 targets of mice 2–24. The tissue side turns a
whole-slide scan into delivery per target, two ways in the same hand-placed ROIs: **pixel
coverage (pipeline A)** and the **fraction of NeuN cells tagged GFP+ (pipeline B**, StarDist,
GFP+ above 4.05× that slice's background). The two agree at r = 0.97 over 204 ROIs. Mice 1–6
are oriented by hand, finetuned, measured and joined to dose: **35 sections, 36 targets**.
Within every animal, higher-dose targets carry more delivery (ΣGFP+÷ΣNeuN vs dose r = +0.63 to
+0.85 in all six). Pooled across animals that falls to r = +0.52, down from +0.66 at five mice,
because Mouse 6's two highest doses do not deliver proportionally more. Imaging is complete at
the lab, so **six animals is probably the final n** unless mice 7–24 get stained. The
remaining gaps are validation, not throughput: no hand count has scored pipeline B, the new
sections have no review ✓, and the cross-animal measure is still unsettled.

**Deliverable plan: analysis finishes 2026-11-30; from 2026-12-01 the work is the journal-style
write-up.**

## Since 11/13

- **Mouse 5 re-downloaded and Mouse 6 ingested (11/16).** The old 3.2 GB `Mouse_05/` was
  deleted and replaced with the lab's full re-upload (5.4 GB, 23 files, CRCs verified), which
  adds **well 11 slide 2** (`Image_Cntr 5.11.2_01.vsi`, two sections in one file). The updated
  Mouse 3/4/5 decks (ordered roughly dorsal to ventral) replaced the 11/01 ones in
  `data/histology/`. Andrew clicked and finetuned the new sections; both pipelines were re-run.
  - **Mouse 5 now has 3 kept sections** (was 2): `w10sl1_02_s0`, `w10sl2_02_s0`,
    `w11sl2_01_s2`. `w11sl2_01_s3` is excluded.
  - **Mouse 6 has 7 kept sections of 12**: `w10sl1_01_s3`, `_s4`, `w11sl1_01_s2`, `_s4`,
    `_s5`, `w11sl2_01_s2`, `_s3`. Five are excluded, but the lab's QC sheet flags only one
    ("only one half of brain"). **Andrew's reasons for the other four aren't written down** —
    worth noting before the meeting.
  - **Well 10 is still short.** The sheet lists 5 imaged sections for Mouse 5 well 10; we have
    3 scans. The re-upload didn't add them. Ask on 11/18.
- **Bug fixed: template-placed sections were silently dropped.** `measure` leaves the `target`
  column blank when it uses plain template circles instead of finetuned shapes, and
  `scripts/slices_dose_delivery.py` dropped those sections. Mouse 5's new `w11sl2_01_s2` was
  the first to hit it. The script now maps slot → target from the notch click, like the ingest
  notebook. **Any new reader of `mouseNN_template_slots.csv` must do the same**
  (`docs/histology-pipeline.md`, Known problems 0).
- **All analyses re-run on six mice** (`MICE` in `slices_dose_delivery.py` is now the single
  place the mouse list lives; `region_delivery.py` and `reference_target.py` follow it).
- **New figure**: `scripts/dose_delivery_by_region.py` →
  `results/histology/region/dose_delivery_by_region.png`. Every target on one axis, dose
  against %GFP+ neurons, coloured by region. Human-checked 11/16. It shows the regions mixed
  throughout: thalamus sits further right (more dose), not higher at matched dose.

### Results at six mice

| Pooled over 36 targets | r | ρ |
|---|---|---|
| ΣGFP+ ÷ ΣNeuN | **+0.52** (was +0.66 at five mice) | +0.63 |
| GFP+ area ÷ animal mean | +0.56 | +0.67 |
| GFP+ cells ÷ animal mean | +0.47 | +0.63 |

Per animal, ΣGFP+÷ΣNeuN vs dose: M1 +0.75, M2 +0.82, M3 +0.85, M4 +0.64, M5 +0.75, M6 +0.63.
Within-slice ρ is positive in 21/21 clear slices (A, median +0.60) and 20/21 (B, +0.71).

**Why pooling got weaker:** without Mouse 6 the pooled r is +0.68. Mouse 6 ranks its own
targets fine, but its two big doses (T3 2.02, T6 3.33) read 24.6% and 17.1% — below a straight
line. That is **assumption A15**: delivery is being fitted as linear in dose, and it looks like
it saturates. A fit against log dose, scored by leaving out one mouse at a time, is the next
model to try and is the "baseline predictive model" the proposal promised.

**Mouse 6 as an out-of-sample test** of what the lab suggested on 11/09. Predictions were
written in chat before measuring but **were never committed**, so this is weaker than a real
pre-registration — say so if it goes in the write-up.
- ✓ dose ranks delivery within sections (3/3 clear slices, 5/7 overall)
- ✓ T5 lowest (dose 0.43 → 1.2% GFP+)
- ✗ T6 highest: T3 is (dose 2.02 → 24.6%) with T6 second (17.1%)
- ✓ T1 < T4 on cells (9.6% vs 15.6%) despite T1 getting *more* dose (1.04 vs 0.92); a tie on pixels

**Region, at six mice** (`results/histology/region/README.md`):
- **No region effect at matched dose.** Thalamus +0.05, p 0.89 (it has the most raw delivery
  and the most dose).
- **T1 < T4 holds in 5/5 mice**, both measures, but weaker than at five mice: cells p 0.047
  (was 0.018), pixels p 0.083 (was 0.034). Mouse 5's new section narrowed its gap.
- **New**: T3 (right thalamus) > T6 (left) in 5/5 mice, p 0.067. So it is not "the right side
  delivers less" — T1 vs T4 is specific.
- Still confounded: T1 is always fired first, so side and firing order can't be separated.

## Since 11/11

- **Lab update 11/09** (correspondence kept locally in `docs/lab/correspondence.md`). **Meeting
  Wed 11/18, 11:15.** All imaging done except Mouse 5 well 11 slide 1. Virus: mice 1–12 AAV9,
  13–24 AAV.CPP16, so **every imaged mouse is AAV9** and capsid is not a factor. Three analyses
  suggested: by region (done), T1 vs T4 (done), and **same dose from different parameters —
  still to do**.
- **New lab files, in `docs/lab/`**: `AAV_By_TreatmentLocation.pptx` (bursts, MPa, dose per
  target, Mice 1–6; Mouse 6's doses match our recordings) and `IF_Images_Table.xlsx`
  (per-section imaging QC, Mice 2–6). **The sheet lists more sections than we have scans** for
  Mice 3/4/5 (10/9/7 vs 7/5/3). Mapping its rows to our scan names is not done.

## Since 11/04

- **Region analysis, first pass (11/13)** (`scripts/region_delivery.py`). Superseded by the
  six-mouse run above; method unchanged: logit(delivery) against log dose with one offset per
  mouse, residuals compared by region and by left/right pair, 20,000 within-mouse permutations.

## Since 11/01

- **Repo simplified (11/04): aggregate figures only.** Removed ~60 per-mouse, per-section and
  per-session images and `scripts/aligned_sections.py`. The numbers stay in the CSVs; old
  images are in git history.

## Since 10/28 (updated 2026-11-01)

- **Mice 3, 4, 5 imaged and measured.** Format from Mouse 3 on: one scan per `.vsi`, named
  `Image_Cntr_<mouse>.<well>.<slide>_<scan>` (export prefix `w<well>sl<slide>_<scan>`), with
  **three channels, CY5/TRITC/FITC, no DAPI**. Mouse 3's `_Image_Cntr_3.11.2_03_/` also holds a
  stray `_Image_Cntr 3.10.2_01_` stack and a `.vsi` that won't open.
- **Coverage varies a lot between sections, most likely because of cut depth.** Some sections
  show only T2, in Mouse 4 and Mouse 2 alike, suggesting the lateral plumes reach higher into
  cortex. A hypothesis, not tested.

## Since 10/20

- **Lab update 10/23**: sections 30 µm thick, ~360 µm apart (to confirm with Josh); Mouse 2 is
  probably two offset series of 6. Preferred outcome measures: GFP+ cells summed over sections,
  then area above a per-mouse threshold, then background-normalised mean intensity. Ceiling
  hits are across the four surgery days, not time within a session.

## Since 9/30

- **Meeting 10/2.** Delivery is summed over all slices of an animal, not normalised per slice.
- **Acoustic features vs delivery** (`src/fus/features.py`,
  `scripts/acoustic_features_vs_delivery.py`). 15 pre-registered features vs 2f AUC alone, on
  11 labelled targets. **No multi-feature model beats AUC convincingly** (best gain had
  permutation p 0.23 and failed leave-one-mouse-out). Treat the screen as hypotheses. **Not yet
  re-run at six mice** — worth doing now that there are 36 targets.

## What exists and is trusted

| Piece | State |
|---|---|
| `src/fus/extract.py` | Python port of `nt_ExtractHarmonicData.m`. All 138 targets of mice 2–24 match the lab's sheet to 4 dp (`docs/acoustic-extraction.md`) |
| `src/fus/histology.py` | `.vsi` → per-target GFP coverage. ROI r = 1.0 template mm (1.74 mm diameter on tissue). `export --prefix/--series` |
| `src/fus/orientation.py` | napari: `annotate` (centre, front, notch side) and `review` (✓/✗). Each window writes only its own fields |
| `src/fus/rois.py` + `data/roi_locations.csv` | Hand-finetuned T1–T6 shapes per section. **All kept sections of Mice 1–6 finetuned** |
| `data/section_orientation.csv` | Clicks for all 36 sections of Mice 1–6. **Review verdicts are patchy** (see Open 1) |
| `src/fus/cells.py` | **Pipeline B**: per-ROI cell count and GFP+ fraction. Run on every kept section of Mice 1–6; ~1–1.5 min per section on `segmentation_alpha/.venv-stardist-gpu` |
| `src/fus/archive.py` | Reads recordings straight out of `data/acoustic/US_Data.zip` (48.6 GB, all 24 mice). **The only copy — don't unpack it** |
| `src/fus/reconstruct.py` | Loose 3D stack of one animal: clicks for alignment, atlas spacing, Gaussian depth blending |
| `notebooks/ingest_new_histology.ipynb` | Per-animal workflow: export → orient → finetune → measure → review |
| `tests/` | 37 tests, passing |

## Open — blocking or slowing

1. **Review verdicts.** Mice 5 and 6's new sections have no ✓/✗. The analysis only drops
   sections marked *wrong*, so unreviewed ones are currently included. Re-tag in the ingest
   notebook (step 6), then re-run `scripts/slices_dose_delivery.py`.
2. **Hand count for pipeline B.** A blinded neuron and GFP+ count in a plume core, a plume edge
   and background. This is decision 001's pre-registered test and the biggest hole a reviewer
   would find. **Nothing else validates the delivery measure.**
3. **Cross-animal normalisation is still unsettled.** ΣGFP+÷ΣNeuN is the candidate (needs no
   rescale), but it has not been chosen formally. Decide it before the write-up, not during.
4. **Why Mouse 6's five sections were excluded** isn't recorded anywhere.
5. **Section order, spacing and starting depth** are still assumed (filename order, 250 µm).
   Asked 10/20 and 11/01; ask Josh directly.
6. **FITC-only signal at Mouse 1 T6.** Unexplained.

## Findings that exist only here

Computed but **not yet scripted or written into `docs/`**. Re-derive before citing:

- Every mouse has the same reference target, 120 bursts at goal 0.75. Across animals the
  cumulative 2nd harmonic ranges 0.70–1.73 (21% CV) and mean voltage 0.067–0.116 V (14% CV).
- After removing the N × goal trend, ~21% of the remaining variance is animal-level and ~79%
  target-level. That argues for within-animal contrasts.
- `Mean MPa` = 4.05 × `Mean Voltage` and `Mean MI` = MPa / √0.837 — fixed conversions, not
  independent predictors.
- Stored burst count = programmed + 15 for mice 1–6 and + 12 for mice 7–24.
- **Ceiling hits by position** (mice 2–24): P1–P6 had 0, 6, 1, 0, 5, 3 of 23. Lateral P2/P5
  account for 11/46 against 4/92 elsewhere (Fisher p 0.001). Per session: 1, 1, 3, 6.
- **Candidate adjacent Mouse 2 pairs** (NeuN-count profile, r ≥ 0.98): slide01_s2↔slide04_s2,
  slide01_s3↔slide04_s3, slide02_s3↔slide02_s4. Weak; don't set the order from it.

## Traps for whoever picks this up

- **Template-placed sections leave `target` blank** in `mouseNN_template_slots.csv`. Map slot →
  target from the notch click. This silently cost a slice once (fixed 11/16).
- **The spreadsheet's Mouse 1 rows 3–6 are rotated by one run.** Read doses from the recordings.
- **Slot ≠ target.** Mirroring comes from the notch clicks (`orientation.is_mirrored`), or
  `KNOWN_MIRRORED` for Mouse 1.
- **`mouse01_slots.csv` is the GFP-fit reference** every Mouse 1 number in `docs/` comes from.
  Click-based runs write `mouseNN_template_slots.csv`; don't overwrite the reference.
- **Don't run `histology.tissue_mask` on a crop.** On an all-tissue crop Otsu splits bright from
  dim tissue instead of tissue from glass.
- **The segmentation models have their own environments** in `segmentation_alpha/`:
  `.venv-stardist-gpu` (what `fus.cells` uses), `.venv-stardist` (CPU), `.venv` (Cellpose).
- **Pipeline B scores are fold over the slice's background**, not the pilot's 0.0575 percentile
  scale.
- **Disk is tight (~15 GB free).** Exports are ~0.4 GB per section; full-resolution exports are
  ~7 GB each. Delete them when done.
- **Use `.venv/bin/python`**: system `python3` is 3.9.
- **Agents here are text-only**: any regenerated figure needs a human to look at it.

## Next steps, in the order I would do them

**Before the 11/18 meeting:**
- Re-tag review verdicts for the new Mouse 5 and 6 sections (Andrew; napari).
- Rough cut of the lab's third suggestion: **same cumulative dose reached by different
  parameters** (few bursts at high intensity vs many at low). Expect few matched pairs — likely
  a descriptive table, not a test.
- Meeting questions: the missing Mouse 5 well 10 sections; why the QC sheet lists more sections
  than we have; section order/spacing/depth from Josh; whether any of mice 7–24 will be
  stained; which cross-animal measure to use; co-authorship and target journal.

**Analysis to finish by 11/30:**
1. **Saturating fit** (delivery vs log dose), scored leave-one-mouse-out. This is the baseline
   predictive model the proposal promised, and A15 says the linear fit is wrong.
2. **Hand count** (open 2) — the only validation of pipeline B.
3. Re-run the acoustic feature screen at six mice (36 targets, was 11).
4. Freeze the dataset with a git tag and re-run every script from scratch.
5. Decide the cross-animal measure and write it into the scripts as the default.

**From 12/01: write-up.** Draft journal manuscript by end of reading period. On current data
that is a validated pipeline plus a consistent within-animal dose–delivery relationship in six
animals — not a trained predictive model. Say that in the design section rather than presenting
it as a shortfall.
