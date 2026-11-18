# Where the project is — updated 2026-11-17 (evening)

Week ~10 of 13. This file is the moment-in-time picture; the durable method write-ups are in
`docs/`.

---

## One-paragraph summary

Both halves of the pipeline work, and **all six imaged brains are measured**. The acoustic
side reproduces the lab's own numbers on all 138 targets of mice 2–24. The tissue side turns a
whole-slide scan into delivery per target, two ways in the same hand-placed ROIs: **pixel
coverage (pipeline A)** and the **fraction of NeuN cells tagged GFP+ (pipeline B**, StarDist,
GFP+ above 4.05× that slice's background). The two agree at r = 0.97 over 204 ROIs. Mice 1–6
are oriented by hand, finetuned, measured and joined to dose: **34 sections, 36 targets**
(Mouse 3's `w10sl1_02_s0` was dropped as a bad slice on 11/17). Within every animal,
higher-dose targets carry more delivery (ΣGFP+÷ΣNeuN vs dose r = +0.63 to +0.86 in all six).
Pooled across animals that falls to r = +0.52 (ρ +0.65), down from +0.66 at five mice, because
Mouse 6's two highest doses do not deliver proportionally more. **No single mouse drives the
pooled result.** The one influential point is Mouse 6 T6, the highest dose in the study, and
only under a straight-line fit (leave-one-out, Since 11/17 evening). A pre-registered,
leave-one-mouse-out model comparison (`burst_ml/`) found that **per-burst acoustic detail adds
nothing to the cumulative 2nd harmonic**: dose predicts delivery in an unseen mouse, but the
mouse sets the level. Imaging is complete at the lab, so **six animals is probably the final
n** unless mice 7–24 get stained. The remaining gaps are validation, not throughput: **the
11/17 re-finetune of the zones isn't measured yet**, no hand count has scored pipeline B, Mice
1–2 have no review ✓ on record, and the cross-animal measure is still unsettled.

**Deliverable plan: analysis finishes 2026-11-30; from 2026-12-01 the work is the journal-style
write-up.**

## Since 11/17 (evening)

- **Mouse 3 `w10sl1_02_s0` excluded** as a bad slice (Andrew's call; `exclude = yes` with a
  note in `data/section_orientation.csv`). It was low-signal and mostly added NeuN cells. Mouse
  3 now has 4 kept sections, 34 in all. Re-ran `slices_dose_delivery.py`, `reference_target.py`,
  `region_delivery.py` and `dose_delivery_by_region.py`. What moved:
  - pooled ΣGFP+÷ΣNeuN r +0.52 (unchanged), ρ +0.63 → +0.65; GFP+ area r +0.56, ρ +0.67 → +0.66;
  - Mouse 3 within-mouse r +0.85 → +0.86 (cells), +0.82 → +0.84 (pixels);
  - T1 < T4 still 5/5 mice, but cells p 0.047 → **0.052**, pixels 0.083 → 0.071. T3 > T6 p 0.067
    → 0.074. Thalamus at matched dose +0.07 (p 0.85).
  - **`burst_ml/` was not re-run.** It reads `cumulative_dose_delivery.csv`, so its committed
    results are for the 35-section data. Expect small shifts.
- **Finetuning moved out of the ingest notebook** into `scripts/finetune_rois.py`: every kept
  section of every mouse in one napari window (`--mouse N` for one). It prints what moved.
- **Zones re-finetuned (11/17) but not yet measured.** Zones moved in 25 of 34 kept sections, by
  up to 0.5 mm, most in Mice 1 and 2; sizes are unchanged. Mouse 5's `w11sl2_01_s2` was
  finetuned for the first time. **Every number in this file still uses the old zones** (Open 0).
- **Sections ranked for quality, all six mice** (`scripts/rank_sections.py` →
  `data/section_quality.csv`, rank 1 = best, no notes yet). The window shows each section
  straightened with its zones, so it needs orientation clicks: Mouse 6's five excluded sections
  have none and couldn't be ranked. **The ranking wasn't blind**: the thumbnails show the GFP.
- **Top-k sections** (`scripts/top_k_sections.py` → `results/histology/slices/top_k_sections.*`):
  dose vs delivery using only each mouse's k best-ranked sections. Pooled r over 36 targets:

  | k | Sections | A (pixels) | B (ΣGFP+÷ΣNeuN) |
  |---|---|---|---|
  | 1 | 6 | +0.68 | +0.55 |
  | 2 | 12 | +0.67 | +0.55 |
  | 3 | 18 | +0.62 | +0.46 |
  | 6 | 27 | +0.55 | +0.46 |
  | 12 (all) | 34 | +0.56 | +0.52 |

  - Fewer, better sections read a little higher on pixels, but the 90% interval on r is about
    ±0.15–0.2, so it's within noise.
  - Mice disagree. Mouse 2 improves with more sections (B +0.57 → +0.82). Mouse 6 gets worse
    (A +0.84 → +0.60). Mouse 4's pixel r drops from +0.91 to +0.45 when its 3rd and 4th
    sections go in.
  - Before using top-k as a rule: the ranking saw the GFP (circular), and k would be picked
    after seeing this trace.
- **Leave-one-out** (`scripts/leave_one_out.py` → `results/histology/slices/leave_one_out.*`).
  - **No mouse is a clear outlier.** Pixels: dropping any one mouse moves r by at most 0.05
    (0.53–0.61 vs 0.56), inside what dropping 6 random targets does (−0.07 to +0.10).
  - Cells: only Mouse 6 stands out, r 0.52 → 0.69 without it (5% of random 6-target drops move
    r that far). Its slope is about 0.2× the others' and it reads 21 points below their fit,
    but it still ranks its own targets (r 0.63).
  - Mice differ in level on cells (M1 +17, M2 +13, M3 −9, M6 −21 points against the others'
    fit). That's the mouse effect, not an outlier.
  - **One dose dominates: Mouse 6 T6** (3.33, the highest in the study; 17% GFP+). Cook's
    distance 1.2 (A) and 2.0 (B), against a flag of 0.11. Dropping it alone: A 0.56 → 0.64, B
    0.52 → 0.67. No other target moves r by more than 0.07.
  - **It's a straight-line effect.** ρ barely moves without it (A 0.66 → 0.67, B 0.65 → 0.67),
    and on log dose B goes 0.57 → 0.65. That is A15 again, which the pre-registered log-dose
    model already handles. Quote ρ or the log-dose fit rather than linear r.
- **Missing scans traced (lab sheet `docs/lab/IF_Images_Table.xlsx` vs the files).** The sheet
  lists 9 more imaged sections for Mice 3–5 than we have (10/9/7 vs 7/5/5), and all 9 are the
  first scan (`_01`) on their slide; the slide overviews arrived without them. Mapping rule
  (inferred, not confirmed): sheet "section n" = scan `_0n`, or the nth section inside a
  multi-section file. Mouse 3's `3.10.2_01` image data is present but misfiled inside
  `_Image_Cntr_3.11.2_03_/` with no `.vsi` (header read: 16-bit, ~34,800 × 29,200 px);
  `3.11.2_01` has a `.vsi` but no data. Re-downloads of Mice 3 and 4 were byte-identical to what
  we had (zips deleted after checksum checks). Up to 8 usable sections if the lab re-uploads.
  Ask on 11/18.
- **The sheet's "ignore" verdicts are all applied.** It marks 4 sections "Not ok / Not great /
  Not very good": 3 are excluded here, the 4th (Mouse 4 `4.11.2_01`) never arrived. We also
  exclude 8 the sheet doesn't flag (our calls): Mouse 3 `w10sl1_02_s0`, `w10sl1_03_s0`; Mouse 5
  `w11sl2_01_s3`; Mouse 6's five. Mouse 4's two "Pretty good – thalamus lost" sections are kept;
  they are its two lowest-ranked, the ones that drop its pixel r in the top-k trace.
- **Within-mouse trends and the mouse effect** (`scripts/mouse_effect.py` →
  `results/histology/mouse_effect/`; cells, logit vs log dose, 35 sonicated targets).
  - Positive in all six mice. Repeated-measures r +0.67 on cells (p 6e-5) and +0.60 on pixels.
    Separate slopes fit no better than one shared slope (p 0.40): doubling the dose multiplies
    the odds of a neuron being GFP+ by 2.8 (90%: 1.9–4.0).
  - Mice differ a lot in level: at dose 1, M1/M2 36%, M5 22%, M4 17%, M6 10%, M3 7% GFP+.
    That's an odds ratio of 7 between highest and lowest, about 4× in dose. Mouse level is 34%
    of the variance on cells.
  - **Exploratory, six mice:** level doesn't track the drive voltage on the shared treatment (r
    −0.08) and only loosely its dose (r +0.54, p 0.27).
- **Not yet looked at by a person:** `top_k_sections.png`, `leave_one_out.png`,
  `mouse_effect/within_mouse_trend.png`, `mouse_effect/mouse_level_vs_acoustics.png`, and the
  regenerated slices and region figures.

## Since 11/16

- **Burst-level models (`burst_ml/`, merged 11/17).** Question: does modelling each burst
  predict delivery better than the cumulative 2nd harmonic, for a mouse the model hasn't seen?
  35 sonicated targets in Mice 1–6, delivery ΣGFP+÷ΣNeuN. Seven beta-binomial models (numpyro),
  each scored by fitting on five mice and predicting the sixth. The models, priors and the
  rule for winning were committed (d90a650) before any model saw delivery data. Full write-up
  in `burst_ml/README.md`; generated tables in `burst_ml/results/README.md`.
  - **Verdict: no.** Against the AUC model (logit delivery linear in log dose, one offset per
    mouse), Δelpd is M4 timing +0.27 (SE 1.25), M5a per-burst power −0.27 (0.59), M5b power
    + bands −0.54 (0.91). The bar was Δ > 0 and Δ/SE > 2. M4's small gain comes from holding
    out Mouse 1 alone.
  - **The lab's third suggestion** (same dose as few intense bursts vs many weak): per-burst
    power ρ = 1.28 (90%: 0.80–2.04). ρ = 1 means only the total counts, and it is well inside
    the interval. Report it as "no detectable difference at n = 35", not "no difference".
    "Intensity" here is each burst's emitted 2nd harmonic, not drive pressure.
  - **Timing:** γ = −0.78 (−1.74 to +0.09), leaning toward early bursts counting more; 0 is
    inside the interval. **Broadband, sub- and ultraharmonic weights are ≈ 0.**
  - **Dose itself predicts delivery in an unseen mouse:** +7.11 over no dose (SE 1.38, 5.1
    SE). The held-out predicted order within a mouse matches the observed order (mean Spearman
    +0.71). The 50% and 90% held-out intervals cover 49% and 94% of targets.
  - **A15 resolved:** a log-dose curve covers Mouse 6's flattening. A plateau (Emax) predicts no
    better (−0.25, SE 0.42).
  - **Dose for 20% GFP+ neurons in a new mouse: 0.86 (90%: 0.22–3.24).** That's a 15-fold
    range, mostly mouse-to-mouse (SD 0.61 on the logit scale). The mouse, not the burst record,
    is now the limit for a delivery-setpoint controller. An animal-level covariate (drive
    voltage needed, reference-target dose) is the next question, but six mice give it six data
    points.
  - **Robustness (not in the plan):** without Mouse 1 T1 (three sonications), every verdict
    holds.
  - **Four figures in `burst_ml/results/` have not been looked at by a person.**
- **Review ✓ recorded for every kept section of Mice 5 and 6** (11/17).

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

### Results at six mice (34 sections since 11/17 evening)

| Pooled over 36 targets | r | ρ |
|---|---|---|
| ΣGFP+ ÷ ΣNeuN | **+0.52** (was +0.66 at five mice) | +0.65 |
| GFP+ area ÷ animal mean | +0.56 | +0.66 |
| GFP+ cells ÷ animal mean | +0.47 | +0.63 |

Per animal, ΣGFP+÷ΣNeuN vs dose: M1 +0.75, M2 +0.82, M3 +0.86, M4 +0.64, M5 +0.75, M6 +0.63.
Within-slice ρ is positive in 21/21 clear slices (A, median +0.60) and 20/21 (B, +0.71).

**Why pooling got weaker:** without Mouse 6 the pooled r is +0.69, and most of that is its T6 alone (leave-one-out). Mouse 6 ranks its own
targets fine, but its two big doses (T3 2.02, T6 3.33) read 24.6% and 17.1% — below a straight
line. That is **assumption A15**: delivery is being fitted as linear in dose, and it looks like
it saturates. **Done since (Since 11/16):** a log-dose fit scored leave-one-mouse-out handles
it, and a plateau adds nothing. That fit is the "baseline predictive model" the proposal
promised.

**Mouse 6 as an out-of-sample test** of what the lab suggested on 11/09. Predictions were
written in chat before measuring but **were never committed**, so this is weaker than a real
pre-registration — say so if it goes in the write-up.
- ✓ dose ranks delivery within sections (3/3 clear slices, 5/7 overall)
- ✓ T5 lowest (dose 0.43 → 1.2% GFP+)
- ✗ T6 highest: T3 is (dose 2.02 → 24.6%) with T6 second (17.1%)
- ✓ T1 < T4 on cells (9.6% vs 15.6%) despite T1 getting *more* dose (1.04 vs 0.92); a tie on pixels

**Region, at six mice** (`results/histology/region/README.md`):
- **No region effect at matched dose.** Thalamus +0.07, p 0.85 (it has the most raw delivery
  and the most dose).
- **T1 < T4 holds in 5/5 mice**, both measures, but weaker than at five mice: cells p 0.052
  (0.047 before Mouse 3's bad slice was dropped; 0.018 at five mice), pixels p 0.071 (0.034 at
  five mice). Mouse 5's new section narrowed its gap.
- **New**: T3 (right thalamus) > T6 (left) in 5/5 mice, p 0.074. So it is not "the right side
  delivers less" — T1 vs T4 is specific.
- Still confounded: T1 is always fired first, so side and firing order can't be separated.

## Since 11/11

- **Lab update 11/09** (correspondence kept locally in `docs/lab/correspondence.md`). **Meeting
  Wed 11/18, 11:15.** All imaging done except Mouse 5 well 11 slide 1. AAV: mice 1–12 AAV9,
  13–24 AAV.CPP16, so **every imaged mouse is AAV9** and capsid is not a factor. Three analyses
  suggested: by region (done), T1 vs T4 (done), and same dose from different parameters (done
  11/17 in `burst_ml/`; see Since 11/16).
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
  permutation p 0.23 and failed leave-one-mouse-out). Treat the screen as hypotheses. The
  15-feature screen itself hasn't been re-run at six mice, but `burst_ml/` (Since 11/16) has
  now tested bands, timing and per-burst power at 35 targets: none beat the AUC.

## What exists and is trusted

| Piece | State |
|---|---|
| `src/fus/extract.py` | Python port of `nt_ExtractHarmonicData.m`. All 138 targets of mice 2–24 match the lab's sheet to 4 dp (`docs/acoustic-extraction.md`) |
| `src/fus/histology.py` | `.vsi` → per-target GFP coverage. ROI r = 1.0 template mm (1.74 mm diameter on tissue). `export --prefix/--series` |
| `src/fus/orientation.py` | napari: `annotate` (centre, front, notch side) and `review` (✓/✗). Each window writes only its own fields |
| `src/fus/rois.py` + `data/roi_locations.csv` | Hand-finetuned T1–T6 shapes per section. **All kept sections of Mice 1–6 finetuned**, re-finetuned 11/17 (not yet measured) |
| `data/section_orientation.csv` | Clicks for all 36 sections of Mice 1–6. Review ✓ on every kept section of Mice 3–6; **none on record for Mice 1–2** (see Open 1) |
| `src/fus/cells.py` | **Pipeline B**: per-ROI cell count and GFP+ fraction. Run on every kept section of Mice 1–6; ~1–1.5 min per section on `segmentation_alpha/.venv-stardist-gpu` |
| `src/fus/archive.py` | Reads recordings straight out of `data/acoustic/US_Data.zip` (48.6 GB, all 24 mice). **The only copy — don't unpack it** |
| `src/fus/reconstruct.py` | Loose 3D stack of one animal: clicks for alignment, atlas spacing, Gaussian depth blending |
| `notebooks/ingest_new_histology.ipynb` | Per-animal workflow: export → orient → measure → review |
| `scripts/finetune_rois.py` | napari: move and resize T1–T6 on every kept section of every mouse in one window (`--mouse N` for one); writes `data/roi_locations.csv` |
| `scripts/rank_sections.py` + `data/section_quality.csv` | Qt window: rank one mouse's sections by quality from straightened thumbnails with zones. All six mice ranked 11/17. Only `top_k_sections.py` reads it |
| `scripts/top_k_sections.py`, `scripts/leave_one_out.py` | Sensitivity of the pooled result: each mouse's top-k ranked sections; drop each mouse, then each target |
| `scripts/mouse_effect.py` | Within-mouse trends, each mouse's level (shared-slope fit), and levels vs the shared treatment's acoustics → `results/histology/mouse_effect/` |
| `tests/` | 37 tests, passing |
| `burst_ml/` | Pre-registered burst-level models. Own environment `burst_ml/.venv` (jax, numpyro). Per-burst cache in `burst_ml/cache/` (git-ignored; rebuild with `.venv/bin/python -m bml.cache` from `burst_ml/`, reads the zip). Full run `.venv/bin/python -m bml.run` (~20 min); 55 tests |

## Open — blocking or slowing

0. **Re-measure with the 11/17 zones** (25 of 34 sections changed).
   - Pipeline A: the notebook's Measure cell, once per mouse. `fus.histology measure` ignores
     finetuned shapes (see Traps).
   - Pipeline B: `fus.cells --mouse N --section …` on the changed sections, ~25–40 min.
   - Then `slices_dose_delivery` → `reference_target` → `region_delivery` →
     `dose_delivery_by_region` → `top_k_sections` → `leave_one_out` → `burst_ml` (`bml.run`,
     ~20 min), and update every quoted number.
1. **Review verdicts.** Every kept section of Mice 3–6 is reviewed ✓ (Mice 5–6 recorded
   11/17; all placements correct). Mice 1 and 2 have no verdict on record (4 and 12 kept
   sections). The analysis only drops sections marked *wrong*, so they are included either way.
2. **Hand count for pipeline B.** A blinded neuron and GFP+ count in a plume core, a plume edge
   and background. This is decision 001's pre-registered test and the biggest hole a reviewer
   would find. **Nothing else validates the delivery measure.**
3. **Cross-animal normalisation is still unsettled.** ΣGFP+÷ΣNeuN is the candidate (needs no
   rescale), but it has not been chosen formally. Decide it before the write-up, not during.
   Leave-one-out adds a second choice: linear r swings on one point (Mouse 6 T6); ρ and the
   log-dose fit don't.
4. **Why Mouse 6's five sections were excluded** isn't recorded anywhere. They have no clicks,
   so the ranking window can't show them, and `section_quality.csv` has no notes yet.
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
- **The design does contrast few-intense vs many-weak bursts** (11/17, from the `burst_ml`
  cache; 35 targets):
  - bursts per target 55–450, median 100 (450 is Mouse 1 T1's three sonications);
  - mean 2nd harmonic per burst 42–289, median 76;
  - 37 pairs of targets within 15% of each other's dose but at least 1.5× apart in burst count.

  So the ρ null isn't from a lack of variation.

## Traps for whoever picks this up

- **`python -m fus.histology measure` ignores finetuned shapes.** It measures plain template
  circles, and the `AGENTS.md` command writes to `mouseNN_template_slots.csv`, the file that
  holds the finetuned results. Measure through the notebook's Measure cell until the CLI reads
  `data/roi_locations.csv`.
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
- Look at the figures before showing them: the four `burst_ml/results/*.png`, and
  `results/histology/slices/top_k_sections.png` and `leave_one_out.png`.
- Meeting questions:
  - the missing Mouse 5 well 10 sections, and why the QC sheet lists more sections than we have;
  - section order, spacing and depth from Josh;
  - whether any of Mice 7–12 (AAV9, like 1–6) will be stained, since the mouse effect is now
    the limit;
  - per-animal records that might explain the mouse effect (weight, microbubble lot, injection
    timing);
  - whether intensity in the third suggestion meant drive pressure or emitted harmonic;
  - which cross-animal measure to use, and whether to quote ρ / the log-dose fit rather than
    linear r;
  - whether choosing sections by a quality ranking is acceptable (if so: ranked blind to the
    GFP, with k fixed in advance);
  - Andrew's idea: total delivery over each slice, not only inside the ROIs;
  - co-authorship and target journal.

**Analysis to finish by 11/30:**
1. **Re-measure with the 11/17 zones** (open 0) and re-run everything downstream.
2. **Hand count** (open 2) — the only validation of pipeline B.
3. Decide the cross-animal measure, and whether the headline is ρ / the log-dose fit rather than
   linear r. Write it into the scripts as the default.
4. Decide the section rule: every kept section, or top-k by a quality ranking done blind to the
   GFP with k fixed first.
5. Freeze the dataset with a git tag and re-run every script from scratch, `burst_ml` included.
6. Optional: re-run the 15-feature screen at 35 targets. `burst_ml` already covers bands,
   timing and per-burst power, so expect nothing new.

**From 12/01: write-up.** Draft journal manuscript by end of reading period. On current data
that is:
- a validated pipeline;
- a consistent within-animal dose–delivery relationship in six animals;
- a calibrated leave-one-mouse-out baseline (log dose with a per-mouse offset);
- a pre-registered null for burst-level features.

It is not a controller that hits a delivery target in a new mouse. Say that in the design
section rather than presenting it as a shortfall.
