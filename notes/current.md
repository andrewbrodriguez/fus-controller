# Where the project is — 2026-09-30 (evening)

Week ~5 of 13. This file is the moment-in-time picture; the durable method write-ups are in
`docs/`.

---

## One-paragraph summary

Both halves of the pipeline work. The acoustic side reproduces the lab's own numbers, and
the tissue side turns a whole-slide scan into GFP coverage per target. **Two brains are
measured.** All 16 usable sections of Mouse 1 and Mouse 2 were oriented by hand (centre,
front and notch side, plus finetuned target circles in napari) and joined to acoustic dose.
Within each slice, higher-dose targets tend to carry more GFP in both animals: 11 of 11
slices with clear signal show a positive rank correlation (median ρ +0.54). A **cell-level
pilot** (StarDist on NeuN, cells tagged GFP+ by mean normalised anti-GFP) works on one
3 mm crop. Neither has been checked against a hand count. Tissue remains the limit: 2 of
24 brains imaged.

## What exists and is trusted

| Piece | State |
|---|---|
| `src/fus/extract.py` | Python port of `nt_ExtractHarmonicData.m`. Run on 19 recordings from mice 1–3; all 15 that can be matched to a spreadsheet row unambiguously match it exactly (`docs/acoustic-extraction.md`) |
| `src/fus/histology.py` | `.vsi` → per-target GFP coverage. ROI r = 1.0 mm. Placement by GFP fit (Mouse 1 reference), clicked template, or finetuned shapes (`analyse_section(..., template=, rois=)`). `export --prefix/--series` |
| `src/fus/orientation.py` | napari: `annotate` (click centre, front, notch side; the view straightens to front-up, notch-left) and `review` (✓/✗). Each window writes only its own fields, so two open windows can't overwrite each other |
| `src/fus/rois.py` + `data/roi_locations.csv` | Hand-finetuned T1–T6 shapes per section, measured exactly as drawn. **All 16 sections finetuned** |
| `data/section_orientation.csv` | Clicks for all 17 sections (Mouse 1 s3 excluded). **No review verdicts**: the ones given were lost to the overwrite bug, since fixed |
| `notebooks/ingest_new_histology.ipynb` | Per-animal workflow: export → orient → finetune → measure → review → coverage per target |
| `notebooks/view_slice.ipynb` | Any section at full resolution (0.325 µm/px) in napari with ROIs. Exports on first open: ~7 min, ~7 GB per section |
| `scripts/slices_dose_delivery.py` | Dose vs coverage per slice and pooled; within-slice z-scores. Output and README in `results/histology/slices/` |
| `scripts/benchmark_mouse01.py` | Scores a placement method against Mouse 1's answer key (control, known front, validated fit) |
| `segmentation_alpha/` | Cell-level pilot: crops, StarDist/Cellpose environments, GFP tagging, napari notebooks. See `docs/gfp-cell-tagging.md` |
| `tests/` | 29 tests, passing |

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
higher (Nick, 9/29).

**Cell-level pilot** (`docs/gfp-cell-tagging.md`):
- **Model choice:** on 500 µm of NeuN, StarDist finds 428 cells in 2 s and Cellpose-SAM 308
  in 22 s; they agree on 271. StarDist was chosen for speed.
- **3 mm crop around Mouse 2 slide04_s3 T3:** 16,723 cells, 28% GFP+ at the log-Otsu
  threshold. That's 97–100% within 0.5 mm of the target centre and under 1% beyond 1.25 mm.

**Mouse 2 doses** (from the recordings; match the sheet): T1 0.62, T2 1.02, T3 0.77, T4 1.07,
T5 0.37, T6 1.60 (interlock fired).

## Settled with Nick

- **Recording → position (9/18).** `TargetN` = position N for every animal except Mouse 1 (first
  three sonications all at position 1, then 4, 5, 6, then a 7th at 2; position 3 never
  sonicated).
- **Notch (9/18, 9/29).** Right hemisphere on Mouse 1, left from Mouse 2 on. Our clicks agree:
  Mouse 1's notch clicks plus its empty target 3 put targets 1–3 on the animal's right,
  consistently across all four sections.
- **Delivery measure (9/18).** The anti-GFP stain (TRITC), not native GFP.
- **Focal spot (9/18).** ~2 × 2 mm in x/y, ~3 mm in z; precise dimensions to come.
- **Feature priority (9/18).** Cumulative 2nd-harmonic AUC first.

## Open — blocking or slowing

1. **Re-tag the review verdicts** in the ingest notebook (step 6), so bad slices are excluded
   by your judgement rather than the post hoc low-signal cutoff. Then re-run
   `scripts/slices_dose_delivery.py`; `scripts/mouse02_dose_vs_coverage.py` also needs ✓s.
2. **Hand count for the cell pilot.** A blinded neuron and GFP+ count in a plume core, a
   plume edge and background. This is decision 001's pre-registered test.
3. **Ask Nick:** the capsid per mouse (AAV9 vs AAV.CPP16), **which marker CY5 is on Mouse 2**
   (NeuN, GFAP or S100β), section order and spacing, imaging exposure settings, and a resend
   of the Mouse 2 summary PDF and the template/mask attachments.
4. **Tissue throughput.** 2 of 24 brains.
5. **Remaining acoustic sessions.** 3 of 4 days are still on Dropbox.
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
- **The segmentation pilot has its own environments** (`segmentation_alpha/.venv` for Cellpose,
  `.venv-stardist` for StarDist). The project `.venv` has neither.
- **Full-resolution exports are 7 GB each**, in `data/processed/histology/Mouse_NN/full/`.
  Delete them when done.
- **`data/histology/Mouse_02.zip` (17.8 GB)** is still there after extraction.
- **`docs/figures/histology-sections-overview.png` still shows 0.75 mm circles**, and no
  script makes it.
- **Use `.venv/bin/python`**: system `python3` is 3.9.
- **Agents here are text-only**: any regenerated figure needs a human to look at it.

## Next steps, in the order I would do them

1. Re-tag review verdicts, re-run the slice analysis, and send Nick an update with the
   two-animal result and the questions above.
2. Hand-count the cell pilot (core, edge, background) and score StarDist plus the GFP+ calls.
   That settles decision 001.
3. Look at Mouse 1 T4 in the benchmark (the one failing check) and decide whether 0.10 is
   the right tolerance.
4. Script the spreadsheet EDA into `docs/controller-data-eda.md`.
5. When new tissue arrives: `notebooks/ingest_new_histology.ipynb`.
6. Later: port the feature extraction to MATLAB for the lab.

## Deliverable and timing

Draft journal manuscript by end of reading period (week 13). On current data that report is a
validated pipeline plus a consistent within-animal dose–delivery pattern in two animals, not a
trained predictive model. If only three or four animals are imaged by December, say so in
the design section rather than presenting it as a shortfall.
