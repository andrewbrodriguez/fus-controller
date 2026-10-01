# AGENTS.md — FUS controller project

Mapping acoustic emissions recorded during focused-ultrasound BBB opening to how much AAV
actually reached the tissue. Background and study design: [README.md](README.md). Charter and
weekly plan: [ARCH.md](ARCH.md).

**Read [`notes/current.md`](notes/current.md) first** — it is the moment-in-time state: what
works, what is blocked, what to do next.

## Environment

Use `.venv/bin/python`. It is Python 3.12 with the project installed editable, so no
`PYTHONPATH` is needed. System `python3` is 3.9 and will not work.

```bash
.venv/bin/python -m pytest tests -q          # 33 tests: geometry, orientation, ROIs, cells
```

Jupyter kernel: **Python 3.12 (fus-research .venv)**.

## Commands

```bash
# acoustic: one recording, or a whole day to CSV
.venv/bin/python -m fus data/acoustic/20260611/Mouse_Cntr_01_Target1.mat
.venv/bin/python -m fus data/acoustic/20260611/*_Target*.mat --csv results/day1.csv --quiet

# tissue: export sections (needs QuPath 0.7), then measure
.venv/bin/python -m fus.histology export data/histology/Mouse_01/Image.vsi \
    data/processed/histology/Mouse_01/ds4
.venv/bin/python -m fus.histology measure data/processed/histology/Mouse_01/ds4/section_s*.ome.tif \
    --csv results/histology/mouse01_slots.csv --figures results/histology/figures

# Mouse 2+: use notebooks/ingest_new_histology.ipynb (export, napari clicks, measure,
# napari review). The pieces, for reference:
.venv/bin/python -m fus.histology export data/histology/Mouse_02/Control_02_Slide_01.vsi \
    data/processed/histology/Mouse_02/ds4 --prefix slide01
.venv/bin/python -m fus.orientation data/processed/histology/Mouse_02/ds4/*.ome.tif \
    data/processed/histology/Mouse_01/ds4/section_s*.ome.tif
.venv/bin/python -m fus.histology measure data/processed/histology/Mouse_02/ds4/*.ome.tif \
    --orientation data/section_orientation.csv --csv results/histology/mouse02_template_slots.csv

# join dose to delivery, and the README example figure
.venv/bin/python scripts/mouse01_dose_vs_coverage.py
.venv/bin/python scripts/mouse02_dose_vs_coverage.py   # stops until orientation is recorded
.venv/bin/python scripts/benchmark_mouse01.py          # hand placement vs Mouse 1 answer key
.venv/bin/python scripts/slices_dose_delivery.py       # per-slice + pooled dose vs coverage

# pipeline B: ROI -> StarDist cells -> GFP-tagged cells (needs segmentation_alpha/.venv-stardist-gpu)
.venv/bin/python -m fus.cells --mouse 2 --section slide04_s3     # ~1-1.5 min per section
.venv/bin/python -m fus.cells --mouse 2 --all                    # every clicked section

# cell-level pilot (own environments; see docs/gfp-cell-tagging.md)
.venv/bin/python segmentation_alpha/make_crop.py 3000
segmentation_alpha/.venv-stardist/bin/python segmentation_alpha/run_stardist.py segmentation_alpha/<crop>.ome.tif
.venv/bin/python segmentation_alpha/gfp_tagging.py segmentation_alpha/<crop>.ome.tif
.venv/bin/python scripts/mouse01_example_figure.py
```

## Layout

| Path | What |
|---|---|
| `src/fus/extract.py` | Acoustic feature extraction; port of the lab's MATLAB script |
| `src/fus/histology.py` | Whole-slide scan → GFP coverage per target |
| `scripts/` | One-off analyses and figure generators; `scripts/qupath/` holds Groovy |
| `docs/` | Method write-ups; `docs/decisions/` holds decision records |
| `notes/` | Meeting notes and the current-state file |
| `data/` | Only the summary sheet is tracked; raw data is git-ignored (see `data/README.md`) |
| `results/` | Generated CSVs and figures |
| `src/fus/orientation.py`, `src/fus/rois.py` | napari tools: orientation clicks, review, finetuned target shapes |
| `src/fus/cells.py` | Pipeline B: per-ROI cell counts and GFP+ fraction (StarDist via `scripts/stardist_segment.py`) |
| `notebooks/` | `ingest_new_histology.ipynb` (per-animal workflow), `cell_pipeline_one_slice.ipynb` (A vs B on one slice), `view_slice.ipynb` (full-resolution viewer) |
| `segmentation_alpha/` | Cell segmentation pilot and the model environments: `.venv-stardist-gpu` (used by `fus.cells`), `.venv-stardist` (CPU), `.venv` (Cellpose). See its README |

## Ground rules

- **Raw data is the source of truth, not the summary sheet.** `Mouse_Controller_Data.xlsx`
  has rotated rows for Mouse 1 (rows 3–6 carry other runs' values), and `Mean MPa` / `Mean MI`
  are fixed multiples of voltage rather than measurements.
- **Recording → target:** `TargetN` was fired at position N for every animal except Mouse 1.
  `data.niscope.tracktime` gives each recording's acquisition time, so firing order is always
  recoverable. See `data/README.md`.
- **Orientation comes from a person, not the fit.** From Mouse 2 on, sections are mounted at
  any rotation and the target pattern nearly matches itself turned 180°, so a person clicks each
  section's centre, front and notch side in `notebooks/ingest_new_histology.ipynb`
  (→ `data/section_orientation.csv`), and reviews the placements there too.
  The tool is a GUI; an agent can't run it for the user. Don't
  infer them from GFP — that would bake the dose–delivery answer into the measurement.
- **Delivery is measured on the anti-GFP stain (TRITC)**, not native GFP (FITC).
- **Two delivery measures, same ROIs.** Pipeline A (pixels, `fus.histology`) and pipeline B
  (GFP+ neurons, `fus.cells`). B scores each cell as fold over *its slice's* background
  (median anti-GFP over tissue), because backgrounds differ several-fold between animals. A
  fixed absolute threshold tagged 5–14% of Mouse 1's no-FUS control.
- **Figures are generated by scripts.** Regenerate; never hand-edit an image.
- **Numbers quoted in `docs/` come from committed CSVs.** If you change a parameter — the ROI
  radius, say — re-run and update every table that quotes a number, in the same pass.
- **Flag assumptions.** `docs/mouse01-dose-delivery.md` keeps an assumptions table with what
  each one rests on; add to it rather than quietly relying on something new.
- Agents here are **text-only**: no reading PDFs or images directly (use `pdftext`, or the
  `.pdf.md` sidecar), and any regenerated figure needs a human to look at it.
