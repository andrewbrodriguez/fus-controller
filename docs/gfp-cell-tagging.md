# Cell-level GFP tagging — pilot

**Status:** first pass, 2026-09-30. One crop of one section (Mouse 2 `slide04_s3`), and
**not yet checked against a hand count**. Code and outputs are in
[`segmentation_alpha/`](../segmentation_alpha/). This is the cell-level alternative (option B)
in [decision 001](decisions/001-delivery-endpoint.md). The project's delivery endpoint is still
area coverage.

The question: instead of the fraction of an ROI's *pixels* that are GFP+, how many
*neurons* took up the virus?

<p align="center">
  <img src="figures/gfp-tagging-stains.png" alt="400 micrometre square of brain tissue at the edge of the T3 target: NeuN in grey and the anti-GFP stain in green, bright on the right where the plume is and dim on the left." width="49%">
  <img src="figures/gfp-tagging-cells.png" alt="The same square with the stains hidden: each segmented neuron outlined, yellow where tagged GFP-positive (mostly on the right) and magenta where negative (mostly on the left)." width="49%">
</p>

*Left: NeuN (grey) and the anti-GFP stain (green) in a 400 µm window on the edge of the T3
ROI. Right: the same window with the stains off. Each StarDist cell is outlined **yellow**
(GFP+) or **magenta** (GFP−) at the first-pass threshold: 269 cells, 77 GFP+. Both images
come from `segmentation_alpha/export_figures.py`, which reproduces the napari view.*

---

## Pipeline

```
full-resolution export of the section (0.325 µm/px)       notebooks/view_slice.ipynb
        │ make_crop.py            3 × 3 mm square centred on the finetuned T3 ROI
        ▼
crop (DAPI, native GFP, anti-GFP, NeuN)
        │ run_stardist.py         StarDist 2D_versatile_fluo on NeuN (CY5), defaults
        ▼
one outline per neuron
        │ gfp_tagging.py          mean normalised anti-GFP per cell → threshold
        ▼
GFP+ / GFP− per cell, by distance to the target           gfp_tagging.ipynb (napari)
```

| Step | Choice | Why |
|---|---|---|
| Crop | 3 × 3 mm around T3 (T3 reads 86% pixel coverage) | Covers the whole plume and reaches tissue > 1.5 mm from every target, giving a GFP-free reference |
| Cells | StarDist on NeuN | ~10× faster than Cellpose-SAM; see the comparison below |
| Normalisation | anti-GFP rescaled over the **whole section's** tissue: 0 = 1st percentile (70 counts), 1 = 99.9th (15,671). Linear, not clipped | The same scale whatever crop is analysed |
| Score | mean normalised anti-GFP over the cell's pixels | Your call (2026-09-30): an intensity-weighted average, one threshold |
| Threshold | Otsu on log₁₀ of the cell scores: **0.058** | Scores span two decades and are bimodal on a log scale: a background peak near 0.014 and a plume hump near 0.2–0.4. Log-Otsu sits in the valley |

## Results (3 × 3 mm crop)

StarDist found **16,723** NeuN cells in 84 s. At the log-Otsu threshold, **4,721 (28.2%)**
are GFP+, concentrated on the target:

| Distance to nearest target | Cells | GFP+ (log-Otsu, 0.058) | GFP+ (reference, 0.021) | GFP+ (linear Otsu, 0.261) |
|---|---|---|---|---|
| 0–0.25 mm | 308 | **100%** | 100% | 75% |
| 0.25–0.5 mm | 1,123 | **97%** | 99% | 64% |
| 0.5–0.87 mm (ROI edge) | 3,948 | **60%** | 88% | 24% |
| 0.87–1.25 mm | 6,595 | **14%** | 50% | 2% |
| 1.25–1.5 mm | 3,234 | **0.6%** | 16% | 0% |
| 1.5–2.25 mm | 1,515 | **0.2%** | ~1% | 0% |

The two alternative thresholds bracket the first pass:
- **Reference (0.021)** is the 99th percentile of random cell-shaped patches in tissue more
  than 1.5 mm from every target. It's permissive: it also takes in the plume's dim halo.
- **Linear Otsu (0.261)** is strict: it cuts into the plume population itself.

`THRESHOLD` in the notebook switches between them, or takes any number.

Source files, in `segmentation_alpha/results/`:
`mouse02_slide04_s3_T3_3000um_gfp_threshold.json`, `_gfp_cells.csv` and `_gfp_profile.csv`.

## Choosing the segmentation model

Both models ran with default settings on the NeuN channel of a 500 µm crop at the same
centre (`results/model_comparison_500um.csv`):

| | StarDist | Cellpose-SAM |
|---|---|---|
| Cells | 428 | 308 |
| Median diameter | 15.5 µm | 14.8 µm |
| 10th–90th percentile diameter | 9.2–18.4 µm | 10.5–19.2 µm |
| Run time | 2.3 s | 21.7 s (Apple GPU) |

They agree on 271 cells (overlap IoU ≥ 0.5; agreement F1 0.74). StarDist finds 157 that
Cellpose doesn't, and Cellpose 37 that StarDist doesn't. StarDist's smaller 10th
percentile suggests some of its extras are fragments or dim out-of-focus blobs.
**StarDist was chosen for speed.** This compares the models with each other; neither has
been scored against a hand count.

Each model runs in its own environment inside `segmentation_alpha/`: `.venv` holds
PyTorch and Cellpose, and `.venv-stardist` holds TensorFlow and StarDist. Both are
git-ignored.

## Caveats

- **No ground truth yet.** Neither the outlines nor the GFP+ calls have been checked against
  a blinded hand count. Decision 001's pre-registered test (GFP+ F1 ≥ 0.85, equal in plume
  core and edge) has not been run.
- **The threshold was chosen after seeing the histogram.** The plan was the reference
  threshold; the switch to log-Otsu came from the bimodal log histogram. All three stay
  visible in the notebook.
- **Bright surroundings can carry a cell over the threshold.** A neuron sitting in
  GFP-filled processes scores high even if its own body isn't transduced, and NeuN outlines
  include some cytoplasm. An earlier version subtracted a ring around each cell to remove
  the surroundings; only about a third of core neurons stood out from their surroundings.
  It was replaced by the plain mean, as asked.
- **One focal plane through a thick section.** DAPI and NeuN line up within about 1 px
  overall, but many individual cells don't match, which fits each colour being in focus at
  a different depth. This was measured once in the 500 µm crop and isn't scripted.
  Segmentation therefore uses NeuN alone.
- **What CY5 shows on Mouse 2 is unconfirmed.** Our docs call it NeuN on the basis of
  Mouse 1, but the lab's protocol uses CY5 for NeuN, GFAP or S100β depending on the panel.
  Ask Nick before reading these as neurons.
- **One crop, one section, one animal.** Nothing here has been run across targets or mice.

## Reproduce

```bash
.venv/bin/python segmentation_alpha/make_crop.py 3000
segmentation_alpha/.venv-stardist/bin/python segmentation_alpha/run_stardist.py \
    segmentation_alpha/mouse02_slide04_s3_T3_3000um.ome.tif
.venv/bin/python segmentation_alpha/gfp_tagging.py segmentation_alpha/mouse02_slide04_s3_T3_3000um.ome.tif
.venv/bin/python segmentation_alpha/export_figures.py        # the two images above
.venv/bin/python segmentation_alpha/compare_models.py        # needs run_cellpose.py output
```

`segmentation_alpha/gfp_tagging.ipynb` does the first three steps when anything is missing,
and opens the result in napari. Setting up the two environments:
`python3.12 -m venv segmentation_alpha/.venv && segmentation_alpha/.venv/bin/pip install cellpose packaging tifffile`,
and the same with `.venv-stardist` and `stardist tensorflow tifffile`.
