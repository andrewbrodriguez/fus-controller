# segmentation_alpha — cell segmentation pilot

The pilot work for counting GFP+ neurons (pipeline B), done on crops of one Mouse 2 slice
(`slide04_s3`, centred on target T3). The production version is
[`src/fus/cells.py`](../src/fus/cells.py); the method, numbers and caveats are in
[`docs/gfp-cell-tagging.md`](../docs/gfp-cell-tagging.md).

## Files

| File | What |
|---|---|
| `make_crop.py` | Cuts a square (500 µm, 1 mm or 3 mm) at full resolution around T3 from the full-section export |
| `run_stardist.py`, `run_cellpose.py` | Segment NeuN in a crop with each model, writing to `results/` |
| `compare_models.py` | StarDist vs Cellpose-SAM on the 500 µm crop → `results/model_comparison_500um.csv` |
| `gfp_tagging.py` | Scores every cell's anti-GFP and sets the thresholds, including the production 4.05× background → `results/*_gfp_*.csv/json` |
| `export_figures.py` | The tagged-cell figures in `docs/figures/gfp-tagging-*.png` |
| `view_crop.ipynb` | 500 µm crop in napari, with the StarDist vs Cellpose outlines |
| `gfp_tagging.ipynb` | Any crop in napari with every cell tagged, a threshold you can change, and distance-from-target plots |
| `*.json` | Where each crop came from (section, pixel bounds, pixel size) |
| `*.ome.tif` | The crops (git-ignored; remake with `make_crop.py`) |

## Environments

Segmentation models live in their own environments here, not in the project `.venv`:

| Environment | Holds | Used by |
|---|---|---|
| `.venv-stardist-gpu` | TensorFlow 2.18 + tensorflow-metal (Apple GPU) + StarDist | `fus.cells` (production), by default |
| `.venv-stardist` | TensorFlow 2.21 (CPU) + StarDist | `run_stardist.py`; the fallback for `fus.cells` |
| `.venv` | PyTorch + Cellpose 4 | `run_cellpose.py` |

All three are git-ignored. To rebuild them, with Python 3.12:

```bash
python3.12 -m venv segmentation_alpha/.venv-stardist-gpu
segmentation_alpha/.venv-stardist-gpu/bin/pip install "tensorflow==2.18.*" tensorflow-metal stardist tifffile
python3.12 -m venv segmentation_alpha/.venv-stardist
segmentation_alpha/.venv-stardist/bin/pip install stardist tensorflow tifffile
python3.12 -m venv segmentation_alpha/.venv
segmentation_alpha/.venv/bin/pip install "cellpose>=4" packaging tifffile
```

tensorflow-metal 1.2.0 fails to load with TensorFlow 2.21, so the GPU environment pins 2.18.
On one 5,539 px crop, GPU and CPU give the same 5,717 cells.

## Running the pilot again

```bash
.venv/bin/python segmentation_alpha/make_crop.py 3000
segmentation_alpha/.venv-stardist/bin/python segmentation_alpha/run_stardist.py \
    segmentation_alpha/mouse02_slide04_s3_T3_3000um.ome.tif
.venv/bin/python segmentation_alpha/gfp_tagging.py segmentation_alpha/mouse02_slide04_s3_T3_3000um.ome.tif
.venv/bin/python segmentation_alpha/export_figures.py
```

`make_crop.py` needs the full-resolution export of `slide04_s3`, which
`notebooks/view_slice.ipynb` makes (~7 min, ~7 GB). `gfp_tagging.ipynb` runs any missing step
itself.
