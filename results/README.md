Generated figures and model outputs. Image files here are git-ignored; regenerate them from the pipeline.

| Path | What | Made by |
|---|---|---|
| `histology/mouse01_slots.csv` | Mouse 1, GFP-fit placement: every Mouse 1 number in `docs/` | `python -m fus.histology measure` |
| `histology/mouseNN_template_slots.csv` | **Pipeline A**: pixel coverage per section × target, in the hand-placed ROIs | `notebooks/ingest_new_histology.ipynb` |
| `histology/cells/mouseNN_cell_rois.csv` | **Pipeline B**: cells, density, and fraction GFP+ per section × target, same ROIs | `python -m fus.cells` |
| `histology/slices/` | Both pipelines against dose, per slice and pooled; its own README | `scripts/slices_dose_delivery.py` |
| `histology/mouse01_benchmark.csv` | Hand placement scored against Mouse 1's answer key | `scripts/benchmark_mouse01.py` |
| `histology/mouse01_dose_vs_coverage.csv` | Mouse 1 dose join (GFP fit) | `scripts/mouse01_dose_vs_coverage.py` |
