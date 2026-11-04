Generated CSVs and figures. Only **aggregate** figures (all mice pooled) are committed;
per-mouse and per-section numbers live in the CSVs. Regenerate with the script named below,
never by hand.

| Path | What | Made by |
|---|---|---|
| `histology/slices/` | **Headline.** Both pipelines against dose, cumulative per target and pooled across mice; its own README | `scripts/slices_dose_delivery.py`, `scripts/reference_target.py` |
| `histology/mouseNN_template_slots.csv` | **Pipeline A**: pixel coverage per section × target, in the hand-placed ROIs | `notebooks/ingest_new_histology.ipynb` |
| `histology/cells/mouseNN_cell_rois.csv` | **Pipeline B**: cells, density, and fraction GFP+ per section × target, same ROIs | `python -m fus.cells` |
| `histology/mouse01_slots.csv` | Mouse 1, GFP-fit placement: every Mouse 1 number in `docs/` | `python -m fus.histology measure` |
| `histology/mouse01_dose_vs_coverage.csv` | Mouse 1 dose join (GFP fit) | `scripts/mouse01_dose_vs_coverage.py` |
| `histology/mouse01_benchmark.csv` | Hand placement scored against Mouse 1's answer key | `scripts/benchmark_mouse01.py` |
| `histology/reconstruction/*.csv` | Mouse 2 section spacing from the atlas (the GIF it also makes is git-ignored) | `scripts/mouse02_reconstruction_3d.py` |
| `acoustic_eda/` | Design, controller behaviour and noise floor over all 24 mice; its own README | `scripts/sonication_eda.py` |
| `acoustic_features/` | Acoustic features vs delivery, model comparison; its own README | `scripts/acoustic_features_vs_delivery.py` |
