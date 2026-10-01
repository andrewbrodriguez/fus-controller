Analysis notebooks. Keep exploratory work here; promote anything reusable into src/fus/.

| Notebook | What it shows |
|---|---|
| [`ingest_new_histology.ipynb`](ingest_new_histology.ipynb) | **Run this for each new animal.** Export → click centre/front/notch side in napari → measure → review placements in napari → coverage per target |
| [`view_slice.ipynb`](view_slice.ipynb) | One section at full resolution (0.325 µm/px) in napari: each stain as a layer, the saved T1–T6 ROIs on top. Exports on first open (~7 min, ~7 GB) |
| [`01_mouse01_target_regions.ipynb`](01_mouse01_target_regions.ipynb) | Opens the Mouse 1 sections and draws the six FUS targets on them, labelled with exposure and GFP coverage |
