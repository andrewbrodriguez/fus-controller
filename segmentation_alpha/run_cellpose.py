"""Cellpose-SAM on the NeuN channel of the pilot crop.

    segmentation_alpha/.venv/bin/python segmentation_alpha/run_cellpose.py

Runs in its own environment (PyTorch). Writes results/cellpose_neun_labels.tif
(one integer per cell, 0 = background) and results/cellpose_neun.json.
"""

import json
import time
from pathlib import Path

import tifffile
from cellpose import models

HERE = Path(__file__).resolve().parent
CROP = HERE / "mouse02_slide04_s3_T3_500um.ome.tif"
CHANNEL = 3  # CY5 = NeuN, order DAPI, FITC, TRITC, CY5

image = tifffile.imread(CROP)[CHANNEL]
model = models.CellposeModel(gpu=True)  # default model: cpsam (Cellpose-SAM)
t = time.time()
labels, _, _ = model.eval(image, flow_threshold=0.4, cellprob_threshold=0.0)
seconds = time.time() - t

out = HERE / "results"
out.mkdir(exist_ok=True)
tifffile.imwrite(out / "cellpose_neun_labels.tif", labels.astype("uint32"), compression="zlib")
info = {"model": "Cellpose-SAM (cpsam, defaults)", "channel": "CY5 (NeuN)",
        "flow_threshold": 0.4, "cellprob_threshold": 0.0, "n_cells": int(labels.max()),
        "seconds": round(seconds, 1), "device": str(model.device)}
(out / "cellpose_neun.json").write_text(json.dumps(info, indent=2) + "\n")
print(info)
