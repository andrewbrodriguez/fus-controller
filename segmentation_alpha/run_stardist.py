"""StarDist (2D_versatile_fluo) on the NeuN channel of the pilot crop.

    segmentation_alpha/.venv-stardist/bin/python segmentation_alpha/run_stardist.py [crop.ome.tif]

Runs in its own environment (TensorFlow). With no argument it segments the 500 um
crop into results/stardist_neun_labels.tif; given a crop, into
results/<crop name>_stardist_labels.tif. One integer per cell, 0 = background,
plus a .json beside it.
"""

import json
import sys
import time
from pathlib import Path

import tifffile
from csbdeep.utils import normalize
from stardist.models import StarDist2D

HERE = Path(__file__).resolve().parent
CROP = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE / "mouse02_slide04_s3_T3_500um.ome.tif"
STEM = "stardist_neun" if len(sys.argv) == 1 else CROP.name.split(".")[0] + "_stardist"
CHANNEL = 3  # CY5 = NeuN, order DAPI, FITC, TRITC, CY5

image = tifffile.imread(CROP)[CHANNEL]
model = StarDist2D.from_pretrained("2D_versatile_fluo")
t = time.time()
n_tiles = tuple(max(1, s // 1024) for s in image.shape)  # keep memory bounded on big crops
labels, _ = model.predict_instances(normalize(image, 1, 99.8), n_tiles=n_tiles)  # default thresholds
seconds = time.time() - t

out = HERE / "results"
out.mkdir(exist_ok=True)
tifffile.imwrite(out / f"{STEM}_labels.tif", labels.astype("uint32"), compression="zlib")
info = {"model": "StarDist 2D_versatile_fluo (defaults)", "channel": "CY5 (NeuN)",
        "crop": CROP.name, "n_tiles": list(n_tiles),
        "normalize_percentiles": [1, 99.8], "prob_thresh": float(model.thresholds.prob),
        "nms_thresh": float(model.thresholds.nms), "n_cells": int(labels.max()),
        "seconds": round(seconds, 1)}
(out / f"{STEM}.json").write_text(json.dumps(info, indent=2) + "\n")
print(info)
