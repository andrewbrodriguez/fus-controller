"""Segment NeuN cells with StarDist, one label image per input crop.

Runs in the StarDist environment, not the project's (TensorFlow lives only there):

    segmentation_alpha/.venv-stardist/bin/python scripts/stardist_segment.py \\
        --channel 3 crop1.ome.tif labels1.tif [crop2.ome.tif labels2.tif ...]

`fus.cells` calls this; the model loads once for all crops. Settings match the
pilot (docs/gfp-cell-tagging.md): 2D_versatile_fluo, default thresholds, input
normalised to its 1st-99.8th percentile.
"""

import argparse
import json
import time

import numpy as np
import tifffile
from csbdeep.utils import normalize
from stardist.models import StarDist2D

p = argparse.ArgumentParser()
p.add_argument("--channel", type=int, default=3, help="NeuN channel index (CY5 = 3)")
p.add_argument("pairs", nargs="+", help="crop, labels, crop, labels, ...")
args = p.parse_args()
if len(args.pairs) % 2:
    p.error("give input/output pairs")

model = StarDist2D.from_pretrained("2D_versatile_fluo")
for src, dst in zip(args.pairs[::2], args.pairs[1::2]):
    with tifffile.TiffFile(src) as tf:
        image = tf.series[0].levels[0].asarray()[args.channel]
    t = time.time()
    n_tiles = tuple(max(1, s // 1024) for s in image.shape)
    labels, _ = model.predict_instances(normalize(image, 1, 99.8), n_tiles=n_tiles, show_tile_progress=False)
    dtype = np.uint16 if labels.max() < 65535 else np.uint32
    tifffile.imwrite(dst, labels.astype(dtype), compression="zlib")
    print(json.dumps({"crop": src, "labels": dst, "cells": int(labels.max()),
                      "seconds": round(time.time() - t, 1)}), flush=True)
