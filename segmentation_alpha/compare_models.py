"""StarDist vs Cellpose-SAM on the 500 um crop: counts, sizes, run time, agreement.

    .venv/bin/python segmentation_alpha/compare_models.py

Reads the label images written by run_stardist.py and run_cellpose.py and writes
results/model_comparison_500um.csv. Agreement is model vs model (no hand count yet):
a cell counts as found by both when the outlines overlap at IoU >= 0.5.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

HERE = Path(__file__).resolve().parent
PX_UM = 0.325
RUNS = {"StarDist": "stardist_neun", "Cellpose-SAM": "cellpose_neun"}


def matched(a, b, iou_min=0.5):
    na, nb = int(a.max()) + 1, int(b.max()) + 1
    inter = np.bincount(a.ravel().astype(np.int64) * nb + b.ravel(), minlength=na * nb).reshape(na, nb)
    iou = inter / (inter.sum(1)[:, None] + inter.sum(0)[None, :] - inter + 1e-9)
    iou[0, :] = iou[:, 0] = 0
    return int((iou >= iou_min).any(axis=1).sum())


labels = {k: tifffile.imread(HERE / f"results/{v}_labels.tif") for k, v in RUNS.items()}
both = matched(*labels.values())
rows = []
for name, lab in labels.items():
    area = np.bincount(lab.ravel())[1:] * PX_UM**2
    d = 2 * np.sqrt(area[area > 0] / np.pi)
    info = json.loads((HERE / f"results/{RUNS[name]}.json").read_text())
    rows.append({"model": name, "cells": int(lab.max()), "median_diameter_um": round(float(np.median(d)), 1),
                 "diameter_p10_um": round(float(np.percentile(d, 10)), 1),
                 "diameter_p90_um": round(float(np.percentile(d, 90)), 1),
                 "seconds": info["seconds"], "found_by_both": both, "only_this_model": int(lab.max()) - both})
total = sum(r["cells"] for r in rows)
for r in rows:
    r["agreement_f1"] = round(2 * both / total, 3)
out = pd.DataFrame(rows)
out.to_csv(HERE / "results/model_comparison_500um.csv", index=False)
print(out.to_string(index=False))
