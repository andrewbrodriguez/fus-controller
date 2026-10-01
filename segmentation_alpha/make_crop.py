"""Cut the segmentation pilot crop: 500 x 500 um at the centre of T3, Mouse 2 slide04_s3.

    .venv/bin/python segmentation_alpha/make_crop.py            # 500 um
    .venv/bin/python segmentation_alpha/make_crop.py 1000       # 1 x 1 mm, same centre

Reads the full-resolution export (0.325 um/px; make it with notebooks/view_slice.ipynb
if missing) and the finetuned T3 ROI from data/roi_locations.csv. Writes, next to
this script:

  mouse02_slide04_s3_T3_<side>um.ome.tif   all four channels, full resolution,
                                           uint16, lossless
  mouse02_slide04_s3_T3_<side>um.json      where the crop came from
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import tifffile
import zarr

from fus import histology, rois

ROOT = histology.REPO_ROOT
HERE = Path(__file__).resolve().parent
MOUSE, SECTION, TARGET = 2, "slide04_s3", 3
FULL = ROOT / f"data/processed/histology/Mouse_{MOUSE:02d}/full/{SECTION}.ome.tif"
DS4 = ROOT / f"data/processed/histology/Mouse_{MOUSE:02d}/ds4/{SECTION}.ome.tif"


def main(side_um: float = 500.0) -> None:
    name = f"mouse{MOUSE:02d}_{SECTION}_T{TARGET}_{side_um:.0f}um"
    roi = rois.read_locations()[(MOUSE, SECTION)][TARGET]
    full_um, ds4_um = rois.export_pixel_um(FULL), rois.export_pixel_um(DS4)
    r, c = (v * ds4_um / full_um for v in roi.centre)
    half = round(side_um / full_um / 2)
    r0, c0 = round(r) - half, round(c) - half
    r1, c1 = r0 + 2 * half, c0 + 2 * half

    with tifffile.TiffFile(FULL) as tf:
        meta = tf.ome_metadata
    channels = [ch.split('Name="', 1)[1].split('"', 1)[0] for ch in meta.split("<Channel ")[1:]]
    level0 = zarr.open(tifffile.imread(FULL, aszarr=True), mode="r")["0"]
    crop = level0[:, r0:r1, c0:c1]

    out = HERE / f"{name}.ome.tif"
    tifffile.imwrite(out, crop, ome=True, compression="zlib", metadata={
        "axes": "CYX", "Channel": {"Name": channels},
        "PhysicalSizeX": full_um, "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": full_um, "PhysicalSizeYUnit": "µm",
    })
    info = {
        "mouse": MOUSE, "section": SECTION, "target": TARGET, "side_um": side_um,
        "pixel_um": full_um, "shape_cyx": list(crop.shape), "channels": channels,
        "source": str(FULL.relative_to(ROOT)),
        "rows": [r0, r1], "cols": [c0, c1],
        "roi_centre_fullres_px": [round(r, 1), round(c, 1)],
        "roi": "finetuned T3 from data/roi_locations.csv, crop centred on it",
    }
    (HERE / f"{name}.json").write_text(json.dumps(info, indent=2) + "\n")
    print(f"wrote {out.relative_to(ROOT)}  {crop.shape} {crop.dtype}, "
          f"{out.stat().st_size / 1e6:.1f} MB, rows {r0}:{r1} cols {c0}:{c1}")


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 500.0)
