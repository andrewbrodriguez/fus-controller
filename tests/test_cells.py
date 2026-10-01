"""Pipeline B pieces that need neither QuPath nor StarDist: boxes and tagging."""

import numpy as np
import pytest
import tifffile

from fus import cells, rois


def test_box_covers_roi_plus_margin():
    roi = rois.Roi(1, 1, (1000.0, 2000.0), 100.0, 100.0)  # ds4 px
    x, y, w, h = cells.box(roi, full_um=0.325)
    margin = cells.MARGIN_UM / 0.325
    assert x == pytest.approx(2000 * 4 - 400 - margin, abs=1)
    assert y == pytest.approx(1000 * 4 - 400 - margin, abs=1)
    assert w == h == pytest.approx(2 * (400 + margin), abs=1)


def test_box_clipped_at_image_origin():
    x, y, _, _ = cells.box(rois.Roi(1, 1, (10.0, 10.0), 100.0, 100.0), full_um=0.325)
    assert x == 0 and y == 0


def test_tag_counts_centroids_inside_roi_and_thresholds(tmp_path):
    image = np.zeros((4, 200, 200), np.uint16)
    labels = np.zeros((200, 200), np.uint16)
    labels[95:105, 95:105] = 1    # centre: bright
    labels[60:70, 95:105] = 2     # inside: dim
    labels[5:15, 5:15] = 3        # outside the ROI
    image[2][labels == 1] = 1000
    image[2][labels == 2] = 120
    image[2][labels == 3] = 1000
    tifffile.imwrite(tmp_path / "crop.ome.tif", image, ome=True, metadata={"axes": "CYX"})
    tifffile.imwrite(tmp_path / "labels.tif", labels)
    roi = rois.Roi(1, 1, (100.0, 100.0), 50.0, 50.0)
    out = cells.tag(tmp_path / "crop.ome.tif", tmp_path / "labels.tif", roi, lo=100.0, hi=1100.0)
    assert sorted(out.label) == [1, 2]
    got = out.set_index("label")
    assert got.gfp_mean[1] == pytest.approx(0.9) and got.gfp_positive[1]
    assert got.gfp_mean[2] == pytest.approx(0.02) and not got.gfp_positive[2]


def test_source_maps_section_names():
    assert cells.source(1, "section_s4")[1] == 4
    assert cells.source(1, "section_s4")[0].name == "Image.vsi"
