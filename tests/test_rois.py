"""Finetuned ROI geometry: napari corners <-> ellipse, mask area, file round trip."""

import numpy as np
import pytest

from fus import orientation, rois


def test_corners_round_trip_rotated_ellipse():
    roi = rois.Roi(4, 1, (300.0, 400.0), 80.0, 50.0, 30.0)
    back = rois.Roi.from_corners(4, 1, roi.corners(scale=4.0), scale=4.0)
    assert back.centre == pytest.approx(roi.centre)
    assert (back.semi_a, back.semi_b) == pytest.approx((80.0, 50.0))
    assert back.angle_deg % 180 == pytest.approx(30.0)


def test_mask_area_is_pi_ab():
    roi = rois.Roi(1, 1, (200.0, 200.0), 60.0, 30.0, 45.0)
    _, m = roi.mask((400, 400))
    assert m.sum() == pytest.approx(np.pi * 60 * 30, rel=0.02)


def test_mask_clipped_at_image_edge():
    roi = rois.Roi(1, 1, (5.0, 5.0), 20.0, 20.0)
    (rs, cs), m = roi.mask((100, 100))
    assert rs.start == 0 and cs.start == 0 and m.shape == (rs.stop, cs.stop)


def test_template_rois_follow_mirroring():
    a = orientation.Annotation(2, "s", (500.0, 500.0), (500.0, 100.0), (700.0, 500.0))
    plain = rois.from_template(a, 1.3, mirrored=False)
    flipped = rois.from_template(a, 1.3, mirrored=True)
    assert plain[1].centre == pytest.approx(flipped[4].centre)
    assert plain[1].semi_a == pytest.approx(1000 * 0.87 / 1.3)


def test_locations_file_round_trip(tmp_path):
    loc = {(2, "slide01_s3"): {1: rois.Roi(1, 1, (10.0, 20.0), 5.0, 4.0, 10.0)}}
    path = tmp_path / "loc.csv"
    rois.write_locations(loc, path)
    back = rois.read_locations(path)[(2, "slide01_s3")][1]
    assert back.centre == (10.0, 20.0) and back.semi_b == 4.0 and back.angle_deg == 10.0
