"""Geometry checks for fus.histology on synthetic sections (no image files needed)."""

import numpy as np
import pytest
from scipy import ndimage as ndi

from fus import histology as h

PIXEL_UM = 10.4  # 32x-downsampled scan


def synthetic_density(angle, scale, centre, missing=(), shape=(900, 900)):
    """GFP density map with a blob at each planned target except ``missing``."""
    img = np.zeros(shape, np.float32)
    pts = h._plan_points(h.PLAN_MM, *centre, angle, scale, PIXEL_UM)
    for k, (r, c) in pts.items():
        if k not in missing:
            img[int(round(r)), int(round(c))] = 1.0
    return ndi.gaussian_filter(img, 250 / PIXEL_UM), pts


@pytest.mark.parametrize("angle, scale", [(90, 0.82), (75, 0.95), (110, 1.1)])
def test_full_fit_recovers_transform(angle, scale):
    density, truth = synthetic_density(angle, scale, (450, 430))
    fit = h.fit_plan(density, PIXEL_UM)
    assert abs(((fit.angle_deg - angle + 180) % 360) - 180) < 2
    assert fit.scale == pytest.approx(scale, abs=0.03)
    for k, (r, c) in fit.points().items():
        assert np.hypot(r - truth[k][0], c - truth[k][1]) < 5  # ~50 um


def test_leave_one_out_places_empty_target():
    density, truth = synthetic_density(90, 0.85, (460, 440), missing={3})
    full = h.fit_plan(density, PIXEL_UM)
    loo = h.fit_plan(density, PIXEL_UM, leave_out=3, start=full)
    r, c = loo.points()[3]
    assert np.hypot(r - truth[3][0], c - truth[3][1]) < 5


def test_recentre_round_trip():
    used = {k: v for k, v in h.PLAN_MM.items() if k != 2}
    r, c = h._recentre(100.0, 200.0, 63.0, 0.9, PIXEL_UM, h.PLAN_MM, used)
    assert (r, c) != (100.0, 200.0)
    back = h._recentre(r, c, 63.0, 0.9, PIXEL_UM, used, h.PLAN_MM)
    assert back == pytest.approx((100.0, 200.0))


def test_rescaled_fit_lands_on_same_tissue():
    fit = h.PlanFit(400.0, 300.0, 90.0, 0.8, PIXEL_UM)
    fine = fit.rescaled(8)
    for k, (r, c) in fit.points().items():
        rf, cf = fine.points()[k]
        assert (rf + 0.5) / 8 - 0.5 == pytest.approx(r)
        assert (cf + 0.5) / 8 - 0.5 == pytest.approx(c)


def test_anterior_faces_left_at_90_degrees():
    pts = h._plan_points(h.PLAN_MM, 0.0, 0.0, 90.0, 1.0, PIXEL_UM)
    assert pts[1][1] < pts[3][1]  # target 1 (y=5.0) left of target 3 (y=2.5)
    assert pts[1][0] > pts[4][0]  # +x (target 1) below -x (target 4)


def test_slot_to_target():
    assert [h.slot_to_target(s, False) for s in range(1, 7)] == [1, 2, 3, 4, 5, 6]
    assert [h.slot_to_target(s, True) for s in range(1, 7)] == [4, 5, 6, 1, 2, 3]
