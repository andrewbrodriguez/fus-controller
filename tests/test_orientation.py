"""Clicks -> plan angle and notch side, without opening napari."""

import numpy as np
import pytest

from fus import histology as h
from fus import orientation as o


@pytest.mark.parametrize("front_offset, angle", [
    ((0, -100), 90.0),   # front to the left: Mouse_01's layout
    ((100, 0), 0.0),     # down
    ((0, 100), 270.0),   # right
    ((-100, 0), 180.0),  # up
])
def test_anterior_angle_matches_plan_convention(front_offset, angle):
    centre = (500.0, 500.0)
    got = o.anterior_angle(centre, np.add(centre, front_offset))
    assert got == pytest.approx(angle)
    # and the plan's +y really points at the click
    dr, dc = h._mm_to_px(0.0, 1.0, got, 1.0, 10.0)
    assert np.sign(np.round([dr, dc], 6)).tolist() == np.sign(front_offset).tolist()


def test_notch_side_follows_plan_x():
    centre = (500.0, 500.0)
    # anterior left (90 deg): plan +x (targets 1-3) is the lower half of the image
    assert o.notch_on_plus_x(90.0, centre, (650.0, 700.0))
    assert not o.notch_on_plus_x(90.0, centre, (350.0, 700.0))
    # turned 180 degrees: +x is now the upper half
    assert o.notch_on_plus_x(270.0, centre, (350.0, 300.0))


def test_sheet_round_trip(tmp_path):
    a = o.Annotation(2, "slide01_s3", (10.0, 20.0), (5.0, 20.0), (15.0, 30.0), False, "x")
    b = o.Annotation(1, "section_s3", exclude=True)
    path = tmp_path / "sheet.csv"
    o.write_sheet({(2, a.section): a, (1, b.section): b}, path)
    back = o.read_sheet(path)
    assert back[(2, "slide01_s3")] == a
    assert back[(1, "section_s3")].exclude and not back[(1, "section_s3")].front
    assert back[(2, "slide01_s3")].angle_deg == pytest.approx(180.0)


def test_mouse_from_path():
    assert o.mouse_of("data/processed/histology/Mouse_02/ds4/slide01_s2.ome.tif") == 2


def test_notch_offset_right_angle():
    a = o.Annotation(2, "x", (500.0, 500.0), (500.0, 300.0), (700.0, 500.0))
    assert a.notch_offset_deg == pytest.approx(90.0)


def _m1(section, notch_row):
    # anterior left; notch clicked below (+x side) or above (-x side) the midline
    return o.Annotation(1, section, (500.0, 500.0), (500.0, 100.0), (notch_row, 500.0))


def test_plus_x_side_from_mouse01_clicks():
    # Mouse 1 notch is on the right. s2 is as scanned, s5 mirrored: if +x is the
    # animal's right, the notch shows below the midline in s2 and above in s5.
    sheet = {(1, "section_s2"): _m1("section_s2", 700.0),
             (1, "section_s5"): _m1("section_s5", 300.0)}
    assert o.plus_x_side(sheet) == "right"
    sheet[(1, "section_s5")] = _m1("section_s5", 700.0)  # contradicts s2
    with pytest.raises(ValueError):
        o.plus_x_side(sheet)
    assert o.plus_x_side({}) is None


def test_mouse02_mirroring_uses_left_notch():
    # +x = animal right; Mouse 2 notch is on the left. Notch below the midline
    # (on +x) means slots 1-3 are the animal's left, i.e. mirrored.
    below = o.Annotation(2, "s", (500.0, 500.0), (500.0, 100.0), (700.0, 500.0))
    above = o.Annotation(2, "s", (500.0, 500.0), (500.0, 100.0), (300.0, 500.0))
    assert o.is_mirrored(below, "right") and not o.is_mirrored(above, "right")


@pytest.mark.parametrize("front, notch, flips", [
    ((500.0, 100.0), (700.0, 500.0), False),  # front left, notch down: a pure rotation
    ((500.0, 100.0), (300.0, 500.0), True),   # front left, notch up: needs a mirror
    ((100.0, 500.0), (500.0, 100.0), False),  # already front up, notch left
])
def test_canonical_affine_puts_front_up_notch_left(front, notch, flips):
    a = o.Annotation(2, "s", (500.0, 500.0), front, notch)
    m = o.canonical_affine(a, scale=2.0)

    def world(p):
        return (m @ np.array([p[0] / 2, p[1] / 2, 1.0]))[:2]

    assert world(a.centre) == pytest.approx([0, 0])
    f, n = world(front), world(notch)
    assert f[1] == pytest.approx(0, abs=1e-9) and f[0] < 0   # straight up
    assert n[1] < 0 and abs(n[0]) < 1e-9                      # straight left (notch is square here)
    assert (np.linalg.det(m[:2, :2]) < 0) == flips
    assert abs(np.linalg.det(m[:2, :2])) == pytest.approx(1)  # no scaling: pixels untouched


def test_update_sheet_leaves_other_fields_and_sections(tmp_path):
    path = tmp_path / "sheet.csv"
    a = o.Annotation(2, "s1", (1.0, 2.0), (0.0, 2.0), (1.0, 5.0), review="wrong")
    b = o.Annotation(2, "s2", notes="keep me")
    o.write_sheet({(2, "s1"): a, (2, "s2"): b}, path)
    # an orientation window saving s1's clicks must not touch its verdict or s2
    o.update_sheet((2, "s1"), path, centre=(9.0, 9.0), exclude=True)
    back = o.read_sheet(path)
    assert back[(2, "s1")].review == "wrong" and back[(2, "s1")].centre == (9.0, 9.0)
    assert back[(2, "s2")].notes == "keep me"
