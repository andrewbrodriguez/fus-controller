"""The design built from the real cache and labels (orchestrator's check; skipped without cache).

The key check: each target's burst sum reproduces the labels' dose, which only holds if every
recording was joined to the right target (Mouse 1's T1 = three recordings, T2 = the repeat).
"""

import numpy as np
import pytest

from bml import CACHE

pytestmark = pytest.mark.skipif(not (CACHE / "recordings.csv").exists(),
                                reason="burst cache not built")


@pytest.fixture(scope="module")
def d():
    from bml.data import load_design
    return load_design()


def test_35_sonicated_targets(d):
    assert len(d["k"]) == 35
    assert d["n_animals"] == 6
    assert not ((d["mouse"] == 1) & (d["target"] == 3)).any()


def test_burst_sums_reproduce_dose(d):
    s = np.where(d["mask"], d["h2"], 0).sum(1) / 1e4
    # the baseline bursts (~1 each, 5 per recording) are the only part of the dose left out
    assert np.allclose(s, d["dose"], atol=0.002), np.c_[d["mouse"], d["target"], s, d["dose"]]


def test_mouse1_t1_is_three_recordings(d):
    i = np.flatnonzero((d["mouse"] == 1) & (d["target"] == 1))[0]
    t = d["t"][i][d["mask"][i]]
    assert (t == 0).sum() == 3          # time restarts in each of the three recordings
    assert d["mask"][i].sum() == d["mask"].sum(1).max()


def test_padding_is_inert(d):
    assert (d["h2"][~d["mask"]] == 0).all()
    assert (d["Z"][~d["mask"]] == 0).all()
    assert np.isfinite(d["Z"]).all()
    assert 50 < d["h_ref"] < 120        # median sonicating h2 over all mice was 77.5
