"""The burst cache against the dose and against the earlier per-burst table.

Skipped when ``burst_ml/cache/`` hasn't been built. ``data/processed/acoustic/bursts.csv`` is
the per-burst table ``scripts/sonication_eda.py`` built independently (lab-window h2 only).
"""

import numpy as np
import pandas as pd
import pytest

from bml import CACHE, REPO

BURSTS, RECS = CACHE / "bursts.csv", CACHE / "recordings.csv"
EDA = REPO / "data/processed/acoustic/bursts.csv"

pytestmark = pytest.mark.skipif(not RECS.exists(), reason="burst cache not built")


@pytest.fixture(scope="module")
def cache():
    return pd.read_csv(BURSTS), pd.read_csv(RECS)


def test_columns(cache):
    b, r = cache
    assert list(b.columns) == ["mouse", "recording", "t", "h2", "wb", "f0", "h2_wf", "sub",
                               "ultra", "bb", "volt"]
    assert list(r.columns) == ["mouse", "recording", "n_baseline", "n_spec", "n_wf", "n_kept",
                               "dose", "h2_baseline_sum"]


def test_one_row_per_recording_and_burst(cache):
    b, r = cache
    assert not r.duplicated(["mouse", "recording"]).any()
    assert not b.duplicated(["mouse", "recording", "t"]).any()
    n = b.groupby(["mouse", "recording"]).size()
    assert (n == r.set_index(["mouse", "recording"])["n_kept"].loc[n.index]).all()


def test_t_counts_from_zero(cache):
    b, _ = cache
    for _, g in b.groupby(["mouse", "recording"]):
        assert (g["t"].to_numpy() == np.arange(len(g))).all()


def test_dose_identity(cache):
    b, r = cache
    s = b.groupby(["mouse", "recording"])["h2"].sum().rename("s")
    j = r.set_index(["mouse", "recording"]).join(s)
    full = j[j.n_kept == j.n_spec]
    assert len(full) > 0
    assert np.allclose((full.s + full.h2_baseline_sum) / 1e4, full.dose, rtol=0, atol=1e-9)


def test_bands_positive(cache):
    b, _ = cache
    for c in ["h2", "wb", "f0", "h2_wf", "sub", "ultra", "bb"]:
        assert (b[c] > 0).all(), c


@pytest.mark.skipif(not EDA.exists(), reason="sonication_eda cache not present")
def test_h2_matches_earlier_table(cache):
    b, _ = cache
    eda = pd.read_csv(EDA)
    eda = eda[~eda.baseline].copy()
    eda["t"] = eda.groupby(["mouse", "recording"]).cumcount()
    j = b.merge(eda[["mouse", "recording", "t", "h2_norm"]], on=["mouse", "recording", "t"])
    assert len(j) > 0
    assert np.allclose(j.h2, j.h2_norm, rtol=1e-9, atol=1e-9)
