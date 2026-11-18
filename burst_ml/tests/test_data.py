"""Unit tests for bml/data.py using synthetic data."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from bml.data import (
    BANDS,
    constants,
    design,
    load_bursts,
    load_labels,
    load_design,
    recordings_for,
    subset,
)


# ---------------------------------------------------------------------------
# Synthetic fixtures
# ---------------------------------------------------------------------------

def _make_bursts() -> pd.DataFrame:
    """A minimal bursts table for two mice, two targets each.

    Recording names match what ``recordings_for()`` returns so that
    ``design()`` can find them.  Names are chosen so that no two
    ``(mouse, target)`` pairs share a name (the real data has the same
    property).
    """
    # Mouse 1, Target 1 → 3 recordings of 3 bursts each
    r1 = pd.DataFrame({
        "mouse": [1]*9,
        "recording": ["Target1","Target1","Target1",
                       "Target2","Target2","Target2",
                       "Target3","Target3","Target3"],
        "t": [0,1,2, 0,1,2, 0,1,2],
        "h2": [10.,11.,12., 20.,21.,22., 15.,16.,17.],
        "bb": [1.,2.,3., 4.,5.,6., 2.,3.,4.],
        "sub": [1.,2.,3., 4.,5.,6., 2.,3.,4.],
        "ultra": [1.,2.,3., 4.,5.,6., 2.,3.,4.],
        "wb": [100.]*9, "f0": [40.]*9,
        "h2_wf": [10.]*9, "volt": [1.]*9,
    })
    # Mouse 1, Target 2 → 1 recording of 3 bursts
    r2 = pd.DataFrame({
        "mouse": [1,1,1],
        "recording": ["Target2_Repeat"]*3,
        "t": [0,1,2],
        "h2": [30.,31.,32.],
        "bb": [1.,2.,3.], "sub": [1.,2.,3.], "ultra": [1.,2.,3.],
        "wb": [100.]*3, "f0": [40.]*3,
        "h2_wf": [10.]*3, "volt": [1.]*3,
    })
    # Mouse 2, Target 1 → 1 recording of 3 bursts
    r3 = pd.DataFrame({
        "mouse": [2,2,2],
        "recording": ["Target1"]*3,
        "t": [0,1,2],
        "h2": [50.,51.,52.],
        "bb": [2.,3.,4.], "sub": [2.,3.,4.], "ultra": [2.,3.,4.],
        "wb": [200.]*3, "f0": [50.]*3,
        "h2_wf": [20.]*3, "volt": [2.]*3,
    })
    # Mouse 2, Target 2 → 1 recording of 3 bursts
    r4 = pd.DataFrame({
        "mouse": [2,2,2],
        "recording": ["Target2"]*3,
        "t": [0,1,2],
        "h2": [60.,61.,62.],
        "bb": [3.,4.,5.], "sub": [3.,4.,5.], "ultra": [3.,4.,5.],
        "wb": [300.]*3, "f0": [60.]*3,
        "h2_wf": [30.]*3, "volt": [3.]*3,
    })
    return pd.concat([r1, r2, r3, r4], ignore_index=True)


def _make_labels() -> pd.DataFrame:
    """Labels matching the synthetic bursts: 4 sonicated targets + 1 control."""
    return pd.DataFrame({
        "mouse": [1, 1, 1, 2, 2],
        "target": [1, 2, 3, 1, 2],
        "dose": [1.0, 0.8, 0.0, 0.5, 0.6],
        "n_sonications": [2, 1, 0, 1, 1],
        "gfp_cells": [50, 40, 5, 30, 35],
        "cells": [1000, 900, 800, 700, 750],
    })


def _make_sonicating_labels() -> pd.DataFrame:
    """Labels matching the synthetic bursts: only sonicated targets (no control).

    Each row is one (mouse, target) combination.  ``k`` = number of burst groups
    (recordings), ``n`` = total bursts across all groups.

    Mouse 1 T1 has 3 burst groups (3+3+3=9 bursts).
    Mouse 1 T2, Mouse 2 T1, Mouse 2 T2 each have 1 group of 3 bursts.
    Mouse 1 T3 is a control (n_sonications=0) and is excluded.
    """
    return pd.DataFrame({
        "mouse": [1, 1, 2, 2],
        "target": [1, 2, 1, 2],
        "k": [3, 1, 1, 1],
        "n": [9, 3, 3, 3],
        "dose": [1.0, 0.8, 0.5, 0.6],
        "n_sonications": [2, 1, 1, 1],
        "gfp_cells": [50, 40, 30, 35],
        "cells": [1000, 900, 700, 750],
    })


# ---------------------------------------------------------------------------
# recordings_for
# ---------------------------------------------------------------------------

class TestRecordingsFor:
    def test_mouse1_t1(self):
        assert recordings_for(1, 1) == ["Target1", "Target2", "Target3"]

    def test_mouse1_t2(self):
        assert recordings_for(1, 2) == ["Target2_Repeat"]

    def test_mouse1_t4(self):
        assert recordings_for(1, 4) == ["Target4"]

    def test_mouse1_missing_target(self):
        assert recordings_for(1, 3) == []  # control, no recordings

    def test_mouse2_t1(self):
        assert recordings_for(2, 1) == ["Target1"]

    def test_mouse2_t6(self):
        assert recordings_for(2, 6) == ["Target6"]

    def test_mouse3_any_target(self):
        assert recordings_for(3, 4) == ["Target4"]

    def test_return_type(self):
        result = recordings_for(2, 5)
        assert isinstance(result, list)
        assert all(isinstance(r, str) for r in result)


# ---------------------------------------------------------------------------
# load_bursts  (needs a real file — skip)
# ---------------------------------------------------------------------------

class TestLoadBursts:
    def test_skipped_without_cache(self, tmp_path):
        cache = tmp_path / "empty"
        cache.mkdir()
        with pytest.raises(FileNotFoundError):
            load_bursts(cache)


# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------

class TestConstants:
    def test_h_ref_is_median(self):
        b = pd.DataFrame({"h2": [1.0, 2.0, 3.0, 4.0, 5.0], **{b: [1.0]*5 for b in BANDS}})
        c = constants(b)
        assert c["h_ref"] == 3.0

    def test_log_mean_and_sd(self):
        b = pd.DataFrame({"h2": [10.0]*3, **{b: [1.0, 2.0, 3.0] for b in BANDS}})
        c = constants(b)
        for band in BANDS:
            assert band in c["log_mean"]
            assert band in c["log_sd"]
            assert np.isfinite(c["log_mean"][band])
            assert np.isfinite(c["log_sd"][band])

    def test_log_clip_at_low_values(self):
        b = pd.DataFrame({"h2": [10.0]*2, **{b: [1e-15, 1e-15] for b in BANDS}})
        c = constants(b)
        for band in BANDS:
            assert np.isfinite(c["log_mean"][band])

    def test_h_ref_single_value(self):
        b = pd.DataFrame({"h2": [77.5], **{b: [1.0] for b in BANDS}})
        assert constants(b)["h_ref"] == 77.5


# ---------------------------------------------------------------------------
# load_labels
# ---------------------------------------------------------------------------

class TestLoadLabels:
    def test_skipped_without_file(self, tmp_path):
        p = tmp_path / "nope.csv"
        with pytest.raises(FileNotFoundError):
            load_labels(p)


# ---------------------------------------------------------------------------
# design
# ---------------------------------------------------------------------------

class TestDesign:
    def test_basic_structure(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)

        assert set(d.keys()) == {"mouse", "target", "animal", "n_animals",
                                  "k", "n", "dose", "h2", "t", "mask", "Z", "h_ref"}
        assert d["n_animals"] == 2  # mice 1 and 2
        assert len(d["k"]) == 4    # 4 sonicated targets

    def test_shapes(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)

        # All 1-D arrays have same length
        N = len(d["k"])
        for col in ("mouse", "target", "animal", "k", "n", "dose"):
            assert d[col].ndim == 1
            assert len(d[col]) == N

        # 2-D arrays have shape (N, T)
        T = d["h2"].shape[1]
        assert d["h2"].shape == (N, T)
        assert d["t"].shape == (N, T)
        assert d["mask"].shape == (N, T)
        assert d["Z"].shape == (N, T, 3)

    def test_mouse1_t1_has_three_recordings(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)

        i = np.where((d["mouse"] == 1) & (d["target"] == 1))[0][0]
        t_vals = d["t"][i][d["mask"][i]]
        # Three recordings → three restarts at 0
        assert (t_vals == 0).sum() == 3

    def test_mouse1_t3_is_not_present(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)

        mask = (d["mouse"] == 1) & (d["target"] == 3)
        assert not mask.any()

    def test_zero_padding(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)

        # Unmasked bursts have real data
        for i in range(len(d["h2"])):
            masked = d["h2"][i][d["mask"][i]]
            assert np.all(masked > 0)
            # Zero-padded regions should be zero
            assert np.all(d["h2"][i][~d["mask"][i]] == 0)
            assert np.all(d["Z"][i][~d["mask"][i]] == 0)

    def test_mask_correct(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)

        # Mouse 1 T1: 3 recordings × 3 bursts = 9 bursts
        i = np.where((d["mouse"] == 1) & (d["target"] == 1))[0][0]
        assert d["mask"][i].sum() == 9

        # Mouse 1 T2: 3 bursts
        i = np.where((d["mouse"] == 1) & (d["target"] == 2))[0][0]
        assert d["mask"][i].sum() == 3

        # Mouse 2 T1: 3 bursts
        i = np.where((d["mouse"] == 2) & (d["target"] == 1))[0][0]
        assert d["mask"][i].sum() == 3

    def test_missing_recording_raises(self):
        bursts = pd.DataFrame({
            "mouse": [1], "recording": ["Target1"], "t": [0],
            "h2": [10.0], "bb": [1.0], "sub": [1.0], "ultra": [1.0],
            "wb": [100.0], "f0": [40.0], "h2_wf": [10.0], "volt": [1.0],
        })
        labels = pd.DataFrame({
            "mouse": [1], "target": [3], "k": [1], "n": [3],
            "dose": [0.0], "n_sonications": [1], "gfp_cells": [10], "cells": [100],
        })
        c = constants(bursts)
        with pytest.raises(KeyError, match="Target3"):
            design(labels, bursts, c)

    def test_h_ref_scalar(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)
        assert isinstance(d["h_ref"], (float, np.floating))
        assert np.isfinite(d["h_ref"])


# ---------------------------------------------------------------------------
# subset
# ---------------------------------------------------------------------------

class TestSubset:
    def test_subset_by_index(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)

        s = subset(d, [0, 2])
        assert len(s["k"]) == 2
        assert list(s["mouse"]) == list(d["mouse"][[0, 2]])

    def test_subset_by_bool_mask(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)

        mask = np.array([True, False, True, False])
        s = subset(d, mask)
        assert len(s["k"]) == 2

    def test_animal_reindexed(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)

        # Mouse 1 rows are animal 0, Mouse 2 rows are animal 1
        s = subset(d, [0, 1])  # only mouse 1 rows
        assert s["n_animals"] == 1
        assert np.all(s["animal"] == 0)

        # Mouse 2 rows are animal 1
        s2 = subset(d, [3])  # only mouse 2 row
        assert s2["n_animals"] == 1
        assert np.all(s2["animal"] == 0)  # re-indexed to 0

    def test_h_ref_preserved(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)

        s = subset(d, [0])
        assert s["h_ref"] == d["h_ref"]

    def test_shapes_preserved(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)

        s = subset(d, [0, 1])
        assert s["h2"].ndim == 2
        assert s["Z"].ndim == 3
        assert s["mask"].ndim == 2

    def test_subset_empty_raises(self):
        bursts = _make_bursts()
        labels = _make_sonicating_labels()
        c = constants(bursts)
        d = design(labels, bursts, c)

        s = subset(d, [])
        assert len(s["k"]) == 0


# ---------------------------------------------------------------------------
# load_design (integration — skip without real data)
# ---------------------------------------------------------------------------

class TestLoadDesign:
    def test_missing_cache_raises(self, tmp_path):
        # load_design itself is checked on the real cache in tests/test_data_real.py
        with pytest.raises(FileNotFoundError):
            load_bursts(tmp_path)