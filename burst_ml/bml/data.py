"""Labels and per-burst arrays for the 35 sonicated targets of mice 1-6.

Everything a model needs is one ``dict`` of numpy arrays (the "design"), one row per target:

==========  ==========================  ==============================================
key         shape                        meaning
==========  ==========================  ==============================================
mouse       (N,) int                    mouse id (1-6)
target      (N,) int                    target / position (1-6)
animal      (N,) int                    0..J-1, consecutive over the mice present
n_animals   int                         J
k           (N,) int                    ΣGFP+ neurons (``gfp_cells``)
n           (N,) int                    ΣNeuN neurons (``cells``)
dose        (N,) float                  cumulative 2nd harmonic (labels CSV ``dose``)
h2          (N, T) float                per-burst ``h2``, zero-padded to T bursts
t           (N, T) float                per-burst ``t`` (restarts at 0 in each recording)
mask        (N, T) bool                 True where a real burst sits
Z           (N, T, 3) float             z-scored log ``bb``, ``sub``, ``ultra``; 0 where masked
h_ref       float                       median sonicating ``h2`` over all 24 mice
==========  ==========================  ==============================================

T is the longest target's burst count (Mouse 1 T1, three recordings, is the longest).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import CACHE, LABELS

BANDS = ("bb", "sub", "ultra")

#: Mouse 1 is the exception to TargetN = position N (lab, 2026-09-18; recording timestamps).
#: Position 3 was never sonicated (the control) and has no recordings.
MOUSE1 = {1: ["Target1", "Target2", "Target3"], 2: ["Target2_Repeat"],
          4: ["Target4"], 5: ["Target5"], 6: ["Target6"]}


def recordings_for(mouse: int, target: int) -> list[str]:
    """Recording names fired at ``target`` of ``mouse``, in firing order."""
    raise NotImplementedError


def load_bursts(cache=CACHE) -> pd.DataFrame:
    """``cache/bursts.csv`` as a DataFrame."""
    raise NotImplementedError


def constants(bursts: pd.DataFrame) -> dict:
    """Unsupervised constants from every row of ``bursts`` (all mice):
    ``{"h_ref": median h2, "log_mean": {band: mean log(band)}, "log_sd": {band: sd log(band)}}``
    for band in BANDS. Use ``np.log(np.clip(x, 1e-12, None))``."""
    raise NotImplementedError


def load_labels(labels=LABELS) -> pd.DataFrame:
    """Labels for the 35 sonicated targets: columns mouse, target, k, n, dose.
    Drop rows with ``n_sonications == 0`` (Mouse 1 T3, the control)."""
    raise NotImplementedError


def design(labels: pd.DataFrame, bursts: pd.DataFrame, consts: dict) -> dict:
    """Build the design dict (module doc) for the rows of ``labels``.

    For each target, concatenate the bursts of ``recordings_for(mouse, target)`` in that order.
    Raise ``KeyError`` naming the recording if one is missing from ``bursts``."""
    raise NotImplementedError


def subset(d: dict, rows) -> dict:
    """The design restricted to ``rows`` (index array or boolean mask), with ``animal``
    re-indexed to 0..J'-1 in order of first appearance and ``n_animals`` = J'.
    ``h_ref`` is carried over unchanged."""
    raise NotImplementedError


def load_design() -> dict:
    """``design(load_labels(), bursts, constants(bursts))`` with ``bursts = load_bursts()``."""
    bursts = load_bursts()
    return design(load_labels(), bursts, constants(bursts))
