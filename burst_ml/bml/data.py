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
    if mouse == 1:
        recs = MOUSE1.get(target, [])
        return list(recs)
    return [f"Target{target}"]


def load_bursts(cache=CACHE) -> pd.DataFrame:
    """``cache/bursts.csv`` as a DataFrame."""
    path = cache / "bursts.csv"
    df = pd.read_csv(path)
    df["mouse"] = df["mouse"].astype(int)
    df["t"] = df["t"].astype(int)
    return df


def constants(bursts: pd.DataFrame) -> dict:
    """Unsupervised constants from every row of ``bursts`` (all mice):
    ``{"h_ref": median h2, "log_mean": {band: mean log(band)}, "log_sd": {band: sd log(band)}}``
    for band in BANDS. Use ``np.log(np.clip(x, 1e-12, None))``."""
    h_ref = float(np.median(bursts["h2"]))
    log_mean: dict[str, float] = {}
    log_sd: dict[str, float] = {}
    for band in BANDS:
        x = np.log(np.clip(bursts[band].values, 1e-12, None))
        log_mean[band] = float(np.mean(x))
        log_sd[band] = float(np.std(x, ddof=0))
    return {"h_ref": h_ref, "log_mean": log_mean, "log_sd": log_sd}


def load_labels(labels=LABELS) -> pd.DataFrame:
    """Labels for the 35 sonicated targets: columns mouse, target, k, n, dose.
    Drop rows with ``n_sonications == 0`` (Mouse 1 T3, the control)."""
    df = pd.read_csv(labels)
    df = df[df["n_sonications"] != 0].copy()
    return df[["mouse", "target", "gfp_cells", "cells", "dose"]].rename(
        columns={"gfp_cells": "k", "cells": "n"}
    )


def design(labels: pd.DataFrame, bursts: pd.DataFrame, consts: dict) -> dict:
    """Build the design dict (module doc) for the rows of ``labels``.

    For each target, concatenate the bursts of ``recordings_for(mouse, target)`` in that order.
    Raise ``KeyError`` naming the recording if one is missing from ``bursts``."""
    rows: list[dict] = []
    for _, row in labels.iterrows():
        m = int(row["mouse"])
        t = int(row["target"])
        recs = recordings_for(m, t)
        parts: list[pd.DataFrame] = []
        missing: list[str] = []
        for rec in recs:
            chunk = bursts[(bursts["recording"] == rec) & (bursts["mouse"] == m)]
            if chunk.empty:
                missing.append(rec)
            else:
                parts.append(chunk)
        if not recs:
            raise KeyError(f"mouse {m} Target{t} has no recordings")
        if missing:
            raise KeyError(f"mouse {m}: {', '.join(missing)} not in the burst cache")
        target_bursts = pd.concat(parts, ignore_index=True)
        rows.append({
            "mouse": m,
            "target": t,
            "k": int(row["k"]),
            "n": int(row["n"]),
            "dose": float(row["dose"]),
            "h2": target_bursts["h2"].values.astype(np.float64),
            "t": target_bursts["t"].values.astype(np.float64),
            "bb": target_bursts["bb"].values.astype(np.float64),
            "sub": target_bursts["sub"].values.astype(np.float64),
            "ultra": target_bursts["ultra"].values.astype(np.float64),
        })

    # Find T (longest burst count)
    T = max(len(r["h2"]) for r in rows)

    h2_arr = np.zeros((len(rows), T), dtype=np.float64)
    t_arr = np.zeros((len(rows), T), dtype=np.float64)
    mask_arr = np.zeros((len(rows), T), dtype=bool)
    Z_arr = np.zeros((len(rows), T, 3), dtype=np.float64)

    for i, r in enumerate(rows):
        n = len(r["h2"])
        h2_arr[i, :n] = r["h2"]
        t_arr[i, :n] = r["t"]
        mask_arr[i, :n] = True
        log_bb = np.log(np.clip(r["bb"], 1e-12, None))
        log_sub = np.log(np.clip(r["sub"], 1e-12, None))
        log_ultra = np.log(np.clip(r["ultra"], 1e-12, None))
        vals_bb = np.zeros(T, dtype=np.float64)
        vals_sub = np.zeros(T, dtype=np.float64)
        vals_ultra = np.zeros(T, dtype=np.float64)
        vals_bb[:n] = (log_bb - consts["log_mean"]["bb"]) / consts["log_sd"]["bb"]
        vals_sub[:n] = (log_sub - consts["log_mean"]["sub"]) / consts["log_sd"]["sub"]
        vals_ultra[:n] = (log_ultra - consts["log_mean"]["ultra"]) / consts["log_sd"]["ultra"]
        Z_arr[i, :, 0] = vals_bb
        Z_arr[i, :, 1] = vals_sub
        Z_arr[i, :, 2] = vals_ultra

    mouse_arr = np.array([r["mouse"] for r in rows], dtype=np.int64)
    target_arr = np.array([r["target"] for r in rows], dtype=np.int64)
    k_arr = np.array([r["k"] for r in rows], dtype=np.int64)
    n_arr = np.array([r["n"] for r in rows], dtype=np.int64)
    dose_arr = np.array([r["dose"] for r in rows], dtype=np.float64)

    # Assign animal ids: consecutive 0..J-1 by mouse id order
    unique_mice = sorted(set(mouse_arr))
    mouse_to_animal = {m: i for i, m in enumerate(unique_mice)}
    animal_arr = np.array([mouse_to_animal[m] for m in mouse_arr], dtype=np.int64)

    return {
        "mouse": mouse_arr,
        "target": target_arr,
        "animal": animal_arr,
        "n_animals": len(unique_mice),
        "k": k_arr,
        "n": n_arr,
        "dose": dose_arr,
        "h2": h2_arr,
        "t": t_arr,
        "mask": mask_arr,
        "Z": Z_arr,
        "h_ref": consts["h_ref"],
    }


def subset(d: dict, rows) -> dict:
    """The design restricted to ``rows`` (index array or boolean mask), with ``animal``
    re-indexed to 0..J'-1 in order of first appearance and ``n_animals`` = J'.
    ``h_ref`` is carried over unchanged."""
    if d["mouse"].ndim == 1:
        d = {k: np.atleast_1d(v) if isinstance(v, np.ndarray) and v.ndim == 1 else v for k, v in d.items()}

    if isinstance(rows, (list, np.ndarray)) and np.issubdtype(np.asarray(rows).dtype, np.integer):
        idx = np.asarray(rows)
    else:
        idx = np.flatnonzero(np.asarray(rows))

    out: dict = {}
    for k in ("mouse", "target", "animal", "k", "n", "dose", "h2", "t", "mask", "Z"):
        val = d[k]
        if isinstance(val, np.ndarray) and val.ndim >= 2:
            out[k] = val[idx]
        else:
            out[k] = np.asarray(val)[idx]

    # Re-index animal by first appearance order in the subset
    animal_map: dict[int, int] = {}
    for a in out["animal"]:
        if a not in animal_map:
            animal_map[a] = len(animal_map)
    out["animal"] = np.array([animal_map[a] for a in out["animal"]], dtype=np.int64)
    out["n_animals"] = len(animal_map)
    out["h_ref"] = d["h_ref"]
    return out


def load_design() -> dict:
    """``design(load_labels(), bursts, constants(bursts))`` with ``bursts = load_bursts()``."""
    bursts = load_bursts()
    return design(load_labels(), bursts, constants(bursts))
