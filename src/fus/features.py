"""Acoustic features beyond the cumulative 2nd-harmonic AUC, from the raw waveforms.

Each recording stores the received time-domain waveform for every burst
(``data.niscope.wf[i].data0``: 150,000 samples at 5 MHz, 30 ms). The lab's script
and :mod:`fus.extract` only use the scope's precomputed spectrum in two narrow
windows. Here every burst is re-transformed and integrated over the bands that
carry distinct physics:

========  ===============  =============================================
band      frequency        what it indicates
========  ===============  =============================================
``f0``    837 kHz          received drive (transmission through skull)
``h2``    2 f0             stable, nonlinear bubble oscillation
``sub``   f0 / 2           subharmonic: stronger, period-doubled oscillation
``ultra`` 3 f0/2, 5 f0/2   ultraharmonics: as above
``bb``    between lines    broadband: inertial cavitation (bubble collapse)
========  ===============  =============================================

Each band's per-burst energy is divided by its mean over the zero-drive
baseline bursts (the receiver's noise floor), as :mod:`fus.extract` does for the
harmonic. Target-level features are then fixed in :data:`FEATURES` -- chosen from
the physics before looking at delivery, so they can be screened honestly.

>>> from fus.features import target_features
>>> target_features(["data/acoustic/20260611/Mouse_Cntr_02_Target4.mat"])
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from .extract import extract

#: Half-width of each harmonic line's integration band. The line is a few bins
#: wide at 33 Hz resolution; AWG2 runs 31 Hz above AWG1, so both sit inside.
LINE_HALF_WIDTH_HZ = 500.0

#: Broadband is everything in this range that is at least this far from any
#: multiple of f0/2, so no harmonic or its skirt leaks in.
BROADBAND_RANGE_HZ = (0.3e6, 2.4e6)
BROADBAND_CLEARANCE_HZ = 20e3

#: Target-level features, with what each is meant to capture. ``h2_cum`` is the
#: same physics as the lab's cumulative AUC but recomputed from the waveforms
#: (energy, not magnitude); ``cum_2nd_harmonic`` is the lab's number itself.
FEATURES = {
    "cum_2nd_harmonic": "lab's cumulative 2nd-harmonic AUC (the baseline predictor)",
    "n_bursts": "sonicating bursts",
    "harmonic_goal": "controller setpoint, burst-weighted over recordings",
    "mean_voltage": "mean drive voltage over sonicating bursts",
    "cum_voltage": "summed drive voltage",
    "f0_cum": "summed received fundamental (fold over noise)",
    "h2_cum": "summed 2nd-harmonic energy (fold over noise)",
    "sub_cum": "summed subharmonic energy",
    "ultra_cum": "summed ultraharmonic energy (3f0/2 + 5f0/2)",
    "bb_cum": "summed broadband energy (inertial cavitation)",
    "h2_per_f0": "median 2nd harmonic / fundamental: nonlinearity per unit received drive",
    "h2_per_volt": "median 2nd-harmonic amplitude per volt of drive",
    "h2_late_early": "2nd harmonic, last third / first third of the run (bubble clearance)",
    "h2_cv": "burst-to-burst CV of the 2nd harmonic (stability)",
    "interlock": "1 if the wideband interlock lowered the goal in any recording",
}


@dataclass
class BurstBands:
    """Per-burst band energies for one recording, sonicating bursts only."""

    f0: np.ndarray
    h2: np.ndarray
    sub: np.ndarray
    ultra: np.ndarray
    bb: np.ndarray
    voltage: np.ndarray


def _masks(freqs: np.ndarray, f0: float) -> dict[str, np.ndarray]:
    def line(mult: float) -> np.ndarray:
        return np.abs(freqs - mult * f0) <= LINE_HALF_WIDTH_HZ

    halves = np.arange(1, int(freqs[-1] / (f0 / 2)) + 2) * f0 / 2
    clear = np.min(np.abs(freqs[:, None] - halves[None, :]), axis=1) > BROADBAND_CLEARANCE_HZ
    lo, hi = BROADBAND_RANGE_HZ
    return {"f0": line(1), "h2": line(2), "sub": line(0.5),
            "ultra": line(1.5) | line(2.5), "bb": clear & (freqs >= lo) & (freqs <= hi)}


def burst_bands(rec) -> BurstBands:
    """Band energies per sonicating burst, as fold over the zero-drive baseline bursts."""
    wf = [np.asarray(b.data0, dtype=float) for b in np.atleast_1d(rec.raw.niscope.wf)]
    n = min(len(w) for w in wf)
    x = np.vstack([w[:n] for w in wf])
    x -= x.mean(axis=1, keepdims=True)
    power = np.abs(np.fft.rfft(x * np.hanning(n), axis=1)) ** 2
    freqs = np.fft.rfftfreq(n, 1.0 / rec.sample_rate_hz)
    bands = {k: power[:, m].sum(axis=1) for k, m in _masks(freqs, rec.carrier_hz).items()}
    nb = rec.n_baseline
    folds = {k: v[nb:] / v[:nb].mean() for k, v in bands.items()}
    volt = rec.drive_voltage[nb:, 0] if rec.drive_voltage is not None else np.full(len(wf) - nb, np.nan)
    return BurstBands(voltage=volt, **folds)


def target_features(paths: Sequence[str | Path]) -> dict[str, float]:
    """Features for one target, from all recordings fired at it.

    Several recordings (Mouse 1 target 1) are combined: sums add, burst-level
    statistics run over the concatenated sonicating bursts in firing order.
    """
    exs = [extract(p) for p in paths]
    bands = [burst_bands(ex.recording) for ex in exs]
    cat = {k: np.concatenate([getattr(b, k) for b in bands]) for k in BurstBands.__dataclass_fields__}
    h2, f0, v = cat["h2"], cat["f0"], cat["voltage"]
    third = len(h2) // 3
    on = v > 0
    n_each = np.array([len(b.h2) for b in bands])
    goals = np.array([ex.recording.harmonic_setpoint for ex in exs], dtype=float)
    return {
        "cum_2nd_harmonic": sum(ex.to_dict()["cum_2nd_harmonic_xlsx"] for ex in exs),
        "n_bursts": float(len(h2)),
        "harmonic_goal": float(np.sum(goals * n_each) / n_each.sum()),
        "mean_voltage": float(v[on].mean()),
        "cum_voltage": float(v.sum()),
        "f0_cum": float(f0.sum()),
        "h2_cum": float(h2.sum()),
        "sub_cum": float(cat["sub"].sum()),
        "ultra_cum": float(cat["ultra"].sum()),
        "bb_cum": float(cat["bb"].sum()),
        "h2_per_f0": float(np.median(h2[on] / f0[on])),
        "h2_per_volt": float(np.median(np.sqrt(h2[on]) / v[on])),
        "h2_late_early": float(h2[-third:].mean() / h2[:third].mean()),
        "h2_cv": float(h2.std() / h2.mean()),
        "interlock": float(any(ex.recording.goal_was_reduced for ex in exs)),
    }
