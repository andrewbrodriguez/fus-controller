"""Per-burst acoustic emission feature extraction for FUS-BBB recordings.

This is a Python port of ``reference/nt_ExtractHarmonicData.m`` (N. Todd).
It reads one ``Mouse_Cntr_XX_TargetYY.mat`` recording and reduces every burst
to four scalar cavitation metrics.

Background
----------
During a sonication the FUS transducer fires a burst roughly once per second.
Microbubbles circulating in the blood oscillate in the acoustic field and
re-radiate sound. Two features of that re-radiated sound matter:

* **Second harmonic** (2x the carrier). Produced by *stable* cavitation --
  bubbles oscillating steadily. This correlates with BBB opening and is what
  the lab's real-time controller regulates.
* **Wideband / broadband** emission. Produced by *inertial* cavitation --
  bubbles collapsing violently. This is the damage signal, and the controller
  backs off when it appears.

The acquisition hardware has already applied the FFT, so each burst arrives as
a magnitude spectrum rather than a time series (see `load_recording`). Feature
extraction is therefore band selection and integration, not transformation.

What the MATLAB script does, in order
-------------------------------------
1. Load the ``.mat``.
2. Pull the carrier frequency, sample rate, frequency axis, and PRF.
3. Place two frequency windows: one on the second harmonic, one on the
   wideband monitor band.
4. For every burst, take peak and integrated magnitude inside each window.
5. Divide each series by its mean over the first few (pre-microbubble) bursts.
6. Cumulatively sum the normalised series to get a "dose".

This module reproduces steps 1-6. Plotting (steps 7-8 of the script) lives in
`fus.plots`.

Correspondence with the MATLAB variables
----------------------------------------
========================  =====================================
MATLAB                    Python
========================  =====================================
``peak_2nd``              ``BurstMetrics.peak_2nd``
``area_2nd``              ``BurstMetrics.area_2nd``
``peak_wb``               ``BurstMetrics.peak_wb``
``area_wb``               ``BurstMetrics.area_wb``
``*_norm``                ``BurstMetrics.normalised()``
``cum_*``                 ``BurstMetrics.cumulative()``
``V1_burst`` / ``V2_burst``  ``Recording.drive_voltage[:, 0/1]``
``t_burst``               ``Recording.burst_times``
========================  =====================================

Usage
-----
>>> from fus.extract import extract
>>> ex = extract("data/acoustic/20260611/Mouse_Cntr_01_Target1.mat")
>>> ex.summary()
>>> ex.metrics.cumulative().area_2nd[-1]      # cumulative harmonic dose
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import scipy.io

# ---------------------------------------------------------------------------
# Defaults, carried over from nt_ExtractHarmonicData.m
# ---------------------------------------------------------------------------

#: Centre of the wideband (broadband / inertial cavitation) monitor band, Hz.
#: Sits ~26 kHz above the second harmonic at the 837 kHz carrier -- close
#: enough to share the transducer's sensitivity, far enough not to pick up
#: the harmonic peak itself.
WIDEBAND_CENTRE_HZ = 1_700_000.0

#: Half-width of the analysis windows, in FFT bins. At df = 19.07 Hz these are
#: +/-381 Hz (harmonic) and +/-1.14 kHz (wideband). The harmonic window is
#: deliberately narrow: the second harmonic is a sharp spectral line, while
#: broadband emission is diffuse and needs a wider window to integrate.
HARMONIC_HALF_WIN_BINS = 20
WIDEBAND_HALF_WIN_BINS = 60

#: MATLAB's ``eps``, used to keep log10 finite at zero.
EPS = np.finfo(float).eps

#: ``Mouse_Controller_Data.xlsx`` stores ``Cumulative 2nd Harmonic`` scaled down
#: by this factor relative to the raw cumulative spectral area computed here.
#: Verified exact against 15 targets from the 2026-06-11 session.
XLSX_HARMONIC_SCALE = 1e4


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def _first(value: Any) -> Any:
    """Return the first element of an array-like, or the value itself.

    ``scipy.io.loadmat(squeeze_me=True)`` collapses 1x1 MATLAB arrays to bare
    scalars but leaves genuine vectors alone, so fields that are sometimes
    scalar and sometimes a vector need normalising.
    """
    arr = np.atleast_1d(value)
    return arr.flat[0] if arr.size else None


def _drive_voltage(posx: Any, n_bursts: int) -> np.ndarray | None:
    """Extract the per-burst drive voltage log, shape ``(n_bursts, 2)``.

    ``data.posx.data`` is a list of sonication log entries. The MATLAB script
    takes ``data.posx.data(end)`` -- the last one. Earlier entries can be stubs
    holding a 2-element placeholder rather than a real per-burst log, so this
    scans backwards for the first entry with a genuine 2-D voltage array.

    Channels are the two amplifier outputs (``AWG1``, ``AWG2``).
    """
    entries = np.atleast_1d(getattr(posx, "data", []))
    for entry in reversed(entries):
        volts = getattr(entry, "V", None)
        if volts is None:
            continue
        volts = np.atleast_2d(np.asarray(volts, dtype=float))
        if volts.ndim == 2 and volts.shape[0] > 1:
            if volts.shape[0] != n_bursts:
                # Baseline (``_BL``) files log fewer voltage rows than they
                # store bursts. Not fatal, but the two series are then not
                # burst-aligned and must not be plotted against each other.
                warnings.warn(
                    f"drive voltage has {volts.shape[0]} rows but there are "
                    f"{n_bursts} bursts; series are not burst-aligned",
                    stacklevel=2,
                )
            return volts
    return None


@dataclass
class Recording:
    """One sonication recording, loaded from a ``.mat`` file.

    Attributes
    ----------
    spectra
        ``(n_bursts, n_freq)`` array of per-burst **magnitude spectra**. Note
        these arrive already transformed -- the acquisition software runs the
        FFT, so the stored values are real, non-negative, and in the frequency
        domain. The MATLAB script's ``abs()`` is therefore a no-op kept for
        safety, and mirrored here.
    freq_axis
        ``(n_freq,)`` frequency axis in Hz, spanning 0 to Nyquist.
    carrier_hz
        Drive frequency of AWG1. 837 kHz for the controller study, so the
        second harmonic lands at 1.674 MHz.
    harmonic_goal
        The controller's ``[lower, upper]`` deadband on normalised second-
        harmonic amplitude. Drive voltage is raised below the lower bound and
        lowered above the upper. **The nominal setpoint is the midpoint** --
        see `harmonic_setpoint`.
    current_harmonic_goal
        The deadband actually in force at the end of the run. If this is below
        `harmonic_goal`, the wideband safety interlock fired mid-sonication and
        reduced the target.
    n_baseline
        Bursts fired before microbubbles reached the brain, used to normalise.
        Read from the file's ``DummyBursts`` field rather than hardcoded.
    """

    path: Path
    spectra: np.ndarray
    freq_axis: np.ndarray
    sample_rate_hz: float
    carrier_hz: float
    prf_hz: float
    n_baseline: int
    harmonic_goal: np.ndarray | None
    current_harmonic_goal: np.ndarray | None
    programmed_bursts: int | None
    control_gain: float | None
    drive_voltage: np.ndarray | None
    raw: Any

    # -- derived quantities -------------------------------------------------

    @property
    def n_bursts(self) -> int:
        return self.spectra.shape[0]

    @property
    def df_hz(self) -> float:
        """FFT bin width. ~19.07 Hz for a 5 MHz sample rate over 131072 bins."""
        return float(self.freq_axis[1] - self.freq_axis[0])

    @property
    def second_harmonic_hz(self) -> float:
        return 2.0 * self.carrier_hz

    @property
    def burst_times(self) -> np.ndarray:
        """Burst onset times in seconds, from the pulse repetition frequency."""
        return np.arange(self.n_bursts) / self.prf_hz

    @property
    def harmonic_setpoint(self) -> float | None:
        """Midpoint of the controller deadband -- the nominal 'Harmonic Goal'.

        This is the value to match against the ``Harmonic Goal`` column of
        ``Mouse_Controller_Data.xlsx``. A deadband of ``[0.7, 0.8]`` in the
        file corresponds to a setpoint of 0.75 in the spreadsheet.
        """
        if self.harmonic_goal is None:
            return None
        goal = np.atleast_1d(self.harmonic_goal).astype(float)
        if not np.all(np.isfinite(goal)):
            return None  # baseline runs store [inf inf]: no control active
        return float(np.mean(goal))

    @property
    def goal_was_reduced(self) -> bool:
        """True if the wideband interlock lowered the goal during the run."""
        if self.harmonic_goal is None or self.current_harmonic_goal is None:
            return False
        a = np.atleast_1d(self.harmonic_goal).astype(float)
        b = np.atleast_1d(self.current_harmonic_goal).astype(float)
        if a.shape != b.shape or not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
            return False
        return bool(np.any(b < a - 1e-9))


def load_recording(path: str | Path) -> Recording:
    """Load a ``Mouse_Cntr_XX_*.mat`` acoustic emission recording.

    These are MATLAB **v5** files despite their size (~0.5 GB), so
    ``scipy.io.loadmat`` reads them directly; ``h5py`` will not. The whole file
    is pulled into memory, which takes a few seconds -- reduce once and cache
    the metrics rather than reloading for every analysis.

    Parameters
    ----------
    path
        Path to the ``.mat`` file.

    Returns
    -------
    Recording
    """
    path = Path(path)
    mat = scipy.io.loadmat(path, struct_as_record=False, squeeze_me=True)
    if "data" not in mat:
        raise KeyError(
            f"{path.name} has no top-level 'data' struct; keys are "
            f"{[k for k in mat if not k.startswith('__')]}"
        )
    data = mat["data"]
    niscope, posx, awg = data.niscope, data.posx, data.AWG

    # data.niscope.ff is an array of per-burst structs; .data0 holds the
    # magnitude spectrum. Stack into (n_bursts, n_freq) once, up front.
    bursts = np.atleast_1d(niscope.ff)
    spectra = np.abs(np.vstack([np.asarray(b.data0, dtype=float) for b in bursts]))

    freq_axis = np.asarray(niscope.freqaxis, dtype=float).ravel()
    if spectra.shape[1] != freq_axis.size:
        raise ValueError(
            f"{path.name}: spectrum length {spectra.shape[1]} does not match "
            f"frequency axis length {freq_axis.size}"
        )

    # DummyBursts is the pre-microbubble burst count. The MATLAB script
    # hardcodes N_baseline = 5 with a note to adjust; every file checked so far
    # reports DummyBursts = 5, so sourcing it from the file is both faithful
    # and more robust.
    n_baseline = int(_first(getattr(posx, "DummyBursts", 5)) or 5)

    return Recording(
        path=path,
        spectra=spectra,
        freq_axis=freq_axis,
        sample_rate_hz=float(_first(niscope.SampleRate)),
        carrier_hz=float(_first(awg.AWG1.Frequency)),
        prf_hz=float(_first(posx.curPRF)),
        n_baseline=n_baseline,
        harmonic_goal=getattr(posx, "HarmGoal", None),
        current_harmonic_goal=getattr(posx, "curHarmGoal", None),
        programmed_bursts=(
            int(_first(posx.NBursts)) if hasattr(posx, "NBursts") else None
        ),
        control_gain=(
            float(_first(posx.PcontrolKp)) if hasattr(posx, "PcontrolKp") else None
        ),
        drive_voltage=_drive_voltage(posx, spectra.shape[0]),
        raw=data,
    )


# ---------------------------------------------------------------------------
# Frequency windows
# ---------------------------------------------------------------------------


@dataclass
class BandWindow:
    """A contiguous slice of FFT bins centred on a frequency of interest."""

    name: str
    centre_hz: float
    indices: np.ndarray
    freq_lo_hz: float
    freq_hi_hz: float

    def __len__(self) -> int:
        return int(self.indices.size)

    def __repr__(self) -> str:
        return (
            f"<BandWindow {self.name}: {self.freq_lo_hz / 1e6:.4f}-"
            f"{self.freq_hi_hz / 1e6:.4f} MHz, {len(self)} bins>"
        )


def band_window(
    freq_axis: np.ndarray, centre_hz: float, half_width_bins: int, name: str = ""
) -> BandWindow:
    """Build an index window of ``2 * half_width_bins + 1`` bins around a frequency.

    Mirrors the MATLAB::

        [~, idx] = min(abs(f_axis - f_centre));
        win = (idx - half_win) : (idx + half_win);

    with two differences worth knowing:

    * MATLAB indexes from 1, NumPy from 0. ``argmin`` handles this, but the
      slice end must be ``+ half + 1`` in Python to keep the same bin count.
    * MATLAB silently produces an invalid index if the window runs off the end
      of the axis. Here it is clipped, and the caller warned -- an asymmetric
      window biases the integrated area, so it should not pass unnoticed.
    """
    centre_idx = int(np.argmin(np.abs(freq_axis - centre_hz)))
    lo, hi = centre_idx - half_width_bins, centre_idx + half_width_bins + 1

    if lo < 0 or hi > freq_axis.size:
        warnings.warn(
            f"{name or 'window'} at {centre_hz / 1e6:.4f} MHz is clipped by the "
            f"edge of the frequency axis; integrated area will be biased low",
            stacklevel=2,
        )
        lo, hi = max(lo, 0), min(hi, freq_axis.size)

    indices = np.arange(lo, hi)
    return BandWindow(
        name=name,
        centre_hz=centre_hz,
        indices=indices,
        freq_lo_hz=float(freq_axis[indices[0]]),
        freq_hi_hz=float(freq_axis[indices[-1]]),
    )


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


@dataclass
class BurstMetrics:
    """Four per-burst scalars, each a ``(n_bursts,)`` array.

    ``peak_*`` is the tallest bin in the window; ``area_*`` is the window sum
    times the bin width -- a Riemann integral of *magnitude* (not power) over
    the band. Peak is the more sensitive detector of a narrow spectral line;
    area is more robust to small frequency drift of that line.
    """

    peak_2nd: np.ndarray
    area_2nd: np.ndarray
    peak_wb: np.ndarray
    area_wb: np.ndarray
    n_baseline: int

    def _baseline_mean(self, series: np.ndarray) -> float:
        return float(np.mean(series[: self.n_baseline]))

    def normalised(self) -> BurstMetrics:
        """Divide each series by its mean over the first ``n_baseline`` bursts.

        Those bursts are fired before microbubbles reach the brain, so they
        capture the response of *this* animal and *this* coupling with no
        cavitation present. Normalising by them removes per-animal differences
        in skull attenuation and transducer coupling, which is what makes
        values comparable across mice.

        A normalised value of 1.0 means "indistinguishable from baseline".
        """
        return BurstMetrics(
            peak_2nd=self.peak_2nd / self._baseline_mean(self.peak_2nd),
            area_2nd=self.area_2nd / self._baseline_mean(self.area_2nd),
            peak_wb=self.peak_wb / self._baseline_mean(self.peak_wb),
            area_wb=self.area_wb / self._baseline_mean(self.area_wb),
            n_baseline=self.n_baseline,
        )

    def cumulative(self) -> BurstMetrics:
        """Running sum of the normalised series -- the acoustic 'dose'.

        The final value of ``area_2nd`` here is the ``Cumulative 2nd Harmonic``
        column of ``Mouse_Controller_Data.xlsx``, and is the quantity plotted
        against delivered AAV.

        Because this is a sum over bursts, a long low-intensity sonication and
        a short intense one can reach the same total. Separating those two
        cases is exactly what the ``N Bursts`` x ``Harmonic Goal`` study design
        was built to test.
        """
        norm = self.normalised()
        return BurstMetrics(
            peak_2nd=np.cumsum(norm.peak_2nd),
            area_2nd=np.cumsum(norm.area_2nd),
            peak_wb=np.cumsum(norm.peak_wb),
            area_wb=np.cumsum(norm.area_wb),
            n_baseline=self.n_baseline,
        )


def compute_metrics(
    recording: Recording, harmonic: BandWindow, wideband: BandWindow
) -> BurstMetrics:
    """Reduce every burst to peak and integrated magnitude in each band.

    The MATLAB script loops burst by burst; here the whole ``(n_bursts,
    n_freq)`` array is already in memory, so both metrics are one vectorised
    reduction along the frequency axis. Results are identical.
    """
    df = float(recording.freq_axis[1] - recording.freq_axis[0])
    seg_2nd = recording.spectra[:, harmonic.indices]
    seg_wb = recording.spectra[:, wideband.indices]

    return BurstMetrics(
        peak_2nd=seg_2nd.max(axis=1),
        area_2nd=seg_2nd.sum(axis=1) * df,
        peak_wb=seg_wb.max(axis=1),
        area_wb=seg_wb.sum(axis=1) * df,
        n_baseline=recording.n_baseline,
    )


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------


@dataclass
class Extraction:
    """Everything `extract` produces for one recording."""

    recording: Recording
    harmonic_window: BandWindow
    wideband_window: BandWindow
    metrics: BurstMetrics

    @property
    def harmonic_to_wideband_ratio(self) -> float:
        """Cumulative harmonic area over cumulative wideband area.

        A crude stable-vs-inertial cavitation index: high means the exposure
        stayed in the useful regime, low means broadband collapse made up more
        of the emission. The MATLAB script prints this as
        ``Harmonic/WB area ratio``.
        """
        cum = self.metrics.cumulative()
        return float(cum.area_2nd[-1] / cum.area_wb[-1])

    def to_dict(self) -> dict[str, Any]:
        """Flat scalar summary -- one row per recording, ready for a DataFrame.

        This is the shape you want when joining many recordings against
        ``Mouse_Controller_Data.xlsx`` and, eventually, GFP measurements.

        Two columns are reported in the spreadsheet's own units so they can be
        compared directly (both verified exact on 15 targets from 2026-06-11):

        * ``cum_2nd_harmonic_xlsx`` -- the cumulative harmonic area divided by
          :data:`XLSX_HARMONIC_SCALE`. The spreadsheet's ``Cumulative 2nd
          Harmonic`` column carries that fixed 1e4 scaling; the unscaled value
          is what the lab's own plots use as "Cumulative 2nd Harmonic AUC".
        * ``mean_voltage`` -- the sum of drive voltage over the **sonicating**
          bursts only, divided by their count. The first ``n_baseline`` bursts
          are fired at zero volts, so including them would drag the mean down.
          ``cum_voltage`` is the plain sum and needs no such correction.
        """
        rec, cum = self.recording, self.metrics.cumulative()
        norm = self.metrics.normalised()
        return {
            "file": rec.path.name,
            "n_bursts": rec.n_bursts,
            "programmed_bursts": rec.programmed_bursts,
            "prf_hz": rec.prf_hz,
            "carrier_hz": rec.carrier_hz,
            "harmonic_setpoint": rec.harmonic_setpoint,
            "goal_was_reduced": rec.goal_was_reduced,
            "n_baseline": rec.n_baseline,
            "mean_2nd_harmonic": float(np.mean(norm.area_2nd)),
            "cum_2nd_harmonic_area": float(cum.area_2nd[-1]),
            "cum_2nd_harmonic_xlsx": float(cum.area_2nd[-1]) / XLSX_HARMONIC_SCALE,
            "cum_2nd_harmonic_peak": float(cum.peak_2nd[-1]),
            "cum_wideband_area": float(cum.area_wb[-1]),
            "cum_wideband_peak": float(cum.peak_wb[-1]),
            "harmonic_wideband_ratio": self.harmonic_to_wideband_ratio,
            "mean_voltage": self.mean_drive_voltage,
            "cum_voltage": (
                float(np.sum(rec.drive_voltage[:, 0]))
                if rec.drive_voltage is not None
                else None
            ),
        }

    @property
    def mean_drive_voltage(self) -> float | None:
        """Mean AWG1 voltage over sonicating bursts, matching the spreadsheet.

        The pre-microbubble bursts are fired at 0 V, so a plain mean over all
        bursts undercounts. Averaging over the non-zero bursts reproduces the
        ``Mean Voltage`` column exactly.
        """
        volts = self.recording.drive_voltage
        if volts is None:
            return None
        v1 = volts[:, 0]
        active = v1[v1 > 0]
        return float(active.mean()) if active.size else 0.0

    def summary(self) -> str:
        """Human-readable report, following the script's fprintf blocks."""
        rec = self.recording
        cum = self.metrics.cumulative()
        lines = [
            "=" * 62,
            f"  {rec.path.name}",
            "=" * 62,
            f"  Sample rate:       {rec.sample_rate_hz / 1e6:.2f} MHz",
            f"  Freq resolution:   {rec.df_hz:.4f} Hz",
            f"  FFT length:        {rec.freq_axis.size} pts",
            f"  Nyquist:           {rec.freq_axis[-1] / 1e6:.2f} MHz",
            f"  Bursts:            {rec.n_bursts}"
            + (
                f"  (programmed {rec.programmed_bursts})"
                if rec.programmed_bursts is not None
                else ""
            ),
            f"  PRF:               {rec.prf_hz:.1f} Hz",
            f"  Duration:          {rec.burst_times[-1]:.0f} s",
            f"  Carrier:           {rec.carrier_hz / 1e6:.4f} MHz",
            f"  2nd harmonic:      {rec.second_harmonic_hz / 1e6:.4f} MHz",
            f"  Wideband monitor:  {self.wideband_window.centre_hz / 1e6:.4f} MHz",
            f"  Baseline bursts:   1 to {rec.n_baseline}",
        ]

        if rec.harmonic_setpoint is not None:
            goal = np.atleast_1d(rec.harmonic_goal).astype(float)
            lines.append(
                f"  Harmonic goal:     {rec.harmonic_setpoint:.2f} "
                f"(deadband {goal[0]:.3f}-{goal[-1]:.3f})"
            )
        else:
            lines.append("  Harmonic goal:     none (control inactive)")
        if rec.control_gain is not None:
            lines.append(f"  Control gain Kp:   {rec.control_gain:g}")
        if rec.goal_was_reduced:
            cur = np.atleast_1d(rec.current_harmonic_goal).astype(float)
            lines.append(
                f"  !! goal reduced mid-run to {cur[0]:.3f}-{cur[-1]:.3f} "
                "(wideband interlock fired)"
            )

        lines += [
            "-" * 62,
            f"  2f window:         {self.harmonic_window.freq_lo_hz / 1e6:.4f}"
            f" - {self.harmonic_window.freq_hi_hz / 1e6:.4f} MHz"
            f" ({len(self.harmonic_window)} pts)",
            f"  WB window:         {self.wideband_window.freq_lo_hz / 1e6:.4f}"
            f" - {self.wideband_window.freq_hi_hz / 1e6:.4f} MHz"
            f" ({len(self.wideband_window)} pts)",
            "-" * 62,
            "  Cumulative metrics (normalised to baseline)",
            f"  2nd harmonic - cumulative area: {cum.area_2nd[-1]:.2f}"
            f"   (xlsx units: {cum.area_2nd[-1] / XLSX_HARMONIC_SCALE:.4f})",
            f"  2nd harmonic - cumulative peak: {cum.peak_2nd[-1]:.2f}",
            f"  Wideband     - cumulative area: {cum.area_wb[-1]:.2f}",
            f"  Wideband     - cumulative peak: {cum.peak_wb[-1]:.2f}",
            f"  Harmonic/WB area ratio:         "
            f"{self.harmonic_to_wideband_ratio:.2f}",
        ]

        if rec.drive_voltage is not None:
            v1 = rec.drive_voltage[:, 0]
            lines.append(
                f"  AWG1 voltage range:             "
                f"{v1.min():.3f} to {v1.max():.3f} V"
            )
            lines.append(
                f"  AWG1 mean / cumulative:         "
                f"{self.mean_drive_voltage:.4f} V / {v1.sum():.3f} V"
            )
        lines.append("=" * 62)
        return "\n".join(lines)


def extract(
    path: str | Path,
    *,
    wideband_centre_hz: float = WIDEBAND_CENTRE_HZ,
    harmonic_half_win: int = HARMONIC_HALF_WIN_BINS,
    wideband_half_win: int = WIDEBAND_HALF_WIN_BINS,
    n_baseline: int | None = None,
) -> Extraction:
    """Load a recording and reduce it to per-burst cavitation metrics.

    Parameters
    ----------
    path
        Path to a ``Mouse_Cntr_XX_TargetYY.mat`` file.
    wideband_centre_hz
        Centre of the broadband monitor band. Defaults to 1.7 MHz, as in the
        MATLAB script.
    harmonic_half_win, wideband_half_win
        Window half-widths in FFT bins.
    n_baseline
        Override the number of pre-microbubble bursts used for normalisation.
        By default this is read from the file's ``DummyBursts`` field.

    Returns
    -------
    Extraction
    """
    rec = load_recording(path)
    if n_baseline is not None:
        rec.n_baseline = int(n_baseline)
    if rec.n_baseline >= rec.n_bursts:
        raise ValueError(
            f"{rec.path.name}: n_baseline={rec.n_baseline} but the recording "
            f"only has {rec.n_bursts} bursts"
        )

    harmonic = band_window(
        rec.freq_axis, rec.second_harmonic_hz, harmonic_half_win, "2nd harmonic"
    )
    wideband = band_window(
        rec.freq_axis, wideband_centre_hz, wideband_half_win, "wideband"
    )

    # The two windows must not touch, or wideband "broadband" energy would
    # partly be the harmonic peak itself and the ratio would be meaningless.
    if np.intersect1d(harmonic.indices, wideband.indices).size:
        warnings.warn(
            "harmonic and wideband windows overlap; the harmonic/wideband "
            "ratio is not interpretable",
            stacklevel=2,
        )

    return Extraction(
        recording=rec,
        harmonic_window=harmonic,
        wideband_window=wideband,
        metrics=compute_metrics(rec, harmonic, wideband),
    )
