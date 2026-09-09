"""Figures for acoustic emission recordings.

Ports sections 7 and 8 of ``reference/nt_ExtractHarmonicData.m``:

* `plot_cavitation_metrics` -- 2x2 grid of the four raw metrics in dB.
* `plot_cumulative_dose` -- 1x2 grid of the cumulative normalised metrics.

Plus one addition not in the MATLAB script, `plot_spectrum`, which draws a
single burst's spectrum with the two analysis windows marked. Use it to sanity
check window placement on a new dataset before trusting any of the numbers.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from .extract import EPS, Extraction

HARMONIC_COLOUR = "tab:blue"
WIDEBAND_COLOUR = "tab:red"


def plot_cavitation_metrics(ex: Extraction) -> Figure:
    """Four raw per-burst metrics against time, in dB, with baseline marked.

    .. note::
       The MATLAB script labels this axis "dB re baseline" but plots
       ``20*log10(peak_2nd)`` -- the **absolute** metric, with a dashed line
       drawn at the baseline level. The values are not baseline-referenced.
       That behaviour is reproduced faithfully here, with the axis relabelled
       "dB (absolute)" so the plot is not misread. For genuinely
       baseline-referenced values, see `plot_cumulative_dose`, which uses the
       normalised series.
    """
    rec, m = ex.recording, ex.metrics
    t = rec.burst_times
    bl = slice(0, rec.n_baseline)

    panels = [
        ("2nd harmonic - peak", m.peak_2nd, HARMONIC_COLOUR),
        ("2nd harmonic - spectral area", m.area_2nd, HARMONIC_COLOUR),
        ("Wideband - peak", m.peak_wb, WIDEBAND_COLOUR),
        ("Wideband - spectral area", m.area_wb, WIDEBAND_COLOUR),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for ax, (title, series, colour) in zip(axes.ravel(), panels):
        ax.plot(t, 20 * np.log10(series + EPS), color=colour, lw=1.5)
        ax.axhline(
            20 * np.log10(np.mean(series[bl]) + EPS),
            color="k",
            ls="--",
            lw=1,
            label="baseline",
        )
        ax.axvspan(
            t[0], t[rec.n_baseline - 1], color="k", alpha=0.06, lw=0
        )  # pre-microbubble window
        ax.set_title(title)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("dB (absolute)")
        ax.set_xlim(t[0], t[-1])
        ax.grid(alpha=0.3)
        ax.legend(loc="lower right", fontsize="small")

    fig.suptitle(
        f"{rec.path.name}  |  {rec.carrier_hz / 1e6:.4f} MHz carrier  |  "
        f"PRF {rec.prf_hz:.0f} Hz  |  {rec.n_bursts} bursts",
        fontsize=11,
    )
    return fig


def plot_cumulative_dose(ex: Extraction) -> Figure:
    """Cumulative normalised harmonic and wideband metrics against time.

    These curves are the acoustic "dose". Their endpoints are the values that
    go into ``Mouse_Controller_Data.xlsx`` and get regressed against delivered
    AAV, so the shape of the curve -- not just where it ends -- is worth
    reading: a straight line means a steady exposure, a knee means the
    controller changed behaviour partway through.
    """
    rec = ex.recording
    cum = ex.metrics.cumulative()
    t = rec.burst_times

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for ax, (title, harm, wide) in zip(
        axes,
        [
            ("Cumulative spectral area", cum.area_2nd, cum.area_wb),
            ("Cumulative peak", cum.peak_2nd, cum.peak_wb),
        ],
    ):
        ax.plot(t, harm, color=HARMONIC_COLOUR, lw=1.5, label="2nd harmonic (stable)")
        ax.plot(t, wide, color=WIDEBAND_COLOUR, lw=1.5, label="Wideband (inertial)")
        ax.set_title(title)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Cumulative sum (normalised)")
        ax.set_xlim(t[0], t[-1])
        ax.grid(alpha=0.3)
        ax.legend(fontsize="small")

    fig.suptitle(
        f"Cumulative cavitation dose  |  {rec.path.name}  |  "
        f"{rec.n_bursts} bursts",
        fontsize=11,
    )
    return fig


def plot_spectrum(ex: Extraction, burst: int = -1, span_khz: float = 120.0) -> Figure:
    """One burst's spectrum around the second harmonic, with windows shaded.

    Not part of the MATLAB script. This is the diagnostic to run first on any
    new recording: it shows whether the harmonic window actually sits on the
    harmonic peak, and whether the wideband window is clear of it. If the
    windows are misplaced, every downstream number is wrong in a way the
    summary statistics will not reveal.

    Parameters
    ----------
    burst
        Index of the burst to draw. Defaults to the last one, where emission
        is usually strongest.
    span_khz
        Half-width of the plotted frequency range around the second harmonic.
    """
    rec = ex.recording
    f_mhz = rec.freq_axis / 1e6
    centre = rec.second_harmonic_hz / 1e6
    span = span_khz / 1e3

    sel = (f_mhz >= centre - span) & (f_mhz <= centre + span)

    fig, ax = plt.subplots(figsize=(10, 4.5), constrained_layout=True)
    ax.semilogy(f_mhz[sel], rec.spectra[burst][sel] + EPS, color="0.3", lw=0.8)

    for window, colour, label in [
        (ex.harmonic_window, HARMONIC_COLOUR, "2nd harmonic window"),
        (ex.wideband_window, WIDEBAND_COLOUR, "wideband window"),
    ]:
        ax.axvspan(
            window.freq_lo_hz / 1e6,
            window.freq_hi_hz / 1e6,
            color=colour,
            alpha=0.25,
            lw=0,
            label=f"{label} ({len(window)} bins)",
        )

    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("Magnitude")
    ax.set_title(
        f"{rec.path.name}  |  burst {burst % rec.n_bursts + 1} of {rec.n_bursts}"
    )
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize="small")
    return fig
