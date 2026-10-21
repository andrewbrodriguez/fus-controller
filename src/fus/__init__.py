"""FUS acoustic feedback controller -- analysis package.

Feature extraction for focused-ultrasound BBB opening experiments: reduces
raw per-burst acoustic emission spectra to cavitation metrics that can be
regressed against delivered AAV.

>>> from fus import extract
>>> from fus import archive
>>> with archive.paths(1, ["Target1"]) as (path,):
...     ex = extract(path)
>>> print(ex.summary())
"""

from .extract import (
    BandWindow,
    BurstMetrics,
    Extraction,
    Recording,
    band_window,
    compute_metrics,
    extract,
    load_recording,
)

__all__ = [
    "BandWindow",
    "BurstMetrics",
    "Extraction",
    "Recording",
    "band_window",
    "compute_metrics",
    "extract",
    "load_recording",
]
