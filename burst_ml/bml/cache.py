"""Per-burst acoustic features for every recording, cached as CSV in ``burst_ml/cache/``.

Run from ``burst_ml/``::

    ../.venv/bin/python -m bml.cache                 # every recording, all 24 mice
    ../.venv/bin/python -m bml.cache --mice 2 --limit 1

Reads recordings with ``fus.archive.recordings()`` and ``fus.archive.local_path(member)``.
The zip is the only copy of the data: never unpack it, extract from it by hand, or write
anywhere under ``data/``.

Outputs (both appended one recording at a time and saved after each, so a killed run resumes;
recordings already in ``recordings.csv`` are skipped):

``cache/bursts.csv`` -- one row per **sonicating** burst (the first ``n_baseline`` bursts are
dropped). Columns:

    mouse, recording, t, h2, wb, f0, h2_wf, sub, ultra, bb, volt

* ``t``: 0-based burst index within the recording, counting sonicating bursts only.
* ``h2``, ``wb``: ``ex.metrics.normalised().area_2nd[nb:]`` and ``.area_wb[nb:]`` where
  ``ex = fus.extract.extract(path)`` and ``nb = ex.recording.n_baseline``.
* ``f0, h2_wf, sub, ultra, bb, volt``: the fields ``f0, h2, sub, ultra, bb, voltage`` of
  ``fus.features.burst_bands(ex.recording)`` (already sonicating bursts only).
* If the two sources differ in length, truncate both to the shorter.

``cache/recordings.csv`` -- one row per recording. Columns:

    mouse, recording, n_baseline, n_spec, n_wf, n_kept, dose, h2_baseline_sum

* ``n_spec`` = ``len(area_2nd) - nb``, ``n_wf`` = ``len(burst_bands(...).h2)``,
  ``n_kept`` = rows written to bursts.csv.
* ``dose`` = ``ex.to_dict()["cum_2nd_harmonic_xlsx"]``.
* ``h2_baseline_sum`` = ``normalised().area_2nd[:nb].sum()``. Then
  ``(bursts.h2.sum() + h2_baseline_sum) / 1e4 == dose`` for each recording when
  ``n_kept == n_spec``.

Suppress warnings while reading, as ``scripts/sonication_eda.py`` does.
"""

from __future__ import annotations

import argparse

from . import CACHE


def read_recording(member) -> tuple["pd.DataFrame", dict]:
    """(bursts rows, recordings row) for one ``fus.archive.Member``."""
    raise NotImplementedError


def build_cache(mice: list[int] | None = None, limit: int | None = None, cache=CACHE) -> None:
    """Fill ``cache/bursts.csv`` and ``cache/recordings.csv`` incrementally (see module doc)."""
    raise NotImplementedError


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mice", type=int, nargs="*")
    ap.add_argument("--limit", type=int, help="read at most this many new recordings")
    a = ap.parse_args(argv)
    build_cache(a.mice, a.limit)


if __name__ == "__main__":
    main()
