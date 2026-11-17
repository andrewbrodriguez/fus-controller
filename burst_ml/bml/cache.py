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
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import fus.archive
from fus.extract import extract
from fus.features import burst_bands

from . import CACHE


def read_recording(member) -> tuple[pd.DataFrame, dict]:
    """(bursts rows, recordings row) for one ``fus.archive.Member``."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with fus.archive.local_path(member) as path:
            ex = extract(path)

    rec = ex.recording
    nb = rec.n_baseline

    # --- per-burst features from extract (normalised area) ---
    norm = ex.metrics.normalised()
    h2 = norm.area_2nd[nb:]
    wb = norm.area_wb[nb:]

    # --- per-burst features from waveform FFT ---
    bands = burst_bands(rec)
    f0 = bands.f0
    h2_wf = bands.h2
    sub = bands.sub
    ultra = bands.ultra
    bb = bands.bb
    volt = bands.voltage

    # Truncate to the shorter source if they differ in length
    n = min(len(h2), len(f0))
    h2 = h2[:n]
    wb = wb[:n]
    f0 = f0[:n]
    h2_wf = h2_wf[:n]
    sub = sub[:n]
    ultra = ultra[:n]
    bb = bb[:n]
    volt = volt[:n]

    t = np.arange(n)

    bursts = pd.DataFrame({
        "mouse": member.mouse,
        "recording": member.recording,
        "t": t,
        "h2": h2,
        "wb": wb,
        "f0": f0,
        "h2_wf": h2_wf,
        "sub": sub,
        "ultra": ultra,
        "bb": bb,
        "volt": volt,
    })

    dose = ex.to_dict()["cum_2nd_harmonic_xlsx"]
    h2_baseline_sum = float(norm.area_2nd[:nb].sum())

    n_spec = len(norm.area_2nd) - nb
    n_wf = len(bands.h2)
    n_kept = int(n)

    rec_row = {
        "mouse": member.mouse,
        "recording": member.recording,
        "n_baseline": nb,
        "n_spec": n_spec,
        "n_wf": n_wf,
        "n_kept": n_kept,
        "dose": dose,
        "h2_baseline_sum": h2_baseline_sum,
    }

    return bursts, rec_row


def build_cache(
    mice: list[int] | None = None,
    limit: int | None = None,
    cache: Path = CACHE,
) -> None:
    """Fill ``cache/bursts.csv`` and ``cache/recordings.csv`` incrementally (see module doc)."""
    cache.mkdir(parents=True, exist_ok=True)

    bursts_path = cache / "bursts.csv"
    recs_path = cache / "recordings.csv"

    # Load existing recordings.csv to know what's already done
    existing: set[tuple[int, str]] | None = None
    if recs_path.exists():
        existing = set(
            zip(pd.read_csv(recs_path)["mouse"].astype(int),
                pd.read_csv(recs_path)["recording"])
        )

    all_members = fus.archive.recordings()

    if mice is not None:
        all_members = [m for m in all_members if m.mouse in mice]

    # Filter out already-processed recordings
    if existing is not None:
        all_members = [m for m in all_members if (m.mouse, m.recording) not in existing]

    if limit is not None:
        all_members = all_members[:limit]

    # Append one recording at a time, save after each
    for i, member in enumerate(all_members, 1):
        print(f"  [{i}/{len(all_members)}] {member.filename}", flush=True)
        bursts, rec_row = read_recording(member)

        if bursts_path.exists():
            existing_df = pd.read_csv(bursts_path)
            combined = pd.concat([existing_df, bursts], ignore_index=True)
            combined.to_csv(bursts_path, index=False)
        else:
            bursts.to_csv(bursts_path, index=False)

        if recs_path.exists():
            existing_df = pd.read_csv(recs_path)
            combined = pd.concat([existing_df, pd.DataFrame([rec_row])], ignore_index=True)
            combined.to_csv(recs_path, index=False)
        else:
            pd.DataFrame([rec_row]).to_csv(recs_path, index=False)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mice", type=int, nargs="*")
    ap.add_argument("--limit", type=int, help="read at most this many new recordings")
    a = ap.parse_args(argv)
    build_cache(a.mice, a.limit)


if __name__ == "__main__":
    main()
