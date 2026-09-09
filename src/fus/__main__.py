"""Command line entry point:  python -m fus <file.mat> [...]

Examples
--------
Summarise one recording::

    python -m fus data/acoustic/20260611/Mouse_Cntr_01_Target1.mat

Show the MATLAB figures plus the window-placement diagnostic::

    python -m fus data/acoustic/20260611/Mouse_Cntr_01_Target1.mat --plot

Reduce a whole experiment day to one CSV row per target::

    python -m fus data/acoustic/20260611/*_Target*.mat --csv results/day1.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from .extract import extract


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m fus",
        description="Extract per-burst cavitation metrics from FUS acoustic "
        "emission recordings (Python port of nt_ExtractHarmonicData.m).",
    )
    parser.add_argument("files", nargs="+", type=Path, help=".mat recordings")
    parser.add_argument(
        "--plot", action="store_true", help="show figures for each recording"
    )
    parser.add_argument(
        "--save-plots",
        type=Path,
        metavar="DIR",
        help="write figures to DIR as PNG instead of showing them",
    )
    parser.add_argument(
        "--csv", type=Path, metavar="PATH", help="write one summary row per file"
    )
    parser.add_argument(
        "--baseline",
        type=int,
        metavar="N",
        help="override the pre-microbubble burst count (default: from file)",
    )
    parser.add_argument(
        "--quiet", action="store_true", help="suppress the per-file summary"
    )
    args = parser.parse_args(argv)

    rows = []
    for path in args.files:
        if not path.exists():
            print(f"skipping {path}: not found", file=sys.stderr)
            continue

        ex = extract(path, n_baseline=args.baseline)
        rows.append(ex.to_dict())
        if not args.quiet:
            print(ex.summary())

        if args.plot or args.save_plots:
            from . import plots

            figures = {
                "metrics": plots.plot_cavitation_metrics(ex),
                "cumulative": plots.plot_cumulative_dose(ex),
                "spectrum": plots.plot_spectrum(ex),
            }
            if args.save_plots:
                args.save_plots.mkdir(parents=True, exist_ok=True)
                for name, fig in figures.items():
                    out = args.save_plots / f"{path.stem}_{name}.png"
                    fig.savefig(out, dpi=150)
                    print(f"  wrote {out}")
            else:
                import matplotlib.pyplot as plt

                plt.show()

    if args.csv and rows:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote {len(rows)} rows to {args.csv}")

    return 0 if rows else 1


if __name__ == "__main__":
    raise SystemExit(main())
