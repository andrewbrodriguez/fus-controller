"""Every target on one plot: acoustic dose against delivery, coloured by brain region.

x is the cumulative 2nd-harmonic AUC from the recordings; y is the fraction of neurons that
are GFP+ (pipeline B), summed over each target's kept slices (ΣGFP+ ÷ ΣNeuN), which needs no
per-animal rescale. Region is the planned target (T1/T4 anterior striatum, T2/T5 lateral
striatum / cortex, T3/T6 thalamus). Mouse 1's T1 (three sonications) and T3 (no FUS) are drawn
hollow.

Run after ``scripts/slices_dose_delivery.py``. Writes
``results/histology/region/dose_delivery_by_region.png``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import region_delivery as rd  # noqa: E402
import slices_dose_delivery as sdd  # noqa: E402

# Categorical slots 1-3 of the reference palette, in fixed order.
COLOR = dict(zip(rd.REGION_ORDER, ["#2a78d6", "#eb6834", "#1baf7a"]))


def main():
    d = rd.table()
    d["special"] = (d.mouse == 1) & d.target.isin([1, 3])
    sdd._style()
    fig, ax = plt.subplots(figsize=(8.4, 6.0))
    fig.subplots_adjust(left=0.09, right=0.97, top=0.86, bottom=0.13)
    for region in rd.REGION_ORDER:
        r = d[d.region == region]
        col = COLOR[region]
        for special, sub in r.groupby("special"):
            ax.scatter(sub.dose, sub.gfp_fraction * 100, s=64, zorder=3,
                       facecolor="none" if special else col, edgecolor=col if special else sdd.SURFACE,
                       linewidth=1.6 if special else 1.2,
                       label=None if special else f"{region} (T{'/T'.join(str(t) for t in sorted(r.target.unique()))}, n {len(r)})")
    ax.grid(color=sdd.GRID, linewidth=0.6, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(-0.06, d.dose.max() * 1.05)
    ax.set_ylim(-1.5, d.gfp_fraction.max() * 100 * 1.08)
    ax.set_xlabel("Acoustic dose (cumulative 2nd-harmonic AUC)")
    ax.set_ylabel("Delivery: % of neurons GFP+ (ΣGFP+ ÷ ΣNeuN)")
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    sdd._header(fig, "Dose vs delivery, every target, by brain region",
                f"{len(d)} targets from Mice {d.mouse.min()}-{d.mouse.max()}; region = planned target")
    sdd._footer(fig, "Hollow: Mouse 1 T1 (three sonications) and T3 (no-FUS control). "
                "Pipeline B, 1.74 mm ROIs, summed over each animal's kept slices.")
    out = rd.OUT / "dose_delivery_by_region.png"
    fig.savefig(out, dpi=170)
    print(f"wrote {out.relative_to(sdd.ROOT)}")
    print(d.groupby("region").agg(n=("dose", "size"), mean_dose=("dose", "mean"),
                                  mean_pct_gfp=("gfp_fraction", lambda x: 100 * x.mean())).round(2))


if __name__ == "__main__":
    main()
