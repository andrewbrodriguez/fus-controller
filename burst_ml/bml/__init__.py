"""Burst-level models of FUS dose -> AAV delivery. Plan and results: burst_ml/README.md."""

from pathlib import Path

HERE = Path(__file__).resolve().parents[1]      # burst_ml/
REPO = HERE.parent                              # research/
CACHE = HERE / "cache"
RESULTS = HERE / "results"
LABELS = REPO / "results/histology/slices/cumulative_dose_delivery.csv"
