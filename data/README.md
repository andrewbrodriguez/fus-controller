# Data

Raw data is **not** in version control. It lives in the lab Dropbox folder shared
by Nick Todd. Only `Mouse_Controller_Data.xlsx` — the experiment summary sheet —
is tracked here, because it is small, text-like, and the key to everything else.

## Layout

```
data/
├── Mouse_Controller_Data.xlsx    tracked — experimental conditions for all targets
├── acoustic/                     git-ignored
│   └── <YYYYMMDD>/               one folder per experiment day, 6 mice each
│       ├── Mouse_Cntr_XX_BL.mat          baseline, before microbubbles
│       ├── Mouse_Cntr_XX_TargetYY.mat    one per target (6 per mouse), ~0.4–0.6 GB
│       └── Mouse_Cntr_XX_*.jpg           per-run spectrogram screenshots
└── histology/                    git-ignored
    └── Mouse_01/
        ├── Image.vsi             Olympus whole-slide scan — open in QuPath
        └── _Image_/              tile stacks (.ets) backing the .vsi
```

Currently synced locally: **one** acoustic day (`20260611`, mice 1–6) and **one**
imaged mouse (`Mouse_01`). The summary sheet covers **24 mice × 6 targets**, so
the remaining acoustic days still need to be pulled from Dropbox.

## Mouse_Controller_Data.xlsx

145 rows, one per treated target. Columns:

| Column | Meaning |
|---|---|
| `Mouse`, `Treatment #` | Animal ID and target index (1–6) |
| `Brain Region` | Target region code (1–6) |
| `N Bursts` | Sonication length — 57 to 255 bursts |
| `Harmonic Goal` | Controller setpoint — 0.55 to 1.05 |
| `Mean Voltage`, `Cumulative Voltage` | Drive voltage the controller settled on |
| `Mean 2nd Harmonic`, `Cumulative 2nd Harmonic` | Measured emission dose |
| `Mean MPa`, `Mean MI` | Derived in-situ pressure and mechanical index |

`N Bursts` × `Harmonic Goal` is a crossed dose design: exposure *duration* varied
against exposure *intensity*.

## Reading the .mat files

These are **MATLAB v5** files, not v7.3 — read them with `scipy.io.loadmat`
(`struct_as_record=False, squeeze_me=True`); `h5py` will not open them. Loading one
takes a few seconds and pulls the whole file into memory, so reduce each recording to
per-burst features once and cache the result rather than reloading raw spectra.

### What is actually stored

`data.niscope.ff(i).data0` is burst *i* — a **real-valued magnitude spectrum of 131072
points already transformed to the frequency domain**, not a raw time series. The FFT has
happened upstream in the acquisition software. `data.niscope.freqaxis` spans 0–2.5 MHz at
df = 19.07 Hz.

| Quantity | Where |
|---|---|
| Carrier frequency | `data.AWG.AWG1.Frequency` — **837 kHz**, so 2f = 1.674 MHz |
| Wideband monitor | 1.7 MHz — deliberately ~26 kHz off the harmonic, to catch broadband |
| PRF | `data.posx.curPRF` |
| Per-burst drive voltage | `data.posx.data(end).V`, shape (n_bursts, 2) |
| Raw time-domain | `data.posx.data(end).rawdata` — usually empty; retention is gated by the `rawdatasave` flags |

`data.posx` also carries the **full controller state**, which is the most useful and least
obvious part of these files:

- `Pcontrol`, `PcontrolKp` — proportional control and its gain
- `HarmGoal`, `curHarmGoal`, `SumHGoal` — the setpoint, including any mid-run change
- `WBthresh`, `WBaction`, `WBtrigger`, `reduceHarmGoalifWB`, `reduceMaxVifWB` — the
  broadband safety interlock and what it did
- `StartV1`, `MaxV1`, `HarmPosInc1`, `HarmNegInc1` — voltage ramp limits and step sizes
- `peakdata0..3`, `areadata0..3`, `Vrmsdata0..3` — the metrics the controller computed
  **online**, for four frequency bands

That last row matters: offline feature extraction can be validated against what the
controller itself recorded in real time.

`Reduced_*.mat` files (~12 KB) are the lab's already-reduced output from
[`reference/nt_ExtractHarmonicData.m`](../reference/nt_ExtractHarmonicData.m), useful as a
regression target when porting that script to Python.

> **Note on frequency:** this controller study runs at 837 kHz. Owusu-Yaw et al. (2024)
> used a 690 kHz transducer. Worth confirming with Nick which rig these came from before
> quoting parameters from that paper.
