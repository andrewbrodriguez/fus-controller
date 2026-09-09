# Data

Raw data is **not** in version control. It lives in the lab Dropbox folder shared
by Dr. Todd. Only `Mouse_Controller_Data.xlsx` — the experiment summary sheet —
is tracked here, because it is small, text-like, and the key to everything else.

## Layout

```
data/
├── Mouse_Controller_Data.xlsx    tracked — experimental conditions for all targets
├── section_orientation.csv       tracked — per-section centre, front and notch-side clicks,
│                                 exclusions and review verdicts (notebooks/ingest_new_histology.ipynb)
├── roi_locations.csv             tracked — hand-finetuned T1–T6 shapes per section (export px),
│                                 from the finetune step of notebooks/ingest_new_histology.ipynb
├── acoustic/                     git-ignored
│   └── US_Data.zip               all 4 sessions as received from Dropbox (48.6 GB, stored, not
│                                 compressed). The ONLY copy: never unpack it; read it with
│                                 fus.archive. Inside, one folder per day (MMDDYY_FUS_AAV_Controller):
│       ├── Mouse_Cntr_XX_BL.mat          baseline, before microbubbles
│       ├── Mouse_Cntr_XX_TargetYY.mat    one per target (6 per mouse), ~0.4–0.6 GB
│       └── Mouse_Cntr_XX_*.jpg           per-run spectrogram screenshots
├── processed/histology/Mouse_NN/  git-ignored, all regenerable
│   ├── ds4/                      4x-downsampled section exports (1.3 µm/px), what pipeline A measures
│   ├── full/                     full-resolution exports (0.325 µm/px), ~7 GB each, made by
│   │                             notebooks/view_slice.ipynb; delete when done
│   └── cells/                    pipeline B per section: <section>_cells.csv (one row per cell),
│                                 <section>_T<k>_labels.tif (StarDist outlines)
└── histology/                    git-ignored
    ├── Mouse_01/
    │   ├── Image.vsi             Olympus whole-slide scan, 5 sections — open in QuPath
    │   └── _Image_/              tile stacks (.ets) backing the .vsi
    └── Mouse_02/
        ├── Control_02_Slide_0{1-4}.vsi   one scan per slide, 3 sections each
        ├── _Control_02_Slide_0{1-4}_/    tile stacks
    ├── Mouse_0{3,4,5,6}_IF_Images.pptx   lab IF slide decks (images only, no text):
    │                         7, 5 and 3 images for Mice 3-5 (2026-11-01); Mouse 6's is 12
    │                         slides, re-sent 2026-11-11 (slides 9-10 label positions with
    │                         bursts and harmonic goal)
    ├── Mouse_0{3,4,5}/           received 2026-11-01
    │   ├── Image_Cntr[_ ]<mouse>.<well>.<slide>_<scan>.vsi   one section each, CY5/TRITC/FITC
    │   │                     (no DAPI); export prefix w<well>sl<slide>_<scan>
    │   └── *_01_Overview.vsi     slide overviews, no sections (skipped on export)
    └── Mouse_06/                 received 2026-11-11 (Mouse_Cntr_06.zip; extraction CRC-checked)
        ├── Image_Cntr 6.<well>.<slide>_01.vsi   ONE FILE PER SLIDE, several sections each
        │                         (per docs/lab/IF_Images_Table.xlsx): 10.1 and
        │                         11.1 have 4 sections, 10.2 and 11.2 have 2, so 12 in all
        │                         (4/2/4/2 large tile stacks). No Overview scans. Channels
        │                         not yet checked. Export like Mouse 2, --prefix w<well>sl<slide>
        └── _Image_Cntr 6.<well>.<slide>_01_/    tile stacks
```

Imaging QC from the lab: `docs/lab/IF_Images_Table.xlsx` (one row per section, Mice 2-6).

Sections: Mouse 3 has 7, Mouse 4 has 5, Mouse 5 has 3. Mouse 3's `_Image_Cntr_3.11.2_03_/` also
holds a stray `_Image_Cntr 3.10.2_01_/` tile stack and a `.vsi` that won't open (asked the lab).

On disk: **all four acoustic sessions** (2026-10-19, `US_Data.zip`: 06-11 mice 1–6, 06-16
mice 7–12, 07-01 mice 13–18, 07-07 mice 19–24; 145 sonication recordings) and **six** imaged
mice (`Mouse_01`; `Mouse_02` from 2026-09-29; Mice 3–5 from 2026-11-01; Mouse 6 from 2026-11-11). The zip's `Mouse_Controller_Data.xlsx` and
`nt_ExtractHarmonicData.m` are byte-identical to the copies here. The previously unpacked
`20260611/` folder was CRC-identical to the zip's 06-11 folder and was deleted on 2026-10-19,
so there is one copy of every recording. Day-3 and day-4 baselines are named `_Baseline`
rather than `_BL`. Mouse 17's `Target5Repeat.mat` is position 6 (`scripts/sonication_eda.py`).

```python
from fus import archive
for m in archive.recordings():            # Member(day, mouse, recording, ...)
    with archive.local_path(m) as path:   # temporary copy (~0.3 s), deleted after
        ex = extract(path)
```

## Which recording went where

`TargetN` was fired at position N for every animal **except Mouse 1**, where the
operator deviated (confirmed by the lab on 2026-09-18 and by the acquisition
timestamps):

| Mouse 1 file | Fired at |
|---|---|
| `Target1`, `Target2`, `Target3` | position 1 — 420 bursts in total |
| `Target4`, `Target5`, `Target6` | positions 4, 5, 6 |
| `Target2_Repeat` | position 2 (7th and last sonication) |
| — | **position 3 was never sonicated** — a no-FUS control |

In the summary sheet, `Treatment #` is the firing order and `Brain Region` is correct.
For Mouse 1 only, the *values* in rows 3–6 are rotated by one row, so that sheet
understates the exposure at regions 4–6. Read doses from the recordings, not the sheet.

Sections are notched for orientation: on the animal's right for Mouse 1, on the left from
Mouse 2 on (lab, 2026-09-29, correcting an earlier report of right). Click the front and
notch of each section with `python -m fus.orientation`; see `docs/histology-pipeline.md#mouse_02`.

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
| Acquisition time | `data.niscope.tracktime` — MATLAB datenum; sorting files by it recovers the firing order of a session |
| Controller's own per-burst readings | `data.posx.data(end).HH` and `.WB` — the harmonic and wideband values the controller compared against its goal |
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
> used a 690 kHz transducer. Worth confirming with the lab which rig these came from before
> quoting parameters from that paper.
