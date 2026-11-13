<p align="center">
  <img src="docs/figures/banner.png" alt="FUS Controller" width="100%">
</p>

<h1 align="center">FUS Controller</h1>

<p align="center">
  <b>Machine learning based acoustic feedback controller for FUS-BBB mediated AAV delivery</b><br>
  <sub>Research for Credit&middot; Fall 2026 &middot; Todd Lab, Brigham and Women's Hospital / Harvard Medical School</sub>
</p>

---

## Overview

The blood-brain barrier (BBB) blocks most gene therapies — including adeno-associated viruses (AAVs) — from reaching the central nervous system. Focused ultrasound (FUS) with circulating microbubbles opens it non-invasively and transiently, but the therapeutic window is narrow: too little acoustic energy and nothing is delivered, too much and tissue is damaged.

The focused ultrasound lab already runs a closed-loop controller during sonication. It ramps drive voltage until the second-harmonic emission from cavitating microbubbles reaches a setpoint — the **harmonic goal** — an approach that traces back to Aryal et al. (2014). What that controller regulates, though, is an *acoustic* quantity. Nobody has yet shown how tightly it maps onto the quantity that actually matters: **how much AAV ends up in the brain**.

That gap is this project. The controller's prescribed dose (`N bursts × harmonic goal`) predicts the measured cumulative second-harmonic emission with r ≈ 0.76 — better than chance, far from deterministic. Skull geometry, vasculature, and microbubble kinetics fill in the rest. If acoustic emissions can be mapped to delivered AAV directly, the controller can be retargeted from an acoustic setpoint to a *delivery* setpoint.

This matters most for translation. Brute force is not an option in a patient: spending four or five minutes of sonication on every target is impractical, and longer exposure raises the risk of tissue damage. Knowing how much was delivered — rather than how much energy was applied — is what makes a shorter, safer exposure defensible.

## A real example: Mouse 1

<p align="center">
  <img src="docs/figures/mouse01-example.png" alt="Mouse 1 brain section with two circled targets. Target 1 was sonicated: its acoustic trace jumps about 40 dB when microbubbles arrive, and 72% of its circle is GFP-positive. Target 3 was a no-FUS control: nothing was recorded, and 2% of its circle is GFP-positive." width="100%">
</p>

One mouse, one AAV injection, six planned targets. At **target 1** the controller drove the transducer while microbubbles circulated. The blue trace is the second-harmonic emission it regulates, about 40 dB above baseline once the bubbles arrive. **Target 3**, 2.5 mm behind it in the same brain, was left unsonicated as a control. In this section, 72% of the target 1 circle expresses GFP from the delivered virus; the control reads 2% (66% vs 2% averaged over four sections).

That is the whole premise: where the acoustic signal says the BBB opened, the virus got in. The open question, and the reason this project exists, is ***how much***.

## Five animals, pooled

Five brains are imaged so far (Mice 1–5, 27 sections). Each target's delivery is summed over every section of its animal, then plotted against acoustic dose (cumulative second harmonic) for all 30 targets on one axis:

<p align="center">
  <img src="results/histology/slices/cumulative_pooled.png" alt="Acoustic dose against cumulative GFP delivery for 30 targets from five mice, three delivery measures side by side" width="100%">
</p>

| Delivery measure (30 targets, 5 mice) | Pearson r with dose |
|---|---|
| GFP+ fraction of neurons, ΣGFP+ ÷ ΣNeuN | **+0.66** |
| GFP+ area ÷ the animal's mean target | +0.61 |
| GFP+ cells ÷ the animal's mean target | +0.56 |

Within single sections the pattern holds too: the higher-dose targets carry more GFP in 17 of 17 sections with clear signal under the pixel measure (median ρ +0.60), and 16 of 17 under the cell count (+0.61). Which cross-animal measure to use is not settled yet, and the p-values treat targets as independent even though they come from five animals. Full tables and caveats: [`results/histology/slices/README.md`](results/histology/slices/README.md).

## From pixels to cells

Coverage counts GFP+ *pixels*, most of which are processes rather than cell bodies. A second pipeline counts *neurons*. It segments every NeuN-stained cell with StarDist, scores each by its mean anti-GFP intensity relative to its slice's background, and tags it GFP+ above 4.05× background.

<p align="center">
  <img src="docs/figures/gfp-tagging-example.png" alt="Three panels of the same 400 micrometre square at the edge of a Mouse 2 target. Left: NeuN in grey and the anti-GFP stain in green, bright on one side where the plume is. Middle: the same with every segmented neuron outlined, yellow where tagged GFP-positive and magenta where negative. Right: the outlines alone on black, yellow cells clustered on the plume side." width="100%">
</p>

*A 400 µm window on the edge of a Mouse 2 target, where the plume fades out. Left: NeuN (grey) and the anti-GFP stain (green). Middle and right: each segmented neuron outlined **yellow** (GFP+) or **magenta** (GFP−), over the stains and alone. 269 cells, 61 GFP+.*

In a 3 × 3 mm square around that target, 16,723 neurons are segmented. 97–100% are tagged within 0.5 mm of the target centre, falling to under 1% beyond 1.25 mm. Normalising to each slice's own background keeps Mouse 1's unsonicated control at about 1% GFP+.

This runs as **pipeline B** (`python -m fus.cells`) on every target ROI, about 1–1.5 min per slice, alongside the pixel pipeline (A). The two agree closely (r = 0.97 across 162 ROIs). It hasn't been checked against a hand count yet. Method, numbers and caveats: [`docs/gfp-cell-tagging.md`](docs/gfp-cell-tagging.md).

## Approach

**1. Feature extraction from acoustic emissions.** Per-burst spectra are reduced to peak and integrated power in the second-harmonic (1.674 MHz) and wideband (1.7 MHz) bands, normalised to the pre-microbubble baseline, then accumulated over the sonication. **Cumulative second-harmonic AUC is the primary feature**, as the lab uses it. Secondary, once that is in place: subharmonic and ultraharmonic bands, the harmonic-to-broadband ratio as a stable-vs-inertial cavitation index, and the temporal shape of the emission trace rather than its sum alone.

**2. Ground truth from tissue.** Brains are sectioned, stained, and imaged on a slide scanner. The delivery measurement is GFP area coverage: the fraction of each target's 2 mm circle that is GFP+ on the anti-GFP stain, following the thresholding approach of Owusu-Yaw et al. (2024). Each section is oriented by hand: its front, its notch side, and optionally each target's circle are set in napari, because the target pattern is too symmetric to place reliably from the GFP alone. A cell-level count (segment neurons, tag the GFP+ ones) is being piloted as an alternative; [decision 001](docs/decisions/001-delivery-endpoint.md) sets how to choose between them.

**3. Mapping acoustics to delivery.** Acoustic features are joined to their tissue measurements per target and used to fit a model predicting delivery from emissions. Six targets per mouse at different exposures means the design carries within-animal contrasts, which the model should exploit rather than ignore.

## Study design

24 mice × 6 targets, exposure duration crossed against controller setpoint:

| Factor | Levels |
|---|---|
| Harmonic goal (setpoint) | 0.55, 0.65, 0.75, 0.85, 0.95, 1.05 |
| Number of bursts | 57 – 255 |
| Brain region | 6 targets per animal |
| Virus | AAV9-GFP, 7.6e12 vg/kg, mice 1–12; AAV.CPP16, mice 13–24 |

The acoustic recordings were made with an **837 kHz** carrier (second harmonic at 1.674 MHz), sampled at 5 MHz and stored as per-burst magnitude spectra. The virus is delivered by tail vein immediately after sonication. Mice 1–12 got AAV9 and 13–24 AAV.CPP16 (lab, 2026-11-09). The capsid is not modelled as a factor: the comparisons are within-animal dose vs delivery.

## Project status

Current state, blockers and next steps: [`notes/current.md`](notes/current.md). Working
conventions for anyone (or any agent) picking this up: [`AGENTS.md`](AGENTS.md).

As of 2026-11-04, week ~10 of 13:
- **Acoustics:** feature extraction reproduces the lab's numbers; all 24 mice (4 sessions) are processed.
- **Tissue:** 5 of 24 brains are imaged, oriented by hand, measured with both pipelines and joined to dose. Mouse 6's scans arrived 2026-11-11 (12 sections, not yet ingested).
- **Still open:** which cross-animal delivery measure to use, a hand count to validate pipeline B, and the section order and depth.

The limiting factor is still tissue.

## Repository layout

```
data/           Experiment summary sheet (tracked); raw acoustic + histology (git-ignored)
src/fus/        Analysis package — feature extraction, quantification, models
scripts/        Helper scripts, incl. QuPath Groovy
tests/          pytest suite
notebooks/      Per-animal ingest (orient, finetune, measure, review), full-resolution slice viewer
segmentation_alpha/  Cell segmentation pilot and its environments (StarDist CPU/GPU, Cellpose)
results/        Generated CSVs and aggregate (all-mice) figures
docs/           Method write-ups, decisions, lab materials (docs/lab/), papers (docs/literature/)
notes/          Meeting notes and working log (`current.md` = state of play)
ARCH.md         Project charter: scope, weekly plan, status
AGENTS.md       Conventions, commands and gotchas for contributors
```

Raw data is not in version control — a single `.mat` acoustic recording is roughly half a gigabyte, and one imaged mouse is several more. See [`data/README.md`](data/README.md) for the layout, the Dropbox source, and how to read the file formats.

## Getting started

```bash
python3.12 -m venv .venv && source .venv/bin/activate   # Python ≥ 3.10 required
pip install -e ".[dev]"
pip install "napari[all]"                               # optional: interactive viewer used in the notebook

# register the venv as a Jupyter kernel for notebooks/
python -m ipykernel install --user --name fus-research --display-name "Python 3.12 (fus-research .venv)"
```

In Jupyter or VS Code, pick the **Python 3.12 (fus-research .venv)** kernel.

Then sync the Dropbox `US_Data` and `IF_Data` folders into `data/acoustic/` and `data/histology/` as described in [`data/README.md`](data/README.md).

Summarise a single recording, or reduce a whole experiment day to one row per target:

```bash
python -m fus Mouse_Cntr_01_Target1
python -m fus Mouse_Cntr_01_Target1 --plot
python -m fus --mouse 1 2 3 4 5 6 --csv results/day1.csv --quiet
```

`src/fus/extract.py` is a documented port of the lab's `nt_ExtractHarmonicData.m`. It reproduces the `Cumulative 2nd Harmonic`, `Mean Voltage`, and `Cumulative Voltage` columns of `Mouse_Controller_Data.xlsx` exactly, verified against 15 targets from the 2026-06-11 session.

`src/fus/histology.py` measures GFP area coverage at each target on the whole-slide scans — export with QuPath, then measure per section:

```bash
python -m fus.histology export data/histology/Mouse_01/Image.vsi data/processed/histology/Mouse_01/ds4
python -m fus.histology measure data/processed/histology/Mouse_01/ds4/*.ome.tif --csv results/histology/mouse01_slots.csv
```

It is a first pass with provisional parameters; see [`docs/histology-pipeline.md`](docs/histology-pipeline.md) for the method, the Mouse_01 results, and the open questions.

For each new animal, run [`notebooks/ingest_new_histology.ipynb`](notebooks/ingest_new_histology.ipynb). It exports the sections, opens napari for the orientation clicks and optional finetuning, measures, and opens a review window. To inspect any section at full resolution, use [`notebooks/view_slice.ipynb`](notebooks/view_slice.ipynb). Pipeline B, the cell count, runs per section or per mouse:

```bash
.venv/bin/python -m fus.cells --mouse 2 --all          # ~1-1.5 min per section
.venv/bin/python scripts/slices_dose_delivery.py       # both pipelines vs dose, per slice and pooled
```

It needs the StarDist environment in `segmentation_alpha/` (setup in [`docs/gfp-cell-tagging.md`](docs/gfp-cell-tagging.md)). [`notebooks/cell_pipeline_one_slice.ipynb`](notebooks/cell_pipeline_one_slice.ipynb) shows both pipelines on one slice, and the pilot is described in [`segmentation_alpha/README.md`](segmentation_alpha/README.md).

The `.mat` files are MATLAB v5 — read them with `scipy.io.loadmat`, not `h5py`. Whole-slide `.vsi` scans open in [QuPath](https://qupath.github.io/); [ImageJ/Fiji](https://imagej.net/software/fiji/) works for tile-level work.

## Deliverable

A written report in the form of a draft journal manuscript, covering the feature extraction methodology, the experimental protocol, and the predictive accuracy of the baseline model. If the model proves accurate, it gets packaged for use in the lab's acquisition workflow.

See [`ARCH.md`](ARCH.md) for the full charter and week-by-week plan.

## People

**Andrew Rodriguez** — student investigator
**Dr. Nicholas Todd** — principal investigator

## Key references

Full index in [`docs/literature/README.md`](docs/literature/README.md).

- Owusu-Yaw et al., *Pharmaceutics* 2024 — FUS-mediated AAV9 delivery in a HD mouse model; the predecessor study and the source of the GFP quantification method.
- Todd et al., *Molecular Therapy* (in preparation) — FUS-mediated AAV gene therapy in the CNS delivery landscape.
- Aryal et al., *Adv Drug Deliv Rev* 2014 — BBB disruption review; origin of harmonic-emission-based control.
- Kofoed et al., *Mol Ther Methods Clin Dev* 2021 — AAV9 vs. engineered PHP capsids.
- Thevenot et al., *Hum Gene Ther* 2012 — early AAV delivery across the FUS-opened BBB.
