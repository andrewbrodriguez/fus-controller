<p align="center">
  <img src="assets/banner.png" alt="Fluorescence micrograph of mouse brain tissue from the study, densely labelled in yellow and red with scattered patches of blue." width="100%">
</p>

<h1 align="center">FUS Controller</h1>

<p align="center">
  <b>Machine learning based acoustic feedback controller for FUS-BBB mediated AAV delivery</b><br>
  <sub>Research for Credit&middot; Fall 2026 &middot; Todd Lab, Brigham and Women's Hospital / Harvard Medical School</sub>
</p>

---

## Overview

The blood-brain barrier (BBB) blocks most gene therapies — including adeno-associated viruses (AAVs) — from reaching the central nervous system. Focused ultrasound (FUS) with circulating microbubbles opens it non-invasively and transiently, but the therapeutic window is narrow: too little acoustic energy and nothing is delivered, too much and tissue is damaged.

The lab already runs a closed-loop controller during sonication. It ramps drive voltage until the second-harmonic emission from cavitating microbubbles reaches a setpoint — the **harmonic goal** — an approach that traces back to Aryal et al. (2014). What that controller regulates, though, is an *acoustic* quantity. Nobody has yet shown how tightly it maps onto the quantity that actually matters: **how much AAV ends up in the brain**.

That gap is this project. The controller's prescribed dose (`N bursts × harmonic goal`) predicts the measured cumulative second-harmonic emission with r ≈ 0.76 — better than chance, far from deterministic. Skull geometry, vasculature, and microbubble kinetics fill in the rest. If acoustic emissions can be mapped to delivered AAV directly, the controller can be retargeted from an acoustic setpoint to a *delivery* setpoint.

## A real example: Mouse 1

<p align="center">
  <img src="docs/figures/mouse01-example.png" alt="Mouse 1 brain section with two circled targets. Target 1 was sonicated: its acoustic trace jumps about 40 dB when microbubbles arrive, and 90% of its circle is GFP-positive. Target 3 was a no-FUS control: nothing was recorded, and 1% of its circle is GFP-positive." width="100%">
</p>

One mouse, one AAV injection, six planned targets. At **target 1** the controller drove the transducer while microbubbles circulated. The blue trace is the second-harmonic emission it regulates, about 40 dB above baseline once the bubbles arrive. **Target 3**, 2.5 mm behind it in the same brain, was left unsonicated as a control. In this section, 90% of the target 1 circle expresses GFP from the delivered virus; the control reads 1% (82% vs 2% averaged over four sections).

That is the whole premise: where the acoustic signal says the BBB opened, the virus got in. The open question — and the reason this project exists — is ***how much***. Plotting all six Mouse 1 targets, acoustic dose against delivery:

<p align="center">
  <img src="docs/figures/mouse01-dose-delivery-scatter.png" alt="Scatter of cumulative second-harmonic dose against GFP coverage for the six Mouse 1 targets, with a least-squares line; r = 0.83" width="560">
</p>

The r = 0.83 comes almost entirely from the two ends: the unsonicated control, and target 1, which was sonicated three times. Among the four targets sonicated once (2, 4, 5 and 6), doses between 0.95 and 1.35 produced anywhere from 19% to 65% coverage, with no trend. That scatter is what a delivery-aware controller has to learn, and it takes more animals than one to learn it.

This is one animal, and several labels rest on assumptions still to be confirmed. The full method, the assumptions, and how each figure was made are in [`docs/mouse01-dose-delivery.md`](docs/mouse01-dose-delivery.md).

## Approach

**1. Feature extraction from acoustic emissions.** Per-burst spectra are reduced to peak and integrated power in the second-harmonic (1.674 MHz) and wideband (1.7 MHz) bands, normalised to the pre-microbubble baseline, then accumulated over the sonication. Beyond the lab's existing metrics: subharmonic and ultraharmonic bands, the harmonic-to-broadband ratio as a stable-vs-inertial cavitation index, and the temporal shape of the emission trace rather than its sum alone.

**2. Ground truth from tissue.** Brains are sectioned, stained, and imaged on a slide scanner. GFP area coverage in the FUS-targeted region — thresholded and expressed as percent of hemisphere, following the method in Owusu-Yaw et al. (2024) — is the delivery measurement the model learns against.

**3. Mapping acoustics to delivery.** Acoustic features are joined to their tissue measurements per target and used to fit a model predicting delivery from emissions. Six targets per mouse at different exposures means the design carries within-animal contrasts, which the model should exploit rather than ignore.

## Study design

24 mice × 6 targets, exposure duration crossed against controller setpoint, in two capsid cohorts:

| Factor | Levels |
|---|---|
| Harmonic goal (setpoint) | 0.55, 0.65, 0.75, 0.85, 0.95, 1.05 |
| Number of bursts | 57 – 255 |
| Brain region | 6 targets per animal |
| Capsid | AAV9 (12 mice), AAV.CPP16 (12 mice) |

The acoustic recordings were made with an **837 kHz** carrier (second harmonic at 1.674 MHz), sampled at 5 MHz and stored as per-burst magnitude spectra. The virus is delivered by tail vein immediately after sonication. The capsid split comes from the lab's controller-study slides (`docs/lab/ControllerStudy_Plots.pptx`); `Mouse_Controller_Data.xlsx` does not yet record which mouse received which, and any delivery model needs that as a factor.

## Repository layout

```
data/           Experiment summary sheet (tracked); raw acoustic + histology (git-ignored)
literature/     Reference papers — index tracked, PDFs git-ignored
reference/      Lab-provided material, incl. the MATLAB extraction script
src/fus/        Analysis package — feature extraction, quantification, models
scripts/        Helper scripts, incl. QuPath Groovy
tests/          pytest suite
notebooks/      Exploratory analysis
results/        Generated figures and model outputs
docs/           Proposal and lab presentations
notes/          Meeting notes and working log
assets/         Images used in documentation
ARCH.md         Project charter: scope, weekly plan, deliverables
```

Raw data is not in version control — a single `.mat` acoustic recording is roughly half a gigabyte, and one imaged mouse is several more. See [`data/README.md`](data/README.md) for the layout, the Dropbox source, and how to read the file formats.

## Getting started

```bash
python3.12 -m venv .venv && source .venv/bin/activate   # Python ≥ 3.10 required
pip install -r requirements.txt -e .
pip install "napari[all]"                               # optional: interactive viewer used in the notebook

# register the venv as a Jupyter kernel for notebooks/
python -m ipykernel install --user --name fus-research --display-name "Python 3.12 (fus-research .venv)"
```

In Jupyter or VS Code, pick the **Python 3.12 (fus-research .venv)** kernel.

Then sync the Dropbox `US_Data` and `IF_Data` folders into `data/acoustic/` and `data/histology/` as described in [`data/README.md`](data/README.md).

Summarise a single recording, or reduce a whole experiment day to one row per target:

```bash
python -m fus data/acoustic/20260611/Mouse_Cntr_01_Target1.mat
python -m fus data/acoustic/20260611/Mouse_Cntr_01_Target1.mat --plot
python -m fus data/acoustic/20260611/*_Target*.mat --csv results/day1.csv --quiet
```

`src/fus/extract.py` is a documented port of the lab's `nt_ExtractHarmonicData.m`. It reproduces the `Cumulative 2nd Harmonic`, `Mean Voltage`, and `Cumulative Voltage` columns of `Mouse_Controller_Data.xlsx` exactly, verified against 15 targets from the 2026-06-11 session.

`src/fus/histology.py` measures GFP area coverage at each target on the whole-slide scans — export with QuPath, then measure per section:

```bash
python -m fus.histology export data/histology/Mouse_01/Image.vsi data/processed/histology/Mouse_01/ds4
python -m fus.histology measure data/processed/histology/Mouse_01/ds4/*.ome.tif --csv results/histology/mouse01_slots.csv --figures results/histology/figures
```

It is a first pass with provisional parameters; see [`docs/histology-pipeline.md`](docs/histology-pipeline.md) for the method, the Mouse_01 results, and the open questions.

The `.mat` files are MATLAB v5 — read them with `scipy.io.loadmat`, not `h5py`. Whole-slide `.vsi` scans open in [QuPath](https://qupath.github.io/); [ImageJ/Fiji](https://imagej.net/software/fiji/) works for tile-level work.

## Deliverable

A written report in the form of a draft journal manuscript, covering the feature extraction methodology, the experimental protocol, and the predictive accuracy of the baseline model. If the model proves accurate, it gets packaged for use in the lab's acquisition workflow.

See [`ARCH.md`](ARCH.md) for the full charter and week-by-week plan.

## People

**Andrew Rodriguez** — student investigator
**Dr. Nicholas Todd** — principal investigator
**Dr. Bernie Owusu-Yaw** — primary mentor

## Key references

Full index in [`literature/README.md`](literature/README.md).

- Owusu-Yaw et al., *Pharmaceutics* 2024 — FUS-mediated AAV9 delivery in a HD mouse model; the predecessor study and the source of the GFP quantification method.
- Todd et al., *Molecular Therapy* (in preparation) — FUS-mediated AAV gene therapy in the CNS delivery landscape.
- Aryal et al., *Adv Drug Deliv Rev* 2014 — BBB disruption review; origin of harmonic-emission-based control.
- Kofoed et al., *Mol Ther Methods Clin Dev* 2021 — AAV9 vs. engineered PHP capsids.
- Thevenot et al., *Hum Gene Ther* 2012 — early AAV delivery across the FUS-opened BBB.
