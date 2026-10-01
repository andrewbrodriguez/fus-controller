# ES91r — Independent Study / Research

**PI:** Dr. Nicholas Todd
**Primary Mentor (day-to-day):** Dr. Bernie Owusu-Yaw
**Semester:** Fall 2026 (13-week independent study)
**Schedule:** Tue 1:00–5:00 PM, Wed 10:00 AM–12:00 PM, Thu 10:00 AM–2:00 PM
**Location:** Brigham and Women's Hospital / HMS (lab safety training required)

## Project Title
Machine Learning Based Acoustic Feedback Controller for FUS-BBB Mediated AAV Delivery

## Status — 2026-09-30 evening (week ~5 of 13)

Moment-in-time detail, including blockers and next steps: [`notes/current.md`](notes/current.md).

| Phase | Plan | Actual |
|---|---|---|
| 1 — Preparation & data (wk 1–4) | Import, filter, feature extraction | **Done and ahead.** `nt_ExtractHarmonicData.m` ported to Python and reproduces the lab's numbers exactly; harmonic and wideband bands, cumulative dose, interlock detection |
| 2 — Experimental (wk 5–8) | In vivo work, sectioning, staining, microscopy | **Partly, out of order.** Two brains measured (16 sections), oriented and finetuned by hand in napari; hand placement passes 4 of 5 checks against Mouse 1's answer key. Two delivery measures run on every target ROI: pixel coverage (A) and GFP+ neuron fraction (B, StarDist on NeuN), which agree at r = 0.98. B is not yet hand-validated |
| 3 — Integration & ML (wk 9–12) | Align ground truth to acoustics, fit a model | **Started early.** Dose and delivery are joined per slice for both mice and both measures: within-slice ρ is positive in 11/11 slices with clear signal (A) and 10/11 (B). No model yet — two animals |
| 4 — Report (wk 13) | Draft manuscript | Not started |

**Ahead of plan:** the computational side. Feature extraction, the histology pipeline with
hand orientation and finetuning, the per-slice dose analysis, a full-resolution viewer, the
cell pipeline (B) on a GPU StarDist environment, 33 tests, and method docs are all in place.

**Behind plan, and the limiting factor:** tissue. Two of 24 brains are imaged, and one of four
acoustic sessions is synced. The subharmonic/ultraharmonic work in phase 1 was deliberately
deferred — the lab asked to lead with cumulative 2nd-harmonic AUC.

**Settled with the PI on 2026-09-18:** recording→target mapping, the notch for section
orientation, the anti-GFP stain as the delivery measure, and a ~2 × 2 × 3 mm focal spot.

**Open:** capsid assignment per mouse; which marker CY5 shows on Mouse 2; review verdicts
(lost to a since-fixed bug); a hand count for the cell pilot (decision 001); tissue
throughput this term; the remaining acoustic sessions.

## Description
An independent study bridging computational modeling with wet-lab biological measurement, applying machine learning and signal processing to a translational medicine problem. The blood-brain barrier (BBB) is a major obstacle for delivering novel gene therapies (AAVs) to the central nervous system. Focused Ultrasound BBB (FUS-BBB) opening is a promising, non-invasive technique that uses targeted acoustic waves and circulating microbubbles to temporarily increase BBB permeability. This project designs a dynamic feedback controller that processes acoustic emissions in real time to predict and control the quantity of AAV delivered to the brain.

## Detailed Project Description
During this 13-week independent study, the project focuses on the first two phases of designing the acoustic feedback controller, heavily utilizing signal processing, machine learning, and quantitative biology:

1. **Data Processing:** Process raw, time-series acoustic emission data from previous FUS experiments (some already collected, more may be collected summer 2026). Use Python to design feature extraction pipelines (Fourier transform shows preliminary promise on the second harmonic) to isolate acoustic signatures correlated with BBB opening.

2. **Experimental Validation:** Assist with in vivo FUS experiments in mouse models. Perform tissue sectioning, staining, and fluorescence microscopy to quantitatively measure the ground truth concentration of AAVs delivered to target tissue.

3. **ML Integration:** Apply machine learning models to map extracted acoustic features (input) to physical AAV delivery outcomes (output).

## Weekly Schedule

### Phase 1: Preparation & Data (Weeks 1–4)

| Week | Focus |
|------|-------|
| 1 | Complete lab safety training (BWH/HMS), review relevant literature on FUS-BBB and acoustic cavitation, finalize computational environment setup |
| 2 | Familiarize with existing raw acoustic emission datasets; begin writing initial Python scripts for data import and basic filtering |
| 3 | Develop and test feature extraction pipelines (e.g., Fourier transforms) on existing acoustic data |
| 4 | Refine feature extraction algorithms to isolate harmonic, subharmonic, and broadband emission signatures |

### Phase 2: Experimental (Weeks 5–8)

| Week | Focus |
|------|-------|
| 5 | Begin assisting with in vivo FUS experiments; observe and learn the focused ultrasound protocols |
| 6 | Continue assisting with in vivo experiments; begin training on wet lab techniques for post-experiment tissue handling |
| 7 | Perform brain sectioning and biological staining on tissue from recent FUS experiments to prepare for AAV quantification |
| 8 | Utilize fluorescence microscopy and imaging software to quantitatively measure physical AAV delivery in stained tissue samples |

### Phase 3: Integration & ML (Weeks 9–12)

| Week | Focus |
|------|-------|
| 9 | Align ground truth AAV measurement data (from Week 8) with corresponding acoustic emission data logs |
| 10 | Begin building and training the baseline machine learning model using the newly structured dataset (mapping acoustic features to AAV delivery) |
| 11 | Test and evaluate model accuracy; perform iterative adjustments to the selected features or model parameters |
| 12 | Finalize data analysis; synthesize findings into visualizations (graphs, performance metrics) for the final report |

### Phase 4: Final Report (Week 13)

| Week | Focus |
|------|-------|
| 13 | Draft, review, and finalize the summative written report for submission by the end of Reading Period |

## Final Deliverable
- **Detailed written report** formatted as a draft journal manuscript, detailing:
  - The feature extraction methodology
  - The experimental protocol
  - The preliminary predictive accuracy of the baseline ML model
- **Model deployment** (dependent on need and accuracy): Make the model easily accessible and deployable to Dr. Todd's group

## Key Concepts
- Blood-brain barrier (BBB) and permeability
- Adeno-Associated Viruses (AAVs) for gene therapy delivery
- Focused Ultrasound BBB (FUS-BBB) mediated opening
- Acoustic cavitation and microbubble dynamics
- Acoustic emission monitoring (piezoelectric detection)
- Signal processing: Fourier transforms, harmonic/subharmonic/broadband analysis
- Feature extraction from time-series data
- Machine learning for biological prediction (acoustic features → AAV delivery outcomes)
- Wet lab techniques: tissue sectioning, biological staining, fluorescence microscopy
- In vivo FUS protocols in mouse models
- Translational medicine: bridging computational modeling with biological measurement

## Resources
- Existing acoustic emission datasets from previous FUS experiments (Dr. Todd's lab)
- May require additional data collection during summer 2026 (Dr. Todd's group)
- Computational environment: Python
- Lab: BWH/HMS — requires lab safety training

## Tools & Software
- Python (feature extraction, ML model development)
- Fourier transform libraries (for spectral analysis)
- Fluorescence microscopy and imaging software (for tissue quantification)
- Standard lab equipment for tissue sectioning and staining

## Grading / Evaluation
<!-- Not specified in the proposal — ES91r is typically graded on the basis of the final report and mentor evaluation -->

## Agent Notes

Working conventions, commands, and gotchas live in [`AGENTS.md`](AGENTS.md); current state in
[`notes/current.md`](notes/current.md).

- Research, not coursework: no lectures, no assignments. Deliverable is a draft journal
  manuscript plus, if it holds up, code the lab can run.
- Wet-lab work happens at BWH/HMS (Tue/Wed/Thu); computational work is done in this repo.
- The weekly plan above is the proposal's. Real sequencing has diverged — see the status table:
  computation ran ahead while tissue lagged.
- Tissue throughput sets the ceiling on everything. Judge scope against how many brains are
  realistically imaged, not against the plan's week numbers.
