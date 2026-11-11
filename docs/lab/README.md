# Lab materials

Files from the Todd lab, kept as received. Agents can't read the images in them; use the
`.pdf.md` sidecars (`pdftext --sidecar`). Correspondence with the lab is kept locally in
`correspondence.md` (git-ignored), not in the repo.

| File | What |
|---|---|
| `nt_ExtractHarmonicData.m` | Dr. Todd's MATLAB script: loads a `.mat` recording and extracts 2nd-harmonic and wideband metrics per burst. `src/fus/extract.py` is its port and reproduces its numbers |
| `ControllerStudy_Plots.pptx` | Controller-study slides (git-ignored, like every `.pptx`) |
| `AAV_By_TreatmentLocation.pptx` | From the lab, last modified 2026-11-06, received 2026-11-11 (git-ignored). One slide per target, Mice 1–6: bursts, mean MPa and cumulative harmonic dose (CHD), then sorted CHDs. Mouse 1 T1 is the three sonications summed (465 bursts, CHD 2.82); Mouse 1 T3 is 0 (never sonicated). Mouse 6 CHDs: T1 1.04, T2 0.47, T3 2.02, T4 0.92, T5 0.43, T6 3.33, matching our recordings (`results/acoustic_eda/recordings.csv`, 2026-11-11) |
| `IF_Images_Table.xlsx` | From the lab, last modified 2026-11-09, received 2026-11-11. Imaging QC, one row per section for Mice 2–6: mouse, well, slide, section, imaged, QC note (e.g. "Not ok - section too torn"). Mouse 5 well 11 slide 1 is "Not Stained?"; Mouse 6 is one file per slide |
| `Mouse_02_Summary.pdf` | **Current** Mouse 2 summary, received 2026-09-29 (it replaced an earlier version after the notch side was corrected to the left). PDF created 2026-09-29 14:49. One page: the six targets with bursts, harmonic goal and cumulative harmonic dose (T1 0.62, T2 1.02, T3 0.76, T4 1.07, T5 0.37, T6 1.60, matching our recordings), NeuN/GFP images, and the label **"FUS + AAV9-GFP @ 7.6e12 vg/kg"** |
| `gABI_template_100um.nii` | Received from the lab 2026-09-18, with the six target locations assuming 2 mm diameter and 0.1 mm voxels. The Allen (ABI) average template at 100 µm, converted from MINC; 114 × 132 × 80 voxels ordered (ML, AP, DV), RAS (+x right, +y anterior, +z dorsal). The header says 1.0 per voxel, but voxels are really 0.1 mm. It matches our Allen CCF annotation (`fus.reconstruct.load_atlas`) at brain-mask IoU 0.83 after flipping AP and DV |
| `Mask_AAVTargets.nii` | Received with the template, same grid. Binary mask of the six targets as 2 mm-diameter sphere **surfaces** (hollow; 5,548 voxels, and the pole voxels come loose as separate pieces). Centres, at template voxel (ML, AP, DV): (72, 86, 38), (92, 66, 38), (72, 61, 38), (42, 86, 38), (22, 66, 38), (42, 61, 38). These are exactly the `posx.xyz` plan, ±1.5/±3.5 mm lateral with rows 2.5 mm and 0.5 mm apart, so they are T1–T6 in that order with +x = the animal's right. All six sit at one depth |

## What the mask says about the targets (2026-10-20)

Read through the Allen CCF annotation (`data/processed/atlas/`), with each target as a solid 1 mm-radius sphere.
This assumes the mask's depth is the planned focal depth.

| Targets | Centre | Sphere contents | Centre depth below the dorsal brain surface |
|---|---|---|---|
| T1 / T4 (front, ±1.5 mm) | caudoputamen | 89% striatum | 3.0 mm |
| T2 / T5 (±3.5 mm) | caudoputamen, lateral | 54% striatum, 28% cortex, 16% fibre tracts | 2.9 mm |
| T3 / T6 (back, ±1.5 mm) | VPM thalamus | 92% thalamus | 3.7 mm |

Two checks against our assumptions:

- **Pattern vs brain centre.** At 1.7–2.7 mm depth, the brain's AP centroid falls at plan
  y = 1.85–2.0 (atlas mm). `histology.TEMPLATE_CENTRE_MM` assumes 1.6, from Mouse 1's GFP
  fits. That agrees within ~0.3 mm, which the hand finetuning of the ROIs absorbs.
- **Section depth.** `fus.reconstruct` matches Mouse 2's tissue areas to atlas depths of
  1.7–2.7 mm. The focal centres here are 2.9–3.7 mm deep, so if both hold, the sections cut the
  upper half of each focus, and T3/T6 (deepest) only its top. Asked of the lab 10/20 and
  again 11/01; still open.
