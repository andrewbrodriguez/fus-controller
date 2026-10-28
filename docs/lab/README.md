# Lab materials

Files from the Todd lab, kept as received. Agents can't read the images in them; use the
`.pdf.md` sidecars (`pdftext --sidecar`).

| File | What |
|---|---|
| `ControllerStudy_Plots.pptx` | Controller-study slides (git-ignored, like every `.pptx`) |
| `Mouse_02_Summary.pdf` | **Current** Mouse 2 summary, from N. Todd's 2026-09-29 4:33 pm email ("left side was notched … updated pdf attached"). PDF created 2026-09-29 14:49. One page: the six targets with bursts, harmonic goal and cumulative harmonic dose (T1 0.62, T2 1.02, T3 0.76, T4 1.07, T5 0.37, T6 1.60, matching our recordings), NeuN/GFP images, and the label **"FUS + AAV9-GFP @ 7.6e12 vg/kg"** |
| `gABI_template_100um.nii` | From N. Todd's 2026-09-18 2:25 pm email ("Mouse brain template and mask for 6 target locations, assuming 2x2x2mm diameter. One voxel = 0.1 x 0.1 x 0.1 mm."). The Allen (ABI) average template at 100 µm, converted from MINC; 114 × 132 × 80 voxels ordered (ML, AP, DV), RAS (+x right, +y anterior, +z dorsal). The header says 1.0 per voxel, but voxels are really 0.1 mm. It matches our Allen CCF annotation (`fus.reconstruct.load_atlas`) at brain-mask IoU 0.83 after flipping AP and DV |
| `Mask_AAVTargets.nii` | Same email and grid. Binary mask of the six targets as 2 mm-diameter sphere **surfaces** (hollow; 5,548 voxels, and the pole voxels come loose as separate pieces). Centres, at template voxel (ML, AP, DV): (72, 86, 38), (92, 66, 38), (72, 61, 38), (42, 86, 38), (22, 66, 38), (42, 61, 38). These are exactly the `posx.xyz` plan, ±1.5/±3.5 mm lateral with rows 2.5 mm and 0.5 mm apart, so they are T1–T6 in that order with +x = the animal's right. All six sit at one depth |
| `Mouse_02_Summary_superseded_2026-09-29.pdf` | Earlier version, created 2026-09-29 11:57 and attached to that morning's email, before the notch correction. Same targets and doses, without the virus and channel labels. Kept only for provenance |

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
  upper half of each focus, and T3/T6 (deepest) only its top. This is asked of N. Todd in
  `notes/next-email-todd.md`.

## Email chain (N. Todd → A. Rodriguez, re "Project update: mouse 1 results and a few questions")

Transcribed from a screenshot on 2026-10-19; the screenshot was then deleted.

- **2026-09-18.** "New data set in the dropbox folder - Mouse_02. Four slides, each with 3
  sections on them for 12 total images." He thought the notch was cut at the bottom
  **right** of the brain, pending confirmation with the student; summary attached on that
  assumption.
- **2026-09-29 11:59 am.** Summary attached (the superseded PDF above).
- **2026-09-29 4:33 pm.** "Update - left side was notched, which actually makes more sense
  with the delivery and harmonic doses for the different targets. Updated pdf attached. Also
  looks like there was some rotation in this mouse, with the right side targets being
  higher and the left side targets being lower."

## Email chain (A. Rodriguez ↔ N. Todd, re "Updates and questions")

Transcribed from the thread Andrew pasted on 2026-10-25. Meeting prep: `notes/meeting_10-28.md`.

- **2026-10-20 8:16 pm, Andrew → Nick.** Sent a shortened version of the draft in
  `notes/next-email-todd.md`, with two figures attached (`cumulative_pooled.png` and
  `dose_vs_exposure.png`). It covered the cumulative result (r = 0.73 area, 0.68 cells; 0.62 /
  0.57 without the control) and dose vs exposure (r = 0.77, ~2.7× spread), and asked about
  normalising across mice, sectioning, ceiling hits and the next mouse. It left out the Mouse 17
  rename, Mouse 6 T2, and the depth mismatch.
- **2026-10-23 5:51 pm, Nick → Andrew.**
  - **Tissue:** "Josh just finished sectioning and staining more brains - I believe 4 more,
    which I am going to image on Monday [10/26]."
  - **Metrics:** three candidate outcome metrics:
    - **Number of GFP+ cells:** "add up across sections within fixed area ROIs and Total
      #cells is comparable across targets and mice."
    - **Area above threshold:** fixed ROIs summed over sections, giving total area or volume.
      The threshold can't be shared across mice: "put an ROI in a non-GFP region, get the
      average signal for that mouse and set the threshold at something like 1 or 2 stds above
      the mean."
    - **Mean GFP signal:** normalise by the same non-GFP ROI's mean. Combining across sections
      is "tricky … We can discuss this one more."
    - He also linked a ChatGPT answer (`chatgpt.com/s/t_6adbc746e6f88191b670ac2f47a896e5`,
      "Normalize GFP intensity"). Not transcribed: the page doesn't render without logging in.
  - **Sectioning:** "Section thickness is 30 microns. For spacing, I believe they are 12 * 30
    um = 360 microns apart. But let me confirm with Josh. Although, for mouse 2 where we you
    have 12 sections, I think you got two pairs of 6, some are adjacent."
  - **Ceiling:** "I wouldn't have predicted it would be related to time in the session - if
    anything I would think intensity would go down as microbubble dose builds up with
    multiple sonication targets." He read "later sessions" as later *within* a session, but the
    email meant the four surgery days.
  - **Next step:** "Let's find a time to chat over zoom next week for more on #1. I can do most
    times M - W."
- **2026-10-24 11:28 am, Andrew → Nick.** Proposed Wednesday 10/28 at 4:30. No confirmation in
  the pasted thread.
