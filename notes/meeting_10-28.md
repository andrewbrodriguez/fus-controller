# Zoom with Nick: normalising delivery across mice, sectioning, new brains

Wed 2026-10-28, 4:30 pm (Andrew proposed it on 10/24; check that Nick confirmed). This follows
the 10/20 email ("Updates and questions") and his 10/23 reply, both transcribed in
`docs/lab/README.md`. He asked for this call mainly to discuss **question 1, normalisation
across mice**.

## Notes from the call

(fill in during the meeting)


virus for the final batch AVcpp16

1) How many cumulative cells do we have across all valid versions of this ROI for this mouse
2) How much ROI intensity do we have across the ROIs as a % so 0-100
Then the score for delivery is something like 2/1 or the ROI score / the total number of cells

quick and dirty 
max delviery for each ROI, circle on it mean signal intesity per segment
mean gfp, at the end we have a number that is the mean gfp in the max slice
for each roi, thats like a top delivery metric

10 refer to the well, the 1 and 2 the 2 refer to the slides,

so mouse.well.slide_section

mouse 6 is a little bit diferent, no indiviudal sections



---

## What he already answered (10/23), so don't re-ask

- **Four more brains.** Josh has sectioned and stained them; Nick images them **Mon 10/26** and
  uploads them after that.
- **Sections are 30 µm thick, probably 12 × 30 = 360 µm apart** (he is confirming with Josh).
  For Mouse 2: "I think you got two pairs of 6, some are adjacent."
- **Outcome metrics, in his order of preference:**
  1. **Number of GFP+ cells:** fixed-area ROIs, summed over sections. "Total #cells is
     comparable across targets and mice."
  2. **GFP+ area above threshold:** fixed ROIs summed over sections, giving total area (or
     volume). The threshold has to be set **per mouse**: put an ROI in a non-GFP region, take
     its mean, and set the threshold at mean + 1–2 SD.
  3. **Mean GFP signal:** divide by the same non-GFP ROI's mean. Combining across sections is
     tricky, so "we can discuss this one more."
- He also sent a ChatGPT link (`chatgpt.com/s/t_6adbc746e6f88191b670ac2f47a896e5`, titled
  "Normalize GFP intensity"). **Not read**: the page doesn't render without logging in.
  Skim it before the call.
- **Ceiling hits.** "I wouldn't have predicted it would be related to time in the session - if
  anything I would think intensity would go down as microbubble dose builds up." **He read
  "later sessions" as later in a session.** See section 3.

---

## 1. Normalising across mice (the main topic)

### What our data say about his three metrics

Computed 10/25 from committed CSVs (`results/histology/slices/cumulative_dose_delivery.csv`),
12 targets from Mice 1–2. Pooled Pearson r with dose; values in brackets leave out the no-FUS
control:

| Delivery measure, summed over slices | Per-animal rescale? | r (no control) | Spearman ρ |
|---|---|---|---|
| GFP+ cells, raw total (his #1 as stated) | none | +0.09 (−0.10) | |
| GFP+ area, raw total (his #2) | none | +0.31 (+0.10) | |
| GFP+ area ÷ number of slices | none | +0.77 (+0.68) | |
| GFP+ cells ÷ number of slices | none | +0.39 (+0.20) | |
| **ΣGFP+ cells ÷ ΣNeuN cells (fraction of neurons)** | **none** | **+0.77 (+0.69)** | **+0.84 (+0.79)** |
| What I emailed: totals ÷ animal mean target, area | yes | +0.73 (+0.62) | +0.79 |
| What I emailed: totals ÷ animal mean target, cells | yes | +0.68 (+0.57) | +0.66 |

Within each mouse, the GFP+ fraction gives r = +0.75 (Mouse 1) and +0.82 (Mouse 2).

**Points to make:**

- **Raw totals aren't comparable while section counts differ.** Mouse 1 has 4 slices and
  Mouse 2 has 12, which is why raw totals barely correlate. Dividing by slice count, or better,
  multiplying by spacing (below), fixes this for area.
- **Raw cell counts also differ ~3× between mice** at the same ROI size: about 1,600 NeuN cells
  per ROI per slice in Mouse 1 and about 5,100 in Mouse 2. That points to NeuN staining,
  imaging, section thickness or segmentation, not biology. So "#GFP+ cells is comparable
  across mice" doesn't hold in our data yet. *Ask: was Mouse 1 cut, stained and imaged the
  same way as Mouse 2?*
- **Proposal: the GFP+ fraction of neurons, summed over slices (Σ GFP+ / Σ NeuN).** It
  normalises itself because the denominator comes from the same ROIs and slices. It needs no
  per-animal rescale, keeps real between-animal differences, and handles the thalamic targets
  (T3/T6), whose ROIs hold about half as many neurons in both mice. It does the best of
  anything above. Caveat: its GFP+ call still uses a per-slice threshold.
- **Area becomes volume once the spacing is known (Cavalieri):** GFP+ volume ≈ Σ area ×
  0.36 mm. That's in mm³ and comparable across mice, *provided every mouse is sampled through
  the same focal column*. Cells convert the same way (Σ count × 12 sections ≈ total count).

### His threshold recipe vs what we do

| | Nick | Pipeline A (pixels) | Pipeline B (cells) |
|---|---|---|---|
| Background | mean of a non-GFP ROI | smooth local background, σ 1 mm, plumes excluded | median anti-GFP over the slice's tissue |
| Scope | per mouse | **per slice** | **per slice** |
| Threshold | mean + 1–2 SD | residual > 5 noise SD | cell mean > 4.05 × background |

- **Argue for per slice rather than per mouse.** Exposures differ between slides (agreed 10/19),
  and Mouse 1's slice backgrounds ranged 435–1,927 counts. A per-mouse threshold would let
  bright slides dominate.
- **1–2 SD is loose.** If the background were Gaussian, 2–16% of pure background would pass.
  We currently hold Mouse 1's no-FUS control at 2% (A) and ~1% (B). Offer a sensitivity sweep
  (k = 1, 2, 3, 5) so we can choose k by where the control stays near zero.
- **Where the "non-GFP ROI" comes from.** From Mouse 2 on, every position is sonicated, so the
  only true control is Mouse 1 T3. Option: a fixed atlas region away from all six foci, chosen
  geometrically rather than by looking at GFP. The 3 mm crop pilot was under 1% GFP+ beyond
  1.25 mm from a target centre. *Ask whether he'd accept that.*
- **Mean intensity.** A already records `mean_residual` (background-subtracted intensity) per
  ROI, so we can show it, but intensity can't be compared across slides with different
  exposures. Ask what he wants it for.

### Reference target as an internal standard

- Every mouse has exactly one 120-burst, goal-0.75 target, at **P1 in 13 mice, P4 in 10 and P5
  in 1**. P1 and P4 are mirror images in anterior striatum, so it's nearly the same anatomy.
- **It fails for Mouse 1**: its reference shot went into position 1 together with two other
  sonications. For Mouse 2 it is T4.
- Its own dose varies 0.70–1.73 across mice (22% CV), so dividing by its delivery also divides
  by its acoustic noise. At best, use delivery per unit dose at the reference. Better as a
  **cross-mouse sanity check** than as the denominator.

**Leave with:** the cross-animal measure (my proposal: Σ GFP+ / Σ NeuN, with Cavalieri area
as the second measure), and whether thresholds are per slice or per mouse, at what k.

## 2. Sectioning, and what "two pairs of 6" means

**30 µm sections, every 12th kept → 360 µm apart.** This conflicts with what we've been
assuming:

- Atlas area matching (`results/histology/reconstruction/mouse02_spacing.csv`) put Mouse 2's 12
  sections at 1.7–2.7 mm depth, about 93 µm apart. `fus.reconstruct` defaults to 250 µm.
- **One series of 12 at 360 µm** would span ~4 mm, more than the whole focal column and not
  consistent with the area matching.
- **Two series of 6 at 360 µm**, one offset from the other so that sections pair up, would
  span ~1.8 mm, with each section 30 µm from its partner. This is my best reading of his
  remark.

**Candidate pairs**, from each section's neuron-count profile across the six ROIs (NeuN only,
not GFP; correlation of the profiles):

- slide01_s2 ↔ slide04_s2 (0.99)
- slide01_s3 ↔ slide04_s3 (0.98)
- slide02_s3 ↔ slide02_s4 (0.99)

This is weak evidence. Use it only to ask Josh, and don't set the order from it.

**Why it matters:**

- Summing over 12 sections that are really 6 pairs counts each depth twice. That's fine within
  an animal, but not across mice cut differently.
- Adjacent pairs are a **free test–retest**: two measurements 30 µm apart should agree.
- The order and spacing turn the 3D model from a sketch into a measurement.

**Ask:**

1. Which 6 sections form each series, and which are the adjacent pairs?
2. What is the cutting order on the slides (Slide 1 → 4, and left → right within a slide)?
3. How deep does the series start?
4. What was the protocol for Mouse 1 (4–5 sections on one slide)?
5. Will the four new brains be cut the same way? Ideally one series per brain, the same
   number of sections, starting at the same landmark, so totals are comparable without
   rescaling.

## 3. Ceiling hits: clarify, then show it's position, not time

Clarify first: by "sessions" I meant the **four surgery days** (6/11, 6/16, 7/1, 7/7; mice 1–6,
7–12, 13–18, 19–24), not time within a session.

Computed 10/25 from `results/acoustic_eda/recordings.csv`, mice 2–24:

| Position (= firing order) | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| Hit 0.125 V (of 23) | 0 | **6** | 1 | 0 | **5** | 3 |

- **11 of 15 ceiling hits are at the lateral targets P2/P5** (±3.5 mm), against 4 of 92 at the
  other positions (Fisher p = 0.001). The pattern by firing order isn't monotone, so it doesn't
  look like microbubble build-up over the session. Order and position are fully confounded in
  this design, though (P1 is always fired first), so this can't be proven.
- **The rise across days is mostly lateral targets:** 1, 1, 3, 6 lateral hits per day.
  Session 4's lateral hits are at goals 0.55, 0.85, 0.95 ×3 and 1.05, so mostly high goals.
  Median peak drive rises only slightly (0.098, 0.095, 0.100, 0.104 V).
- **Likely cause:** more bone and temporal muscle in the beam at the lateral sites.
- **Ask:** did anything change between days (animal age or weight, coupling or shaving,
  transducer, microbubble lot, operator)?
- **For the controller:** lateral targets need voltage headroom. Delivery is still measured
  against the *achieved* dose, so the dose–delivery analysis is unaffected.

## 4. Four new brains

- **Which mice?** Mice 6 (T2 hit the ceiling and the interlock, ending at ~0.09× goal) and 20
  (three ceiling hits) would test the low-dose end.
- **Same layout as Mouse 2?** One `.vsi` per slide, 3 sections each, notch on the left.
- **Plan:** `notebooks/ingest_new_histology.ipynb` per brain. The orientation clicks and ROI
  finetuning are a GUI step that Andrew does himself. Pipeline B takes ~1–1.5 min per section,
  so ~1 h for 48 sections.
- **Before ingesting:** settle the metric and threshold (section 1), so the new mice are
  measured once under the final rules and nothing is picked after looking at them.

## If there's time

- Mouse 17's `Target5Repeat.mat` is really position 6 (left out of the 10/20 email). Rename it
  on Dropbox?
- Mouse 2's sections sit at ~1.7–2.7 mm depth by area matching, above the focal centres the
  mask places at 2.9–3.7 mm (`docs/lab/README.md`). This may resolve once the spacing is known.
- The FITC-only signal at Mouse 1 T6 is still unexplained.

## Before the call

- [ ] Read the ChatGPT link.
- [ ] Optional: script the GFP+ fraction into `scripts/slices_dose_delivery.py` and the
      per-position ceiling table into `scripts/sonication_eda.py`, so the numbers above sit in
      committed READMEs rather than only here.
- [ ] Optional: a pipeline A threshold sweep at k = 1, 2, 3, 5 on Mice 1–2, ready for his
      "mean + 1–2 SD" suggestion.
- [ ] If the brains arrive Monday: export them, but hold off measuring until section 1 is
      settled.
