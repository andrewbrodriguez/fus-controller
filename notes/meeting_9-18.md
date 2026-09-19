# Zoom with Nick — target mapping + next steps

spread sheet should be correct brain region all correct
treatment number is the order that we did them in, by mistake 3 sonications all R1
then 456 and then came back to 2 but never did 3

time of dose is a reason to not do this in a human, cant spend 4 or 5 minutes at each lcoation
and then also potential saftey, 

yes except for mouse 1, target number and location number always match except for mouse 1
treatment number always matches the brain region, except for mess up on mouse one

notch in the bottom right of mouse 1, bottom left going forward

pretty soon 2 more datasets to work with

start with 2nd harmonic AUC, for US data, maybe look at others
keep in mind extra US data, as secondary set of analysis






Follow-up to the 9/16 email. He answered Q1, 4, 5, 6; said "let's discuss" on Q2 (sheet vs deck)
and Q3 (does TargetN mean target N).

## Correct up front

- **I messed up Q2:** The sheet and the deck do *not* disagree. `Brain Region` reads
  1, 1, 1, 4, 5, 6, 2 for Mouse 1 — exactly his account.
- **Real issue:** Mouse 1 sheet rows 3–6 hold whole rows from the wrong runs, rotated by one
  (rows carry Target6, Target3, Target4, Target5; should be Target3, Target4, Target5, Target6).
- **Consequence:** sheet says regions 4/5/6 got 255@0.55, 105@0.95, 75@0.85. Actually
  105@0.95, 75@0.85, 255@0.55.
- **Mice 2–6 are clean** — names, firing order, treatment #, region all line up.

## Already answered — don't spend time

- Mouse 1: first three sonications all at location 1 (0.75/120, 0.85/60, 0.55/240); then 4, 5, 6
  correct; 7th back at location 2 (`Target2_Repeat`); location 3 never sonicated.
- Sections are horizontal. Flipped ones: notch in bottom right + target 3 unsonicated.
- Delivery measure = GFP **stain** (TRITC), not native.
- Focal spot ≈ 2 × 2 mm in x/y, ~3 mm in z (precise dims to follow).

## Tell him

- Every `.mat` stores its acquisition time (`niscope.tracktime`). Mouse 1's firing order from
  timestamps matches his account exactly — so any session's order can be reconstructed.
- Sheet-wide findings (all 24 mice, no raw files needed):
  - every mouse has the same reference target, 120 bursts @ 0.75 → emission varies 0.70–1.73
    (21% CV), voltage 0.067–0.116 V (14% CV) at identical settings
  - after removing the N × goal trend: only 21% of what's left is animal-level, 79% is
    target-level
  - `Mean MPa` = 4.05 × mean voltage exactly, `Mean MI` = MPa / √f — fixed conversions, no skull
    correction, no extra information
  - dose→emission slope varies by region: 88 (region 2) to 183 (region 6)
- Mouse 1 delivery: control 2%, sonicated targets 19–82%. r = 0.83 but carried by the control
  and target 1.

## Questions

**Mapping / bookkeeping**
1. How does the sheet get filled in, and from what source? Want corrected Mouse 1 rows from me?
2. Were the other three sessions transcribed the same way?
3. Any other animal with a re-target, repeat, or missed location like Mouse 1?
4. Proposal: I build `data/target_map.csv` (mouse, file, location, exposure) from timestamps +
   his notes; he confirms per animal; everything downstream joins on it. OK?

**Study design**
5. Which mice are AAV9 vs AAV.CPP16? Not in the sheet, and it has to be a factor in any model.
6. Same GFP reporter for both capsids?
7. Is there an untreated animal, or any other negative reference for setting the GFP threshold?
   Mouse 1's control only exists by accident.

**Tissue**
8. How many of the 24 brains are banked and sectionable this term? Who can section/stain/image,
   and at what rate?
9. Mouse 1 sections: order and spacing? Which are in the 3 mm focal column?
10. Are imaging exposure settings fixed across sessions (needed for intensity comparisons)?
11. The FITC-only signal at location 6 — autofluorescence, blood, or failed staining? Can show
    him the full-res crop.

**Data access**
12. Remaining three acoustic sessions on Dropbox? Any other imaged brains?

**Acoustics**
13. What are the dims of `posx.peakdata0` / `areadata0` (3 × 4 × 8 × bursts)? Which four bands
    do the trackers watch? (screenshot showed one at 2.5 × f₀)
14. Are raw time-domain waveforms saved anywhere, or only spectra?

**Direction**
15. ROI: use 1 mm radius (from the 2 × 2 mm spot) instead of my 0.75 mm placeholder?
16. What does he want first — acoustic EDA across all 24 mice, scaling the tissue pipeline, or
    the MATLAB port?
17. Report expectations: is packaging for lab use in scope this term?

## Leave the call with

- One authoritative recording → location mapping, and who owns it.
- Capsid assignment per mouse.
- A realistic number of brains imaged by December.
- Agreement on the ROI size and the delivery endpoint.
