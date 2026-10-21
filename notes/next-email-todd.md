# Next email to Dr. Todd

## Draft (2026-10-20)

Attach:

1. `results/histology/slices/cumulative_pooled.png`: delivery summed over every slice vs dose
2. `results/acoustic_eda/dose_vs_exposure.png`: 2nd-harmonic dose vs applied acoustic exposure, all 145
   sonications (the more detailed `controller_headline.png` is there if he asks)

> **Subject:** Update: cumulative delivery over all slices, acoustic dose across 24 mice
>
> Hi Dr. Todd,
>
> Two quick results and a few questions.
>
> **1. Delivery summed over every slice (figure 1).** As you suggested, I now sum each
> target's delivery over all slices of the animal, rather than normalising per slice. Across
> the 12 targets in Mice 1 and 2, cumulative 2nd-harmonic dose tracks delivery: r = 0.73 for
> GFP+ area and r = 0.68 for GFP+ neurons (0.62 and 0.57 without Mouse 1's no-FUS control).
> To put both mice on one axis, each animal's totals are divided by its mean target. It's
> only two animals, so I'm treating it as a consistent pattern rather than a fitted model.
>
> **2. Acoustic dose across all 24 mice (figure 2).** I pulled all four sessions from Dropbox
> and reprocessed every recording. All 138 target doses for mice 2–24 match the spreadsheet
> to four decimals. Figure 2 plots each sonication's 2nd-harmonic AUC against the acoustic
> exposure applied (peak pressure summed over bursts). They track each other (r = 0.77 on
> log axes), but loosely: at the same exposure, the harmonic dose varies about 2.7×
> (10th–90th percentile). The applied energy alone doesn't pin down the bubble response,
> which is the case for measuring it.
>
> Two controller observations: it settles at about 0.92× its goal at every setpoint, and it
> hits the 0.125 V ceiling more often in later sessions (2, 1, 4 and 8 of ~36). Mouse 6
> Target 2 hit the ceiling and the interlock, and ended at about 0.09× its goal, so it is
> close to an unsonicated target.
>
> Small thing: Mouse 17's `Target5Repeat.mat` is position 6 (it matches the spreadsheet's
> position 6 exactly, and your reduced file is named Target6). It may be worth renaming
> on Dropbox.
>
> **Questions**
> 1. **Normalising across mice:** have you come up with any good ways to normalise delivery
>    across animals? Right now I divide each animal's totals by its mean target, which
>    removes anything that scales a whole brain (slice count, staining, exposure) but also
>    removes real between-animal differences. One option I'm considering is the reference
>    target every mouse has (120 bursts at goal 0.75), as an internal standard. I'm curious
>    what you'd prefer.
> 2. **Sectioning:** for Mouse 2 and later mice, what is the section thickness, how far
>    apart are the mounted sections (consecutive or every Nth), and do Slide 1 → 4 and the
>    three sections on each slide run dorsal → ventral? Roughly how deep does the series
>    start, and was the protocol the same for every mouse? Your target mask puts the focal
>    centres about 3–3.7 mm below the dorsal surface, but matching Mouse 2's section areas to
>    the Allen atlas puts its 12 sections at roughly 1.7–2.7 mm, above most of the foci. Does
>    that fit how they were cut? This would tell us which slices fall inside the focal column,
>    and whether cuts are comparable across mice.
> 3. **Ceiling hits:** is the rise in voltage-ceiling hits in the later sessions expected
>    (animal size, coupling, transducer)?
> 4. **Next tissue:** when is the next batch of brains likely, how many, and which mice? Same
>    layout as Mouse 2? Mice 6 and 20 would be especially useful because of the
>    capped/failed targets above.
>
> Thanks,
> Andrew

Before sending: look at both figures (no person has viewed them yet).



---

## Background notes

Drafted 2026-10-19. Settled already, so don't ask again: the capsid (assume AAV9 for every
animal; it doesn't matter for our analysis), CY5 = NeuN (his Mouse 2 summary), no contrast
MRI, no confocal, and the exposures (we normalise per slice). The raw acoustic files are on
Dropbox and Andrew is pulling them. Head position per animal is **not** an ask yet, and we
assume nothing about it.

## 1. How the sections were cut

Lower priority for the dose–delivery analysis, but it is the one thing that turns the 3D
reconstruction from a sketch into a measurement, and it answers his own question about
whether cuts are consistent across mice. For Mouse 2, and whether it holds for later mice:

1. **Section thickness** (µm).
2. **Which sections were kept.** Consecutive sections, or every Nth? If every Nth, what is N,
   or equivalently, how far apart are neighbouring mounted sections?
3. **Order.** Do `Slide_01` → `Slide_04`, and the three sections on each slide (left → right?
   top → bottom?), run dorsal → ventral in cutting order?
4. **Starting depth.** Roughly how far below the dorsal surface the first kept section is,
   or which landmark it starts at.
5. **Was the same protocol used for every mouse?** If not, it will be recorded per animal.

Why we're asking: the brain-size estimate from the Allen atlas puts the 12 Mouse 2 sections
within ~1 mm of depth (~85–95 µm apart), which would cover only a third of the ~3 mm focal
column. If they are ~250 µm apart they span the whole column. The answer decides which
slices sit inside the focus and whether summing over slices means the same thing in every
mouse.

## 2. Next tissue

1. **When** is the next batch of brains likely, and **how many**?
2. **Which mice**, in what order? Anything known about each (damage, missing sections)?
3. Same file layout as Mouse 2 (one `.vsi` per slide, 3 sections each, notch on the left)?

The draft manuscript is due at the end of reading period, so whether we have 4 or 12
animals by December changes how the paper is framed.

## To include as updates (not questions)

- Two-animal cumulative result: 2f AUC vs delivery summed over slices, 12 targets,
  r = +0.73 (pixels, A) and +0.68 (GFP+ neurons, B); +0.62 / +0.57 without the no-FUS
  control. Two animals, so the p-values are optimistic. Source:
  `results/histology/slices/README.md`.
- One figure, after a human has looked at it: `results/histology/reconstruction/mouse02_topdown.png`
  or `results/histology/aligned/mouse02_grid.png`.
- Controller behaviour, all 145 recordings (mice 1–24): it settles at ~0.92× the harmonic
  goal in every session, only 34% of recordings average inside the band, and the interlock
  fired in 36/145. Voltage-ceiling hits climb across sessions (2, 1, 4, 8 of ~36); is that
  expected (older mice, thicker skulls, coupling)? Source: `results/acoustic_eda/README.md`.
- Housekeeping: Mouse 17's `Target5Repeat.mat` is really position 6 (it matches the sheet's
  P6 exactly). Worth a rename on Dropbox so nobody else trips on it.
