# Next email to Dr. Todd

> **Draft, 2026-11-01. Not sent.** Andrew's wording, with the results section filled in. The
> email that actually went out on 10/20, and Dr. Todd's 10/23 reply, are transcribed in
> `docs/lab/README.md`. The 10/28 call notes are in `notes/meeting_10-28.md`.

## Attach

1. `results/histology/slices/pooled_gfp_fraction.png`: Σ GFP+ neurons ÷ Σ NeuN neurons in
   each target's ROI, summed over slices, against dose. All 5 mice, 30 targets.
2. `results/histology/slices/reference_normalised.png`: Mice 2–5. Left panel: the same
   measure. Right panel: every target divided by its own mouse's shared target (120 bursts,
   goal 0.75).

## Draft

Andrew's final wording (2026-11-01).

> **Subject:** Mice 3–5 results + Mouse 6 files missing
>
> Hi Dr. Todd,
>
> Thanks for sending the new brains. I've run Mice 3, 4 and 5 through alignment and
> measurement. A few things to sort out first, then results.
>
> Data questions:
>
> 1. Mouse 6: The Mouse_Cntr_06 folder came through empty on my end in Dropbox, although the
>    Mouse 6 slide deck has 2 images. Could you re-upload the scans? I know you also mentioned
>    something about the slices being in the same file.
>
> 2. Mice 3, 4, and 5 are supposed to have 7, 5, and 3 images each correct? I ended up using
>    5, 4, and 2 respectively. With respect to this within mouse 3 inside
>    _Image_Cntr_3.11.2_03_ there's an extra tile folder (_Image_Cntr 3.10.2_01_, 1.1 GB) and
>    a .vsi that won't open. Is a 3.10.2_01 scan missing?
>
> 3. On the virus, the mice after mouse 6 all received AVcpp16 correct?
>
> Some results from Mice 1-5:
>
> I'm using 27 sections across the five mice (4, 12, 5, 4 and 2), with all others being a bit
> too messy to use. For the first result, it is the fraction of neurons that are GFP+ in each
> target's ROI: GFP+ neurons summed over all of a mouse's sections, divided by all neurons in
> the same ROIs. Because it's a fraction of the neurons actually present, it corrects for
> ROIs that run off the tissue edge, for some targets holding fewer neurons, and for
> differences in cell counts between mice.
>
> Pooled across all mice (figure 1): across the 30 targets, delivery tracks acoustic dose at
> r = 0.66 (0.65 without Mouse 1's no-FUS control). Within each mouse, r ranges from 0.64 to
> 0.85. See this attached as pooled_gfp_fraction.png.
>
> The second figure is normalized by the shared target we discussed. Mouse 1 is left out,
> because its shared treatment went into a site that also got two other sonications. The
> shared target received similar doses in Mice 2–5 (0.83–1.07), but its delivery ranges from
> 1.6% to 46% GFP+. The two lowest are both at T1, which reads low in all three new mice.
> Dividing by those near-zero values inflates everything else (up to 24x in Mouse 3), and the
> pooled correlation drops from 0.65 to 0.45. My read is that the new mice's sections may not
> reach the anterior targets, rather than delivery being low there. See this attached as
> reference_normalized.png.
>
> Let me know once mouse 6 is in the Dropbox and if those image counts sound correct to you.
>
> Thanks,
> Andrew

## Before sending

- [ ] Look at both figures. No person has seen them yet.
- [ ] The email names `reference_normalized.png`, but the file is `reference_normalised.png`.
      Rename one so they match.
- [ ] Optional: "drops from 0.65 (Mice 2–5)". Figure 2 leaves out Mouse 1, so its starting r
      is 0.65, not figure 1's 0.66.

## Where each number comes from

| Claim | Source |
|---|---|
| 27 sections: 4, 12, 5, 4, 2 | `data/section_orientation.csv` (of 32 exported; excluded 1, 0, 2, 1, 1) |
| 7, 5, 3 images in the decks | `data/histology/Mouse_0{3,4,5}_IF_Images.pptx`; also the number of scans received |
| r = 0.66 / 0.65, 30 targets | `results/histology/slices/README.md`, `gfp_fraction` rows |
| Per-mouse r 0.64–0.85 | `cumulative_dose_delivery.csv`, Pearson of `gfp_cells / cells` with dose per mouse (M1 0.75, M2 0.82, M3 0.85, M4 0.64, M5 0.70) |
| Shared target 1.6–46%, doses 0.83–1.07 | `results/histology/slices/reference_target.csv` |
| 24× and 0.65 → 0.45 | `results/histology/slices/reference_normalised.csv`, 20 non-shared targets, Mice 2–5 (normalised r on a log scale) |

## Left out on purpose (bring to a call if useful)

- **Summed GFP intensity** instead of cell counts: better within a mouse (median ρ 0.83 vs
  0.77), worse pooled, because exposure differs between slides. Saturation hasn't been
  checked.
- **Mouse 5's "extra data in QuPath":** every QuPath reader finds only 3 sections, and the
  deck has 3 images. Question 2 covers it implicitly.
- **CY5 = NeuN without DAPI** on the new scans: assumed, as for Mouse 2.
- **Mouse 17's `Target5Repeat.mat`** is really position 6. A rename on Dropbox would help.
