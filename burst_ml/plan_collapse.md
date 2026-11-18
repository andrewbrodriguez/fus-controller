# Pre-registered (2026-11-18, before any delivery fit): within-burst microbubble collapse

Motivation: Tsitsos et al. 2025 (bioRxiv 10.1101/2025.01.17.633644) find delivery follows how fast
bubbles collapse/replenish, not frequency-domain dose within a group.

Features per sonicating burst, from `niscope.wf` (raw time series), computed label-free:
- `c2`: log ratio of 2nd-harmonic energy in the second half vs first half of the burst window
  (bandpass at 2*f0 +/- 10%, squared analytic envelope). Negative = collapse within the burst.
- `r2`: log ratio of this burst's 2nd-harmonic energy to the previous burst's (replenishment).
Per target: mean `c2`, mean `r2` over sonicating bursts (nan-safe).

Test (fixed now): within-mouse LOMO RMSE on logit(GFP+/NeuN), `ld` unpenalised, +mean c2, +mean r2
(two separate models and both together = 3 candidates); permutation null shuffles the candidate
within mouse and takes the best of the 3. Bar: gain > 0 and permutation p < 0.05. Run once; report all.
Labels: results/histology/slices/cumulative_dose_delivery.csv as of 2026-11-18.

## Result (run once, 2026-11-18; scripts/exploratory_collapse_test.py)
Held-out RMSE: dose 0.921; +c2 0.926 (p 0.43); +r2 0.898 (gain 0.023, p 0.15); +both 0.900 (p 0.16), p = best of 3.
Bar (gain > 0 and p < 0.05) not met. No evidence that within-record collapse or burst-to-burst replenishment adds to dose.
Caveat: the record is one 30 ms window per burst, so c2 is a coarse proxy for the paper's within-pulse collapse.
