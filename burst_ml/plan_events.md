# Pre-registered (2026-11-18, before any delivery fit): event statistics of the within-record envelope

Per burst (30 ms record, `niscope.wf`), envelope^2 of the band-passed signal (+/-10%), label-free:
- `k2`: log kurtosis of the 2f0 envelope (burstiness of bubble activity)
- `f2`: fraction of time the 2f0 envelope exceeds 2x its median
- `p2`: log(max / median) of the 2f0 envelope
- `ku`: log kurtosis of the 1.5f0 (ultraharmonic) envelope
- `fu`: fraction of time the 1.5f0 envelope exceeds 2x its median
Per target: mean over sonicating bursts. Test (fixed now): within-mouse LOMO RMSE, log dose unpenalised
+ each candidate singly (5 candidates); null shuffles candidates within mouse, best gain over the 5.
Bar: gain > 0 and permutation p < 0.05. Run once; report all. Labels as of 2026-11-18.

## Result (run once, 2026-11-18; scripts/exploratory_events_test.py)
Held-out RMSE vs dose alone 0.921: k2 0.942, f2 0.923, p2 0.926, ku 0.956, fu 0.915 (best gain 0.006, p 0.57, best of 5).
Bar not met. Event statistics of the 2f0 and 1.5f0 envelopes add nothing to dose.
