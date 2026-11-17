# Goal: burst-level ML for FUS dose → delivery

Run ~8 hours as orchestrator. Use sparky agents (`sparky-agent` skill) for most of the legwork:
feature caching, scripts, tests, sweeps. Keep design decisions and verification yourself. All
work lives in `burst_ml/`.

## Question

Does a model of the per-burst acoustic record predict AAV delivery better, out of sample, than
the lab's cumulative 2nd-harmonic AUC (`dose`) alone?

## Data

- Labels: `results/histology/slices/cumulative_dose_delivery.csv`, mice 1–6, 36 targets.
  Outcome = `gfp_cells` / `cells` (ΣGFP+ ÷ ΣNeuN). Never use the `*_rel` columns: they need the
  animal's other targets, which is leakage.
- Recordings: read with `fus.archive` from `data/acoustic/US_Data.zip`. It is the only copy:
  never unpack it, extract from it, or write to `data/`. Per-burst bands come from `fus.extract`.
- Mapping: TargetN = position N, except Mouse 1: T1 = Target1 + Target2 + Target3 (in firing
  order), T2 = Target2_Repeat, T3 = no-FUS control (no recording, dose 0).
- Mice 7–24 have recordings but no labels: use them only unsupervised.

## Models

Write the model list and metric into `burst_ml/README.md` before fitting anything.

Baselines, all scored leave-one-animal-out (6 folds):

1. Animal mean only (no dose).
2. Pooled logit(p) ~ log dose, no animal effect.
3. Hierarchical: beta-binomial (or binomial + observation-level effect), random animal
   intercept, mean function logit-log-dose and Emax (p0 pinned near the control, ~1%).

Burst level (aggregated-output regression): the target's prediction is g(Σ_t w(t)·f(x_t)) with
the same likelihood and animal intercept. Cumulative AUC is the case f = identity on the 2nd
harmonic, w = 1. Try:

4. a learned temporal weighting w(t) (early vs late bursts);
5. a small monotone per-burst nonlinearity f on the 2nd harmonic, then with wideband and
   sub/ultraharmonic bands added. Heavily regularised.

Metric: LOAO log predictive density, with the held-out animal's intercept marginalised, plus
calibration of the held-out GFP+ fraction. A burst model counts as better only if it beats
model 3 out of sample.

## Rules

- New dependencies (numpyro/JAX or PyMC) go in `burst_ml/.venv` with `fus` installed editable,
  not the main `.venv`. `.venv/bin/python -m pytest tests -q` must still pass, and so must the
  tests in `burst_ml/tests/`.
- Cache per-burst features, not spectra, under `burst_ml/cache/`. Git-ignore the cache and the
  venv. Disk has ~15 GB free.
- Sparky: 2–4 agents at once, one- or two-sentence briefs, `--worktree` for repo work, never
  pointed at `data/`. Verify every result yourself.
- Branch `burst-ml`. Commit locally with no Claude attribution; never push.
- Figures come from scripts, pool all mice, and need a human to look at them. Say so.
- n = 36 in 6 animals: report intervals, don't tune on the test folds, and treat a null result as
  a result.

## Done when

`burst_ml/README.md` holds: the model list fixed before fitting; a LOAO table for models 1–5;
the dose that reaches a chosen delivery level (inverted Emax), with an interval for an unseen
animal; a plain verdict on whether the burst-level models beat the AUC; and a log of what sparky
did and what was rejected. If 8 hours pass first, write the README with what is finished and
what is next.
