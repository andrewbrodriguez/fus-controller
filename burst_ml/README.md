# Burst-level models of dose → delivery

**Question.** Does a model of the per-burst acoustic record predict AAV delivery better, out of
sample, than the lab's cumulative 2nd-harmonic AUC (`dose`) alone?

Goal and rules: [`goal.md`](goal.md). This README has two parts: the **plan**, fixed and
committed before any model was fitted to delivery, and the **results**, added afterwards.
The plan is not edited after the first fit. Deviations are listed under *Deviations*, with the
reason.

---

## Plan (pre-registered, committed before fitting)

### Data

- **Outcome.** For each target, `k` = ΣGFP+ neurons and `n` = ΣNeuN neurons over every kept
  section of its animal (`gfp_cells`, `cells` in
  `results/histology/slices/cumulative_dose_delivery.csv`). The `*_rel` columns are never
  used: they divide by the animal's other targets, which a controller can't know.
- **Evaluation set: the 35 sonicated targets of mice 1–6.** Mouse 1 T3 (no-FUS control, dose 0,
  1.0% GFP+) is left out of every fit and score. It informs only the prior on the Emax floor
  `p0`.
- **Bursts.** Per sonicating burst (the `n_baseline` pre-bubble bursts dropped), from
  `burst_ml/cache/` (built by `bml.cache` from `data/acoustic/US_Data.zip` via `fus.archive`):
  - `h2`: the lab's 2nd-harmonic area, normalised to the baseline bursts
    (`fus.extract`, `metrics.normalised().area_2nd`). Σ`h2`/10⁴ is the dose to within the
    baseline bursts' contribution (~0.0005).
  - `bb`, `sub`, `ultra`: broadband, subharmonic and ultraharmonic energy per burst, fold over
    the zero-drive baseline (`fus.features.burst_bands`).
- **Recording → target.** TargetN = position N, except Mouse 1: T1 = `Target1`, `Target2`,
  `Target3` concatenated in firing order; T2 = `Target2_Repeat`; T4–T6 = `Target4`–`Target6`.
  Burst time `t` restarts at 0 in each recording (each sonication gets its own bubbles).
- **Mice 7–24 (no labels) are used only for constants:** `h_ref`, the median sonicating `h2`
  over all 24 mice, and the mean and SD of each log band used to z-score it.

### Likelihood and shared structure

Every model uses the same likelihood, so their scores are comparable:

    k ~ BetaBinomial(n, p·φ, (1 − p)·φ),   log φ ~ Normal(log 10, 1)

The animal effect is non-centred: `a_j = σ_a·z_j`, `z_j ~ Normal(0, 1)`,
`σ_a ~ HalfNormal(1)`. `logit⁻¹` is the logistic function and `L = log D`.

### Models

| # | Name | Mean function | Extra priors |
|---|---|---|---|
| M1 | Animal only | `logit p = α + a_j` | `α ~ N(−1.4, 1.5)` |
| M2 | Pooled dose | `logit p = β0 + β1·log dose` (no animal effect) | `β0 ~ N(−1.4, 1.5)`, `β1 ~ N(0, 2)` |
| M3a | Hierarchical, logit–log dose | `logit p = β0 + a_j + β1·log dose` | as M2 |
| M3b | Hierarchical Emax | `p = p0 + (pmax − p0)·logit⁻¹(h·(log dose − log EC50 − a_j))` | `p0 ~ Beta(1.5, 100)`, `pmax ~ Beta(2, 2)`, `log EC50 ~ N(0, 1)`, `h ~ LogNormal(log 1.5, 0.5)` |
| M4 | Time-weighted bursts | as M3a with `D = Σ_t exp(γ·(t/100 − 0.5))·h2_t / 10⁴` | `γ ~ N(0, 1)` |
| M5a | Per-burst power | as M3a with `D = (h_ref/10⁴)·Σ_t (h2_t/h_ref)^ρ` | `log ρ ~ N(0, 0.5)` |
| M5b | Per-burst power + bands | as M5a, each burst also × `exp(β_bb·z_bb + β_sub·z_sub + β_ultra·z_ultra)` | `β_· ~ N(0, 0.25)` |

`z_·` is the log band energy, z-scored with the all-mice constants. M4 with `γ = 0` and M5a
with `ρ = 1` give `D` = the dose, so each burst model nests M3a. M5a's `ρ` is the lab's
question in one number: `ρ > 1` means intense bursts count for more than their share of
the dose (few bursts at a high goal beat many at a low goal), and `ρ < 1` the reverse.

In M3b the animal effect shifts EC50 (how much dose an animal needs), and the floor and
plateau are shared. M3a is the special case `p0 = 0`, `pmax = 1`.

### Fitting and scoring

- **Inference:** NUTS (numpyro), 4 chains × 1,000 warm-up + 1,000 draws, `target_accept 0.95`,
  fixed seeds. Report R̂ and divergences for every fit. A fit with R̂ > 1.01 or more than 1%
  divergent draws is flagged in the results, not silently dropped.
- **Leave one animal out (LOAO), 6 folds.** Fit on five mice; predict the sixth with a fresh
  animal effect drawn from the fitted population (`a_new ~ N(0, σ_a)`, 200 draws per posterior
  draw). Nothing about the held-out animal enters the fit.
- **Primary score:** `elpd_LOAO = Σ_j log p(y_j | y_−j)`, the log predictive density of each
  held-out animal's targets *jointly* (they share `a_new`), by Monte Carlo over posterior draws
  × `a_new` draws.
- **Comparison:** `Δelpd` against M3a, with SE = √6 · SD of the six per-animal differences.
- **Decision rule, fixed now:** a burst model (M4, M5a, M5b) **beats the AUC** only if
  `Δelpd vs M3a > 0` **and** `Δelpd / SE > 2`. Anything else is reported as "no evidence that
  the burst record adds to the AUC". The same rule decides M3b vs M3a (saturation) and M3a vs
  M1 (whether dose predicts at all).
- **Secondary:** per-target marginal log predictive density; PIT values and 50% / 90%
  interval coverage of the held-out GFP+ fraction; Spearman ρ between predictive median and
  observed fraction, within held-out animals.
- **Dose for a delivery level (full data, no hold-out):** from M3b and M3a fitted on all 35
  targets, the dose that gives **20% GFP+ neurons** in a new animal (posterior predictive over
  `a_new`, median and 90% interval, plus the probability it's unreachable because
  `pmax < 0.2`), and from M3b the dose reaching **90% of the plateau**,
  `EC50·9^(1/h)` for a new animal.

### Not tested here

- Firing order: from Mouse 2 on, position = order, and the AAV goes in after the last target,
  so T1's opening is oldest when the AAV arrives. That is confounded with position and isn't
  modelled.
- Capsid (all six imaged mice are AAV9), section depth, and the choice of delivery measure
  (decision 001) are taken as given.

---

## Results

*Not yet run.*

## Deviations

*None yet.*

## Log

What each sparky agent was asked to do, what came back, and what was kept or rejected.

| Agent | Task | Outcome |
|---|---|---|
