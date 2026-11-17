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

Run 2026-11-16 with the settings above (`bml.run`, ~20 min on CPU). Every table is generated:
[`results/README.md`](results/README.md) (by `bml/report.py`); CSVs beside it. All 49 fits
converged (max R̂ 1.008); the LOAO fits had no divergent transitions, and M5b's full-data fit
had 1 of 4,000.

### Verdict

**No. The per-burst record does not predict delivery better than the cumulative 2nd-harmonic
AUC.** None of the three burst models meets the pre-registered rule, and none comes close:

| vs M3a (the AUC in a hierarchical model) | Δelpd | SE | Δ/SE |
|---|---|---|---|
| M4, learned time weighting | +0.27 | 1.25 | +0.22 |
| M5a, per-burst power | −0.27 | 0.59 | −0.45 |
| M5b, power + broadband/sub/ultraharmonic | −0.54 | 0.91 | −0.60 |

The fitted burst parameters say the same thing. Each one's posterior covers the value that
turns the model back into the plain dose (`burst_params.png`):

- **γ (time weighting) = −0.78, 90% interval −1.74 to +0.09.** It leans towards early bursts
  counting for more, but 0 (flat) is inside the interval, and M4's small gain comes from one
  held-out mouse (Mouse 1, +0.95; the other five sum to −0.68).
- **ρ (per-burst power) = 1.28, 0.80 to 2.04.** This is the lab's question: does the same dose
  delivered as few intense bursts beat many weak ones? ρ = 1 (only the total matters) is well
  inside the interval. The data can't rule out a moderate advantage for intense bursts either,
  so the honest answer is "no detectable difference at n = 35", not "no difference".
- **Broadband, sub- and ultraharmonic weights ≈ 0** (all three intervals straddle 0 within
  ±0.5). Most bursts carry almost no energy in those bands: broadband sits at 1.09× the noise
  floor (median, max 1.20×), sub- and ultraharmonics at ~1.05×.

### What the dose model itself shows

- **Dose predicts delivery in an animal the model has never seen.** M3a beats the no-dose
  model by Δelpd +7.11 (SE 1.38, Δ/SE 5.2). Within a held-out mouse the predicted order of its
  targets matches the observed order with mean Spearman ρ +0.71.
- **The animal effect is real but not decisive here.** Dropping it (M2) costs 3.12 (SE 2.55),
  worse in 5 of 6 held-out mice. The animal SD is 0.61 on the logit scale (0.27–1.21).
- **No evidence that a plateau helps prediction.** Emax (M3b) vs logit–log (M3a): −0.25
  (SE 0.42). M3b puts the plateau at 61% GFP+ (41–86%) and EC50 at dose 1.41 (0.63–3.08), but
  the logit–log curve predicts held-out mice just as well. The flattening at Mouse 6's high doses
  that motivated assumption A15 doesn't need a plateau in the model; a log-dose slope covers it.
- **Calibration is good.** M3a's held-out 50% and 90% intervals contain 49% and 94% of the 35
  targets; mean PIT 0.49.

### Dose for a delivery level

For **20% GFP+ neurons in a new animal**, the required cumulative 2nd harmonic is **0.86 (90%
interval 0.22–3.24)** under M3a and **0.77 (0.19–2.80)** under M3b. The interval spans about
15-fold, and most of that is animal-to-animal spread. The dose reaching **90% of the Emax
plateau** is 7.2 (1.4–53), more than twice the largest dose any target received (3.33), so it
is extrapolation and shouldn't be used.

**What this means for the controller.** Retargeting the controller from an acoustic setpoint
to a delivery setpoint needs more than a better feature of the burst record. On these data the
burst record adds nothing to the cumulative 2nd harmonic. The limit is the animal: the same
dose gives very different delivery in different mice, and nothing in this analysis predicts
which mouse will respond strongly. The useful next question is whether an **animal-level**
acoustic quantity (skull transmission at f0, the drive voltage the controller needed, the
reference target's dose) explains the animal effect.

### Caveats

- n = 35 targets in 6 mice. With six folds, SEs are themselves rough; a Δ/SE of 2 is a
  minimum bar, not proof.
- Monte Carlo error is negligible next to these differences: re-drawing the animal effects
  moves one fold's score by ~0.004, and refitting with a new seed by ~0.03.
- Mouse 1 T1 combines three sonications (dose 2.44) and has leverage in every burst model.
  It is in all fits, as planned.
- Delivery is pipeline B (ΣGFP+ ÷ ΣNeuN), which no hand count has validated yet.
- The firing-order confound, capsid and section depth are untouched (see *Not tested here*).

### Exploratory: without Mouse 1 T1 (not pre-registered)

Mouse 1 T1 is the one target built from three sonications (dose 2.44), and it has leverage in
every burst model. Rerunning everything without it
([`results/exploratory_no_m1t1/`](results/exploratory_no_m1t1/README.md), 34 targets) leaves
the verdict unchanged. M4 +0.20 (SE 1.10), M5a +0.10 (0.82), M5b −0.09 (1.10) against M3a; no
burst model gets past |Δ/SE| = 0.2. ρ = 1.39 (0.83–2.29) and γ = −0.76 (−1.76 to +0.15) both
still cover their no-effect values. Dose still beats no-dose (+5.77, SE 1.78), and the dose for
20% GFP+ in a new animal is unchanged: 0.86 (0.22–3.37).

## Deviations

None from the analysis plan. Process and implementation details the plan didn't fix:

- **Sparky agents ran one at a time**, not 2–4 (see *Log*).
- **float64** throughout: BetaBinomial log-probabilities with n up to 74,550 lose ~0.01–0.05 per
  target in float32.
- **Held-out predictions use 1,000 of the 4,000 posterior draws** (evenly thinned) × 200 animal
  effects each. The plan fixed the 200, not the thinning.

## Log

What each sparky agent was asked to do, what came back, and what was kept or rejected.

Sparky is the local Qwen3.6-35B-A3B box. Every agent gets a one- or two-sentence brief that
points at a skeleton whose docstrings are the spec; the orchestrator (Claude) wrote the
skeletons, the plan and the cache tests, and checks every result.

**Concurrency broke the model.** Four agents started together all failed: turns ended after 2
output tokens (`finish_reason: length`), one agent produced repeated "the the …" text, and one
wrote a corrupt edit to `bml/cache.py` (restored from git). A direct one-line request during
that load also came back as repeated text; with no load the same request was answered
correctly. From then on agents ran **one at a time**, not the 2–4 the goal allowed.

| Agent | Task | Outcome |
|---|---|---|
| env | create `burst_ml/.venv` (jax 0.12, numpyro 0.22, `fus` editable) | **kept**; imports checked by hand |
| cache, data, models, invert | first attempt, four at once | **rejected**: degenerate output under concurrent load (above) |
| cache2 | `bml/cache.py`: per-burst bands for every recording | **kept** after review; removed an unused variable and a stray temp-file line, added progress output. `tests/test_cache.py` (orchestrator's) passes: Σh2/10⁴ + baseline = dose exactly, and h2 matches the earlier per-burst table burst for burst |
| — (orchestrator) | full cache build, all 24 mice | 145 recordings, 18,031 sonicating bursts, none truncated |
| data3 | `bml/data.py`: labels + padded per-burst design | **kept** after review. Its first draft selected bursts by recording name only (every mouse's `Target1` together) and `subset` dropped `animal`; it fixed both itself. Orchestrator made a partly missing target raise instead of silently dropping recordings, and replaced one test that assumed no cache. On the real data every target's burst sum matches its dose to 0.0015 (the left-out baseline bursts) |
| models2 | `bml/models.py`: the seven numpyro models, likelihood, `loglik` | **kept** after review; priors unchanged. Bug fixed by the orchestrator: the non-hierarchical model built its zero animal effect as `jnp.zeros(d["n"])` (an array *shaped by* the neuron counts). Orchestrator's reference tests pass: log-likelihood = scipy `betabinom`, M3b reduces to M3a at p0 = 0 / pmax = 1, burst doses by hand. Stopped after its first turn, once every test passed. Smoke fits on real data: ~1 s each, no divergences |
| loao | `bml/loao.py` (asked); it also wrote `run.py` and `invert.py` unasked, **in the main tree instead of its worktree** (it followed the absolute venv path in its brief) | **stopped.** `loao()`, `compare()`, `calibration()` kept. `fit()` used four numpyro APIs that don't exist, and its R̂ merged draws with vector dimensions; `predict()` looped over 200,000 draw × animal pairs in Python and sampled fractions with an extra dimension, mixing targets in every interval. Orchestrator rewrote `fit`/`predict` (vectorised with `jax.vmap`) and tightened the reference test to check the fraction shape. `run.py`: `fit` not imported, `--out`/`--quick`'s `n_a`/`--force` ignored, host devices set too late for parallel chains; rewritten on its outline. `invert.py`: reported log-doses as doses and called levels below `p0` unreachable; rewritten to the docstring, with orchestrator tests |
| — (orchestrator) | `bml/report.py` (generated `results/README.md`) and the full run | ~20 min; all fits converged |
| figures | `bml/figures.py`, the four aggregate figures | **kept** after review, developed against the quick-run outputs. Orchestrator fixes: `comparison.png`'s grey dots were the base model's own row, not the six per-mouse differences; `dose_response.png` drew its bands at the observed doses in data order (unsorted) instead of on a dose grid; band weights were labelled "band 0/1/2"; two unrequested per-model PNGs removed; units added to axes. **No person has looked at the figures yet** |
