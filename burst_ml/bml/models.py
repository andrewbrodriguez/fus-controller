"""numpyro models M1-M5b, exactly as in burst_ml/README.md "Models". Do not change priors.

Each model is a ``Spec``: a pure-jax ``prob(th, d, a)`` giving delivery probability ``p`` per
target, the prior of every parameter in ``priors``, and whether it has an animal effect.

* ``th``: dict of one draw of every parameter in ``priors`` (scalars or arrays).
* ``d``: the design dict from ``bml.data`` (arrays may be numpy or jax).
* ``a``: (N,) animal effect per target, already indexed by animal (``a_j[d["animal"]]``).
  Models without an animal effect ignore it.

``make_model(spec)`` returns the numpyro model ``model(d, k=None)``:

1. ``th[name] = numpyro.sample(name, dist)`` for each item of ``spec.priors``.
2. If ``spec.hier``: ``sigma_a ~ HalfNormal(1)``, ``z ~ Normal(0, 1).expand([d["n_animals"]])``
   (sites named ``"sigma_a"``, ``"z"``), ``a = (sigma_a * z)[d["animal"]]``. Else ``a = 0``.
3. ``p = clip(spec.prob(th, d, a), 1e-6, 1 - 1e-6)``.
4. ``log_phi ~ Normal(log 10, 1)`` (site ``"log_phi"``), ``phi = exp(log_phi)``.
5. ``numpyro.sample("k", BetaBinomial(p * phi, (1 - p) * phi, total_count=d["n"]), obs=k)``.

``loglik(spec, th, d, a, k)`` returns the (N,) BetaBinomial log-probabilities of ``k`` (same
clipping), for scoring held-out data.

Burst sums run over ``d["mask"]`` only. Where a power or log of padded zeros would appear,
substitute a safe value with ``jnp.where`` *before* the operation (gradients of ``0 ** rho``
are NaN), then zero the padded terms.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpyro
import numpyro.distributions as dist
from numpyro.distributions import BetaBinomial, HalfNormal, Normal
from jax.scipy.special import expit

# float64 everywhere: BetaBinomial log-probabilities with n ~ 75,000 involve log-gamma terms
# near 1e5, where float32 rounding (~0.01-0.05 per target) would add noise to every score.
numpyro.enable_x64()

import jax.numpy as jnp  # noqa: E402

SCALE = 1e4   # the lab's dose = Σ h2 / 1e4

LOGIT_20 = -1.3862943611198906   # logit(0.2), prior centre for intercepts


def _intercept_slope():
    return {"beta0": dist.Normal(LOGIT_20, 1.5), "beta1": dist.Normal(0.0, 2.0)}


@dataclass(frozen=True)
class Spec:
    name: str
    prob: Callable
    priors: dict = field(default_factory=dict)
    hier: bool = True


# --- safe helpers for burst sums -----------------------------------------------------------

def _safe_power(base, exponent):
    """exp(exponent * log(base)), safe at base=0.

    Where base is 0 the result is 0 (for positive exponent).  The gradient of ``0 ** rho``
    is NaN when rho < 1, so we use the log-space form and mask padded terms afterwards.
    """
    safe_base = jnp.where(base > 0, base, 1.0)
    return jnp.where(base > 0, jnp.exp(exponent * jnp.log(safe_base)), 0.0)


# --- dose functions D(th, d) -> (N,) ------------------------------------------------------

def dose_m4(th, d):
    """``Σ_t exp(γ·(t/100 − 0.5))·h2_t / 1e4`` over masked bursts (γ = th["gamma"])."""
    gamma = th["gamma"]
    h2 = d["h2"]
    t = d["t"]
    mask = d["mask"]
    w = jnp.exp(gamma * (t / 100.0 - 0.5))
    result = (w * h2 * mask).sum(axis=1) / SCALE
    return result


def dose_m5a(th, d):
    """``(h_ref/1e4)·Σ_t (h2_t/h_ref)^ρ`` over masked bursts, ``ρ = exp(th["log_rho"])``."""
    log_rho = th["log_rho"]
    rho = jnp.exp(log_rho)
    h_ref = d["h_ref"]
    h2 = d["h2"]
    mask = d["mask"]
    x = jnp.where(mask, h2 / h_ref, 1.0)
    term = _safe_power(x, rho)
    # zero out padded (unmasked) terms
    term = jnp.where(mask, term, 0.0)
    result = (h_ref / SCALE) * term.sum(axis=1)
    return result


def dose_m5b(th, d):
    """As ``dose_m5a`` with each burst term also × ``exp(Z_t · th["beta_bands"])``
    (``Z`` is (N, T, 3), ``beta_bands`` is (3,))."""
    log_rho = th["log_rho"]
    rho = jnp.exp(log_rho)
    beta_bands = th["beta_bands"]
    h_ref = d["h_ref"]
    h2 = d["h2"]
    Z = d["Z"]
    mask = d["mask"]
    x = jnp.where(mask, h2 / h_ref, 1.0)
    term = _safe_power(x, rho)
    # band-weight: exp(Z @ beta_bands), zero where masked
    band_weight = jnp.exp(Z @ beta_bands)
    term = jnp.where(mask, term * band_weight, 0.0)
    result = (h_ref / SCALE) * term.sum(axis=1)
    return result


# --- probability functions prob(th, d, a) -> (N,) ----------------------------------------

def prob_m1(th, d, a):
    """``expit(alpha + a)``."""
    return expit(th["alpha"] + a)


def prob_m2(th, d, a):
    """``expit(beta0 + beta1·log dose)``; ignores ``a``."""
    return expit(th["beta0"] + th["beta1"] * jnp.log(d["dose"]))


def prob_m3a(th, d, a):
    """``expit(beta0 + a + beta1·log dose)``."""
    return expit(th["beta0"] + a + th["beta1"] * jnp.log(d["dose"]))


def prob_m3b(th, d, a):
    """``p0 + (pmax − p0)·expit(h·(log dose − log_ec50 − a))``, ``h = exp(th["log_h"])``."""
    h = jnp.exp(th["log_h"])
    return th["p0"] + (th["pmax"] - th["p0"]) * expit(
        h * (jnp.log(d["dose"]) - th["log_ec50"] - a)
    )


def prob_m4(th, d, a):
    """``expit(beta0 + a + beta1·log dose_m4(th, d))``."""
    dose_val = dose_m4(th, d)
    return expit(th["beta0"] + a + th["beta1"] * jnp.log(dose_val))


def prob_m5a(th, d, a):
    """``expit(beta0 + a + beta1·log dose_m5a(th, d))``."""
    dose_val = dose_m5a(th, d)
    return expit(th["beta0"] + a + th["beta1"] * jnp.log(dose_val))


def prob_m5b(th, d, a):
    """``expit(beta0 + a + beta1·log dose_m5b(th, d))``."""
    dose_val = dose_m5b(th, d)
    return expit(th["beta0"] + a + th["beta1"] * jnp.log(dose_val))


# --- the pre-registered models -----------------------------------------------------------

MODELS: dict[str, Spec] = {
    "m1": Spec("m1", prob_m1, {"alpha": dist.Normal(LOGIT_20, 1.5)}),
    "m2": Spec("m2", prob_m2, _intercept_slope(), hier=False),
    "m3a": Spec("m3a", prob_m3a, _intercept_slope()),
    "m3b": Spec("m3b", prob_m3b, {
        "p0": dist.Beta(1.5, 100.0),
        "pmax": dist.Beta(2.0, 2.0),
        "log_ec50": dist.Normal(0.0, 1.0),
        "log_h": dist.Normal(jnp.log(1.5), 0.5),
    }),
    "m4": Spec("m4", prob_m4, {**_intercept_slope(), "gamma": dist.Normal(0.0, 1.0)}),
    "m5a": Spec("m5a", prob_m5a, {**_intercept_slope(), "log_rho": dist.Normal(0.0, 0.5)}),
    "m5b": Spec("m5b", prob_m5b, {**_intercept_slope(), "log_rho": dist.Normal(0.0, 0.5),
                                   "beta_bands": dist.Normal(0.0, 0.25).expand([3])}),
}


def make_model(spec: Spec):
    """The numpyro model for ``spec`` (module doc, steps 1-5)."""

    def model(d, k=None):
        # 1. Sample each prior
        th = {}
        for name, prior in spec.priors.items():
            th[name] = numpyro.sample(name, prior)

        # 2. Animal effect
        if spec.hier:
            sigma_a = numpyro.sample("sigma_a", HalfNormal(1.0))
            z = numpyro.sample("z", Normal(0, 1).expand([d["n_animals"]]))
            a = (sigma_a * z)[d["animal"]]
        else:
            a = 0.0

        # 3. Compute probability, clip
        p = jnp.clip(spec.prob(th, d, a), 1e-6, 1 - 1e-6)

        # 4. Overdispersion parameter
        log_phi = numpyro.sample("log_phi", Normal(jnp.log(10), 1))
        phi = jnp.exp(log_phi)

        # 5. Likelihood
        numpyro.sample("k", BetaBinomial(p * phi, (1 - p) * phi, total_count=d["n"]), obs=k)

    return model


def loglik(spec: Spec, th, d, a, k):
    """(N,) BetaBinomial log-probabilities of ``k`` under one parameter draw ``th``
    (``th`` must include ``log_phi``)."""
    p = jnp.clip(spec.prob(th, d, a), 1e-6, 1 - 1e-6)
    phi = jnp.exp(th["log_phi"])
    return BetaBinomial(p * phi, (1 - p) * phi, total_count=d["n"]).log_prob(k)