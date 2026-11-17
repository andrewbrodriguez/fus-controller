"""Tests that all MODELS run under numpyro.infer.Predictive and have finite gradients."""

import numpy as np
import jax
import jax.numpy as jnp

import numpyro
import numpyro.infer
from numpyro.distributions import HalfNormal, Normal

from bml.models import MODELS, make_model, loglik


def _design(rng, N=4, T=6):
    h2 = rng.uniform(20, 200, (N, T))
    mask = np.ones((N, T), bool)
    mask[0, 4:] = False
    h2[~mask] = 0.0
    Z = rng.normal(size=(N, T, 3))
    Z[~mask] = 0.0
    t = np.tile(np.arange(T, dtype=float), (N, 1))
    dose = (h2 * mask).sum(1) / 1e4
    return {"animal": np.array([0, 0, 1, 1]), "n_animals": 2, "k": np.array([10, 200, 50, 900]),
            "n": np.array([1000, 1000, 800, 2000]), "dose": dose, "h2": h2, "t": t,
            "mask": mask, "Z": Z, "h_ref": 77.5}


def _numpy_params(spec, jax_key):
    """Return a dict of parameters drawn from each prior in ``spec``."""
    out = {}
    for name, dist_obj in spec.priors.items():
        params = dist_obj.sample(jax_key)
        if hasattr(params, 'numpy'):
            params = params.numpy()
        out[name] = params
    return out


def test_predictive_runs_all_models():
    """Every model in MODELS runs under numpyro.infer.Predictive without error."""
    rng = jax.random.PRNGKey(42)
    d = _design(np.random.default_rng(42))
    k = d["k"]

    for name, spec in MODELS.items():
        model_fn = make_model(spec)
        # Draw one set of params from the priors
        params = _numpy_params(spec, rng)
        # log_phi is sampled inside the model; we need to provide it too
        params["log_phi"] = np.log(10.0)
        if spec.hier:
            params["sigma_a"] = 0.5
            params["z"] = np.array([0.1, -0.1])

        predictive = numpyro.infer.Predictive(
            model_fn, params=params, num_samples=3, return_sites=["k"]
        )
        try:
            result = predictive(rng, d, k=d["k"])
        except Exception as e:
            raise AssertionError(f"Model {name} failed under Predictive: {e}") from e
        assert "k" in result, f"Model {name} missing 'k' in Predictive output"


def test_gradients_finite_padded():
    """Gradients of loglik are finite when some bursts are padded (masked out).

    Models M4, M5a, M5b compute dose from burst sums that involve powers and logs.
    With padded (zero) bursts, naive code would produce NaN gradients from ``0 ** rho``.
    """
    rng = jax.random.PRNGKey(99)
    d = _design(np.random.default_rng(99))

    for name, spec in MODELS.items():
        if not spec.hier:
            continue  # non-hierarchical models don't have animal effects to differentiate

        th = _numpy_params(spec, rng)
        # log_phi is required by loglik
        th["log_phi"] = np.log(12.0)
        a = np.array([0.3, 0.3, -0.2, -0.2])

        # loglik should be finite
        ll = loglik(spec, th, d, a, d["k"])
        assert np.all(np.isfinite(ll)), f"{name}: loglik not finite: {ll}"

        # gradient w.r.t. th params should be finite
        def _ll(th_dict):
            return loglik(spec, th_dict, d, a, d["k"]).sum()

        try:
            grads = jax.grad(_ll)(th)
        except Exception as e:
            raise AssertionError(f"{name}: grad failed: {e}") from e

        for param_name, g in grads.items():
            assert np.all(np.isfinite(g)), (
                f"{name}: grad of {param_name} not finite: "
                f"min={np.nanmin(g)}, max={np.nanmax(g)}"
            )


def test_loglik_finite_padded():
    """loglik returns finite values for all models, even with padded bursts."""
    rng = jax.random.PRNGKey(7)
    d = _design(np.random.default_rng(7))

    for name, spec in MODELS.items():
        th = _numpy_params(spec, rng)
        th["log_phi"] = np.log(12.0)
        a = np.array([0.3, 0.3, -0.2, -0.2])
        ll = loglik(spec, th, d, a, d["k"])
        assert np.all(np.isfinite(ll)), (
            f"{name}: loglik not finite: {ll}"
        )