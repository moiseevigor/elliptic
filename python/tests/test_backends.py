"""Backend smoke tests that make the Torch/JAX CI jobs exercise those paths."""

import numpy as np
import pytest

import elliptic


def _array(xp, values, *, complex_values=False):
    if xp is np:
        dtype = np.complex128 if complex_values else np.float64
        return np.asarray(values, dtype=dtype)
    if xp.__name__ == "torch":
        dtype = xp.complex128 if complex_values else xp.float64
        return xp.tensor(values, dtype=dtype)
    dtype = xp.complex128 if complex_values else xp.float64
    return xp.asarray(values, dtype=dtype)


def _numpy(value):
    if hasattr(value, "detach"):
        return value.detach().cpu().numpy()
    return np.asarray(value)


def test_core_functions_preserve_backend_and_values(xp):
    u = _array(xp, [0.2, 0.7, 1.1])
    m = _array(xp, [0.2, 0.5, 0.8])

    F, E, _ = elliptic.elliptic12(u, m)
    sn, cn, dn, _ = elliptic.ellipj(u, m)
    Pi = elliptic.elliptic3(u, m, 0.2)
    B, D, _ = elliptic.ellipticBD(m)
    Eu, Du, _ = elliptic.jacobiEDJ(u, m)

    for value in (F, E, sn, cn, dn, Pi, B, D, Eu, Du):
        assert _numpy(value).shape == (3,)
        assert np.all(np.isfinite(_numpy(value)))

    np.testing.assert_allclose(_numpy(sn) ** 2 + _numpy(cn) ** 2, 1.0, atol=2e-13)
    np.testing.assert_allclose(_numpy(Eu), _numpy(u) - _numpy(m) * _numpy(Du), atol=2e-12)


def test_auxiliary_functions_preserve_backend(xp):
    v = _array(xp, [0.1, 0.3])
    m = _array(xp, [0.4, 0.7])
    q = elliptic.nomeq(m)
    m_back = elliptic.inversenomeq(q)
    th = elliptic.theta(3, v, m)
    zeta = elliptic.weierstrassZeta(v, 1.0, 0.0, -1.0)
    sigma = elliptic.weierstrassSigma(v, 1.0, 0.0, -1.0)
    arc = elliptic.arclength_ellipse(
        _array(xp, [2.0, 3.0]),
        _array(xp, [3.0, 2.0]),
        0.1,
        _array(xp, [0.5, 0.7]),
    )

    np.testing.assert_allclose(_numpy(m_back), _numpy(m), atol=2e-13)
    for value in (th, zeta, sigma, arc):
        assert np.all(np.isfinite(_numpy(value)))


def test_complex_functions_preserve_backend(xp):
    u = _array(xp, [0.4 + 0.2j, 0.8 - 0.1j], complex_values=True)
    m = _array(xp, [0.3, 0.7])
    F, E, _ = elliptic.elliptic12i(u, m)
    sn, cn, _ = elliptic.ellipji(u, m)
    for value in (F, E, sn, cn):
        assert np.all(np.isfinite(_numpy(value)))


def test_jax_jit_core_paths():
    jax = pytest.importorskip("jax")
    jnp = pytest.importorskip("jax.numpy")
    jax.config.update("jax_enable_x64", True)
    u = jnp.asarray([0.2, 0.7], dtype=jnp.float64)

    outputs = [
        jax.jit(lambda x: elliptic.elliptic12(x, 0.5)[0])(u),
        jax.jit(lambda x: elliptic.elliptic3(x, 0.5, 0.2))(u),
        jax.jit(lambda x: elliptic.theta(3, x, 0.5))(u),
    ]
    # These larger graphs only need a tracing guard here; compiling the full
    # fixed-iteration inverse/theta expansions would make this smoke test a
    # poor CI citizen.
    jax.make_jaxpr(lambda x: elliptic.inverselliptic2(x, 0.5))(u)
    jax.make_jaxpr(lambda x: elliptic.weierstrassZeta(x, 1.0, 0.0, -1.0))(u)
    jax.make_jaxpr(lambda x: elliptic.arclength_ellipse(2.0, 3.0, 0.0, x))(u)
    assert all(np.all(np.isfinite(np.asarray(value))) for value in outputs)
