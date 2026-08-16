"""Regression coverage from the post-0d09740 deep audit."""

import builtins
import math

import mpmath as mp
import numpy as np
import pytest
from scipy import special

import elliptic


def _float(value):
    return float(np.asarray(value))


def test_m1_periods_do_not_hide_the_first_kind_pole():
    for phi in [2.0, math.pi, 4.0, 10.0]:
        F, E, Z = elliptic.elliptic12(phi, 1.0)
        turns = math.floor((abs(phi) + math.pi / 2.0) / math.pi)
        expected_E = (-1.0) ** turns * math.sin(abs(phi)) + 2.0 * turns
        assert math.isinf(_float(F))
        np.testing.assert_allclose(_float(E), expected_E, atol=1e-14)

        Fn, En, Zn = elliptic.elliptic12(-phi, 1.0)
        assert math.isinf(_float(Fn)) and _float(Fn) < 0
        np.testing.assert_allclose(_float(En), -expected_E, atol=1e-14)
        np.testing.assert_allclose(_float(Zn), -_float(Z), atol=1e-14)


def test_ellipticbd_preserves_small_parameter_limits():
    for m in [0.0, 1e-20, 1e-16, 1e-12, 1e-8]:
        B, D, S = elliptic.ellipticBD(m)
        np.testing.assert_allclose(_float(B), math.pi / 4.0, atol=2e-8)
        np.testing.assert_allclose(_float(D), math.pi / 4.0, atol=2e-8)
        np.testing.assert_allclose(_float(S), math.pi / 16.0, atol=2e-8)


def test_theta_exact_and_near_endpoint_parameters():
    v = 0.37
    assert _float(elliptic.theta(1, v, 0.0)) == 0.0
    assert _float(elliptic.theta(2, v, 0.0)) == 0.0
    assert _float(elliptic.theta(3, v, 0.0)) == 1.0
    assert _float(elliptic.theta(4, v, 0.0)) == 1.0

    mp.mp.dps = 50
    for m in [1e-20, 1e-15, np.nextafter(1.0, 0.0)]:
        mm = mp.mpf(float(m))
        q = mp.e ** (-mp.pi * mp.ellipk(1 - mm) / mp.ellipk(mm))
        for j in range(1, 5):
            expected = float(mp.jtheta(j, v, q))
            np.testing.assert_allclose(
                _float(elliptic.theta(j, v, m)),
                expected,
                rtol=2e-13,
                atol=2e-15,
            )


def test_elliptic3_near_pole_uses_full_precision_carlson_form():
    phi = np.array([1.2, 1.5, math.pi / 2.0])
    m = np.array([0.8, 0.95, 0.9])
    n = np.array([0.99, 0.999, 0.9999])
    s = np.sin(phi)
    c = np.cos(phi)
    d2 = 1.0 - m * s**2
    p = 1.0 - n * s**2
    expected = (
        s * special.elliprf(c**2, d2, 1.0)
        + n * s**3 / 3.0 * special.elliprj(c**2, d2, 1.0, p)
    )
    np.testing.assert_allclose(
        elliptic.elliptic3(phi, m, n), expected, rtol=2e-13, atol=2e-13
    )


def test_negative_phase_crossing_third_kind_pole_is_rejected():
    with pytest.raises(ValueError, match="Cauchy principal-value"):
        elliptic.elliptic3(-1.0, 0.5, 2.0)


def test_large_argument_ellipj_near_endpoint_parameter():
    mp.mp.dps = 70
    u = 1_000_000.123
    m = np.nextafter(1.0, 0.0)
    sn, cn, dn, _ = elliptic.ellipj(u, m)
    uu = mp.mpf(float(u))
    mm = mp.mpf(float(m))
    refs = [mp.ellipfun(name, uu, mm) for name in ("sn", "cn", "dn")]
    for got, expected in zip((sn, cn, dn), refs):
        np.testing.assert_allclose(_float(got), float(expected), rtol=2e-8, atol=2e-12)


def test_weierstrass_near_pole_is_finite_until_the_actual_lattice_point():
    values = elliptic.weierstrassP(np.array([0.0, 1e-11, 5e-11]), 1.0, 0.0, -1.0)
    assert np.isinf(values[0])
    assert np.all(np.isfinite(values[1:]))
    np.testing.assert_allclose(values[1:], np.array([1e22, 4e20]), rtol=2e-15)

    with pytest.raises(ValueError, match="real inputs only"):
        elliptic.weierstrassP(0.2 + 0.1j, 1.0, 0.0, -1.0)


def test_public_runtime_does_not_import_scipy(monkeypatch):
    real_import = builtins.__import__

    def reject_scipy(name, *args, **kwargs):
        if name == "scipy" or name.startswith("scipy."):
            raise AssertionError(f"unexpected SciPy runtime import: {name}")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", reject_scipy)

    elliptic.theta(1, 0.2, 0.5)
    elliptic.theta_prime(1, 0.2, 0.5)
    elliptic.jacobiThetaEta(0.2, 0.5)
    elliptic.nomeq(0.5)
    elliptic.inversenomeq(0.04)
    elliptic.inverselliptic2(0.5, 0.5)
    elliptic.elliptic12i(0.7 + 0.2j, 0.5)


def test_arclength_ellipse_is_vectorized_and_validates_axes():
    arcs = elliptic.arclength_ellipse(
        np.array([5.0, 10.0, 3.0]),
        np.array([10.0, 5.0, 3.0]),
    )
    np.testing.assert_allclose(
        arcs,
        np.array([48.44224110273839, 48.44224110273839, 6.0 * math.pi]),
        rtol=2e-14,
    )
    with pytest.raises(ValueError, match="strictly positive"):
        elliptic.arclength_ellipse(0.0, 1.0)
