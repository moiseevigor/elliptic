"""Tests for elliptic12 (F, E, Z) and ellipj (sn, cn, dn, am)."""
import numpy as np
import pytest
from scipy.special import ellipkinc, ellipeinc, ellipj as sp_ellipj, ellipk, ellipe

import elliptic


class TestElliptic12:
    def test_scalar_vs_scipy(self):
        for phi, m in [(np.pi / 4, 0.5), (np.pi / 2, 0.7), (1.0, 0.3), (0.1, 0.99)]:
            F, E, _ = elliptic.elliptic12(phi, m)
            assert abs(F - ellipkinc(phi, m)) < 1e-12, f"F err at phi={phi}, m={m}"
            assert abs(E - ellipeinc(phi, m)) < 1e-12, f"E err at phi={phi}, m={m}"

    def test_array_vs_scipy(self):
        phi = np.linspace(0.01, np.pi / 2, 50)
        m = np.linspace(0.01, 0.99, 50)
        F, E, _ = elliptic.elliptic12(phi, m)
        np.testing.assert_allclose(F, ellipkinc(phi, m), atol=1e-12)
        np.testing.assert_allclose(E, ellipeinc(phi, m), atol=1e-12)

    def test_m0_special(self):
        u = np.array([0.3, 0.7, 1.2])
        F, E, Z = elliptic.elliptic12(u, 0.0)
        np.testing.assert_allclose(F, u, atol=1e-15)
        np.testing.assert_allclose(E, u, atol=1e-15)
        np.testing.assert_allclose(Z, 0.0, atol=1e-15)

    def test_complete_KE(self):
        m = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        F, E, _ = elliptic.elliptic12(np.full_like(m, np.pi / 2), m)
        np.testing.assert_allclose(F, ellipk(m), atol=1e-12)
        np.testing.assert_allclose(E, ellipe(m), atol=1e-12)

    def test_odd_symmetry(self):
        phi = np.array([0.3, 0.7, 1.1])
        m = 0.5
        Fp, Ep, Zp = elliptic.elliptic12(phi, m)
        Fm, Em, Zm = elliptic.elliptic12(-phi, m)
        np.testing.assert_allclose(Fp, -Fm, atol=1e-14)
        np.testing.assert_allclose(Ep, -Em, atol=1e-14)

    def test_z_jacobi_zeta(self):
        """Z = E - E(m)/K(m) * F."""
        phi, m = 0.8, 0.6
        F, E, Z = elliptic.elliptic12(phi, m)
        K_m = ellipk(m)
        E_m = ellipe(m)
        Z_ref = E - (E_m / K_m) * F
        assert abs(Z - Z_ref) < 1e-13

    def test_beyond_half_period(self):
        """F(phi+k*pi|m) = F(phi|m) + 2k*K(m), likewise E (issue #35)."""
        m = 0.4
        phi = np.array([-1.2, 1.2, np.pi / 2, 1.6, 2.7, np.pi, 4.0, 7.0])
        F, E, _ = elliptic.elliptic12(phi, m)
        np.testing.assert_allclose(F, ellipkinc(phi, m), atol=1e-12)
        np.testing.assert_allclose(E, ellipeinc(phi, m), atol=1e-12)


class TestElliptic12iBranch:
    """Issue #35: elliptic12i branched at Re(u) = pi/2."""

    def test_continuous_across_half_period(self):
        # The branch point sits at pi/2 + i*acosh(1/sqrt(m)); below it the
        # function is analytic, so crossing Re(u) = pi/2 must be smooth.
        m = 0.4
        psi = 0.5                            # < acosh(1/sqrt(0.4)) = 1.0317
        u = np.linspace(0.2, 3.0, 141) + 1j * psi
        Fi, Ei, _ = elliptic.elliptic12i(u, m)
        assert np.max(np.abs(np.diff(Fi))) < 0.05   # ~0.03 smooth, ~K if branching
        assert np.max(np.abs(np.diff(Ei))) < 0.05

    def test_half_period_line_closed_form(self):
        """At Re(u) = pi/2 exactly the imaginary part used to collapse to 0.

        For |psi| < acosh(1/sqrt(m)):
            F(pi/2 + i*psi | m) = K(m) + i*F(mu | 1-m),
            sin(mu) = tanh(psi)/sqrt(1-m)
        """
        for m in (0.1, 0.4, 0.75):
            for psi in (0.01, 0.05, 0.3):
                mu = np.arcsin(np.tanh(psi) / np.sqrt(1 - m))
                want = ellipk(m) + 1j * ellipkinc(mu, 1 - m)
                got = np.atleast_1d(elliptic.elliptic12i(np.pi / 2 + 1j * psi, m)[0])[0]
                assert abs(got - want) < 1e-12, f"m={m}, psi={psi}: {got} != {want}"

    def test_real_axis_matches_elliptic12(self):
        m = 0.4
        phi = np.array([0.3, 1.2, np.pi / 2, 2.0, 3.4])
        Fi, Ei, _ = elliptic.elliptic12i(phi + 0j, m)
        np.testing.assert_allclose(np.real(Fi), ellipkinc(phi, m), atol=1e-9)
        np.testing.assert_allclose(np.real(Ei), ellipeinc(phi, m), atol=1e-9)


class TestEllipj:
    def test_scalar_vs_scipy(self):
        for u, m in [(0.5, 0.7), (1.0, 0.3), (2.0, 0.5), (0.1, 0.99)]:
            sn, cn, dn, am = elliptic.ellipj(u, m)
            sn_r, cn_r, dn_r, am_r = sp_ellipj(u, m)
            assert abs(sn - sn_r) < 1e-12, f"sn err u={u} m={m}"
            assert abs(cn - cn_r) < 1e-12
            assert abs(dn - dn_r) < 1e-12

    def test_identity_sn2_cn2(self):
        u = np.linspace(0.0, 2.0, 40)
        m = 0.7
        sn, cn, dn, _ = elliptic.ellipj(u, m)
        np.testing.assert_allclose(sn ** 2 + cn ** 2, 1.0, atol=1e-14)

    def test_identity_dn2_m_sn2(self):
        u = np.linspace(0.0, 2.0, 40)
        m = 0.6
        sn, cn, dn, _ = elliptic.ellipj(u, m)
        np.testing.assert_allclose(dn ** 2 + m * sn ** 2, 1.0, atol=1e-14)

    def test_m0_special(self):
        u = np.array([0.3, 1.0, 2.0])
        sn, cn, dn, am = elliptic.ellipj(u, 0.0)
        np.testing.assert_allclose(sn, np.sin(u), atol=1e-15)
        np.testing.assert_allclose(cn, np.cos(u), atol=1e-15)
        np.testing.assert_allclose(dn, 1.0, atol=1e-15)

    def test_m1_special(self):
        u = np.array([0.5, 1.0, 1.5])
        sn, cn, dn, _ = elliptic.ellipj(u, 1.0)
        np.testing.assert_allclose(sn, np.tanh(u), atol=1e-15)
        np.testing.assert_allclose(cn, 1.0 / np.cosh(u), atol=1e-15)
        np.testing.assert_allclose(dn, 1.0 / np.cosh(u), atol=1e-15)
