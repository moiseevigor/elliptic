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
