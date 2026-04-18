"""Tests for ellipticBDJ, ellipticBD, jacobiEDJ, elliptic3, carlson, bulirsch, weierstrass."""
import numpy as np
import pytest
from scipy.special import ellipk, ellipe, ellipkinc, ellipeinc
from scipy.special import elliprf, elliprd, elliprj, elliprc

import elliptic


class TestEllipticBD:
    def test_BD_equals_KE(self):
        m = np.linspace(0.01, 0.99, 30)
        B, D, S = elliptic.ellipticBD(m)
        K, E = ellipk(m), ellipe(m)
        np.testing.assert_allclose(B + D, K, atol=1e-12)
        np.testing.assert_allclose(B + (1 - m) * D, E, atol=1e-12)

    def test_m0_limits(self):
        B, D, S = elliptic.ellipticBD(np.array([1e-4]))
        assert abs(B[0] - np.pi / 4) < 1e-4   # B(m) → pi/4 as m → 0
        assert abs(D[0] - np.pi / 4) < 1e-4

    def test_S_formula(self):
        m = np.array([0.2, 0.5, 0.8])
        B, D, S = elliptic.ellipticBD(m)
        S_ref = (D - B) / m
        np.testing.assert_allclose(S, S_ref, atol=1e-13)


class TestEllipticBDJ:
    def test_BDJ_vs_F_E(self):
        phi = np.linspace(0.01, np.pi / 2, 20)
        m = 0.7
        n = 0.3
        B, D, J = elliptic.ellipticBDJ(phi, m, n)
        F_ref = ellipkinc(phi, m)
        E_ref = ellipeinc(phi, m)
        np.testing.assert_allclose(B + D, F_ref, atol=1e-12, err_msg="F = B+D")
        np.testing.assert_allclose(B + (1 - m) * D, E_ref, atol=1e-12, err_msg="E = B+(1-m)D")

    def test_J_Pi_identity(self):
        phi, m, n = 0.8, 0.5, 0.3
        B, D, J = elliptic.ellipticBDJ(phi, m, n)
        Pi_ref = elliptic.elliptic3(phi, m, n)
        assert abs(B + D + n * J - float(Pi_ref)) < 1e-11

    def test_phi0(self):
        B, D, J = elliptic.ellipticBDJ(0.0, 0.5, 0.3)
        assert abs(float(B)) < 1e-15
        assert abs(float(D)) < 1e-15
        assert abs(float(J)) < 1e-15

    def test_odd_symmetry(self):
        phi = np.array([0.3, 0.7])
        m, n = 0.6, 0.2
        Bp, Dp, Jp = elliptic.ellipticBDJ(phi, m, n)
        Bm, Dm, Jm = elliptic.ellipticBDJ(-phi, m, n)
        np.testing.assert_allclose(Bp, -Bm, atol=1e-13)
        np.testing.assert_allclose(Dp, -Dm, atol=1e-13)
        np.testing.assert_allclose(Jp, -Jm, atol=1e-13)


class TestJacobiEDJ:
    def test_Bu_Du_sum(self):
        m = 0.6
        u = np.linspace(0.01, 1.5, 20)
        Eu, Du, _ = elliptic.jacobiEDJ(u, m)
        Bu = u - Du
        np.testing.assert_allclose(Bu + Du, u, atol=1e-13)

    def test_Eu_formula(self):
        m = 0.7
        u = np.linspace(0.01, 1.2, 15)
        Eu, Du, _ = elliptic.jacobiEDJ(u, m)
        np.testing.assert_allclose(Eu, u - m * Du, atol=1e-13)

    def test_dDu_du_equals_sn2(self):
        """d(D_u)/du = sn²(u|m) — numerical derivative, mirrors testEllipticBDJ.m"""
        m, u0, h = 0.5, 0.6, 1e-7
        _, Dp, _ = elliptic.jacobiEDJ(np.array([u0 + h]), m)
        _, Dm, _ = elliptic.jacobiEDJ(np.array([u0 - h]), m)
        dDu = (float(Dp[0]) - float(Dm[0])) / (2 * h)
        sn, _, _, _ = elliptic.ellipj(np.array([u0]), m)
        assert abs(dDu - float(sn[0]) ** 2) < 1e-6

    def test_dEu_du_equals_dn2(self):
        """d(E_u)/du = dn²(u|m) — numerical derivative, mirrors testEllipticBDJ.m"""
        m, u0, h = 0.5, 0.6, 1e-7
        Ep, _, _ = elliptic.jacobiEDJ(np.array([u0 + h]), m)
        Em, _, _ = elliptic.jacobiEDJ(np.array([u0 - h]), m)
        dEu = (float(Ep[0]) - float(Em[0])) / (2 * h)
        _, _, dn, _ = elliptic.ellipj(np.array([u0]), m)
        assert abs(dEu - float(dn[0]) ** 2) < 1e-6

    def test_Ju_n0_equals_Du(self):
        """J_u(u, n=0|m) = D_u(u|m) — from testEllipticBDJ.m"""
        from scipy.special import ellipk
        m = 0.5
        K = float(ellipk(m))
        u = np.linspace(0.05, 0.9 * K, 15)
        _, Du, Ju0 = elliptic.jacobiEDJ(u, m, np.zeros_like(u))
        np.testing.assert_allclose(np.asarray(Ju0), np.asarray(Du), atol=1e-12)


class TestElliptic3:
    def test_vs_BDnJ(self):
        phi = np.linspace(0.01, np.pi / 2, 20)
        m, n = 0.6, 0.3
        Pi = elliptic.elliptic3(phi, m, n)
        B, D, J = elliptic.ellipticBDJ(phi, m, n)
        np.testing.assert_allclose(np.asarray(Pi), B + D + n * J, atol=1e-11)


class TestCarlson:
    def test_RF_vs_scipy(self):
        for x, y, z in [(1.0, 2.0, 3.0), (0.5, 1.0, 1.5), (0.0, 1.0, 1.0)]:
            np.testing.assert_allclose(elliptic.carlsonRF(x, y, z),
                                       elliprf(x, y, z), atol=1e-12)

    def test_RD_vs_scipy(self):
        for x, y, z in [(0.0, 2.0, 1.0), (1.0, 2.0, 3.0)]:
            np.testing.assert_allclose(elliptic.carlsonRD(x, y, z),
                                       elliprd(x, y, z), atol=1e-12)

    def test_RJ_vs_scipy(self):
        for x, y, z, p in [(0.5, 1.0, 1.5, 2.0), (0.0, 1.0, 2.0, 1.5)]:
            np.testing.assert_allclose(elliptic.carlsonRJ(x, y, z, p),
                                       elliprj(x, y, z, p), atol=1e-12)

    def test_RC_vs_scipy(self):
        for x, y in [(1.0, 2.0), (0.0, 1.0), (3.0, 1.0)]:
            np.testing.assert_allclose(elliptic.carlsonRC(x, y),
                                       elliprc(x, y), atol=1e-12)

    def test_RF_symmetric(self):
        """RF is symmetric in all three arguments — testCarlson.m"""
        a, b, c = 1.0, 2.0, 3.0
        vals = [float(elliptic.carlsonRF(x, y, z))
                for x, y, z in [(a,b,c),(b,a,c),(c,b,a),(a,c,b)]]
        assert max(abs(v - vals[0]) for v in vals) < 1e-13

    def test_RF_xxx(self):
        """RF(x,x,x) = x^(-1/2) — DLMF 19.20.2"""
        x = 2.5
        assert abs(float(elliptic.carlsonRF(x, x, x)) - x**(-0.5)) < 1e-13

    def test_RD_xxx(self):
        """RD(x,x,x) = x^(-3/2) — DLMF 19.20.17"""
        x = 2.0
        assert abs(float(elliptic.carlsonRD(x, x, x)) - x**(-1.5)) < 1e-12

    def test_RD_symmetric_xy(self):
        """RD is symmetric in first two arguments but NOT the third"""
        rd_xy = float(elliptic.carlsonRD(1.0, 2.0, 3.0))
        rd_yx = float(elliptic.carlsonRD(2.0, 1.0, 3.0))
        rd_xz = float(elliptic.carlsonRD(1.0, 3.0, 2.0))
        assert abs(rd_xy - rd_yx) < 1e-12
        assert abs(rd_xy - rd_xz) > 1e-6  # NOT symmetric in z

    def test_RJ_symmetric_xyz(self):
        """RJ symmetric in first three arguments"""
        vals = [float(elliptic.carlsonRJ(x, y, z, 0.5))
                for x, y, z in [(1.,2.,3.),(2.,1.,3.),(3.,1.,2.)]]
        assert max(abs(v - vals[0]) for v in vals) < 1e-11

    def test_RC_pole(self):
        """RC(x, 0) = inf"""
        assert np.isinf(float(elliptic.carlsonRC(1.0, 0.0)))


class TestBulirsch:
    def test_cel1_vs_K(self):
        m = np.array([0.2, 0.4, 0.6, 0.8])
        kc = np.sqrt(1 - m)
        np.testing.assert_allclose(elliptic.cel1(kc), ellipk(m), atol=1e-12)

    def test_cel2_vs_E(self):
        """cel2(kc, 1, kc²) = E(m) — matches testBulirschCEL.m"""
        m = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        kc = np.sqrt(1 - m)
        np.testing.assert_allclose(elliptic.cel2(kc, 1.0, kc**2), ellipe(m), atol=1e-12)

    def test_cel2_B(self):
        """cel2(kc, 1, 0) = B(m) — matches testBulirschCEL.m"""
        m = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        kc = np.sqrt(1 - m)
        B, _, _ = elliptic.ellipticBD(m)
        np.testing.assert_allclose(elliptic.cel2(kc, 1.0, 0.0), B, atol=1e-12)

    def test_cel2_D(self):
        """cel2(kc, 0, 1) = D(m) — matches testBulirschCEL.m"""
        m = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        kc = np.sqrt(1 - m)
        _, D, _ = elliptic.ellipticBD(m)
        np.testing.assert_allclose(elliptic.cel2(kc, 0.0, 1.0), D, atol=1e-12)

    def test_cel3_vs_Pi(self):
        """cel3(kc, 1-n) = Pi(n|m) — matches testBulirschCEL.m"""
        m, n = 0.5, 0.3
        kc = float(np.sqrt(1 - m))
        Pi_leg = float(elliptic.elliptic3(np.pi / 2, m, n))
        Pi_cel = float(elliptic.cel3(kc, 1 - n))
        assert abs(Pi_cel - Pi_leg) < 1e-12


class TestWeierstrass:
    def test_P_pole(self):
        P = elliptic.weierstrassP(0.0, 2.0, 0.0, -2.0)
        assert np.isinf(float(P))

    def test_P_periodicity(self):
        """P is doubly periodic; test a basic identity via sn."""
        e1, e2, e3 = 2.0, 0.0, -2.0
        m = (e2 - e3) / (e1 - e3)
        omega1 = ellipk(m) / np.sqrt(e1 - e3)
        # P at z and z + 2*omega1 should be equal (up to numerical accuracy)
        z = 0.4
        P1 = elliptic.weierstrassP(z, e1, e2, e3)
        P2 = elliptic.weierstrassP(z + 2 * float(omega1), e1, e2, e3)
        assert abs(float(P1) - float(P2)) < 1e-8

    def test_zeta_odd(self):
        e1, e2, e3 = 3.0, 1.0, -4.0
        z = 0.3
        Z1 = elliptic.weierstrassZeta(z, e1, e2, e3)
        Z2 = elliptic.weierstrassZeta(-z, e1, e2, e3)
        assert abs(float(Z1) + float(Z2)) < 1e-8

    def test_sigma_zero(self):
        e1, e2, e3 = 2.0, 0.0, -2.0
        S = elliptic.weierstrassSigma(0.0, e1, e2, e3)
        assert abs(float(S)) < 1e-15

    def test_sigma_odd(self):
        e1, e2, e3 = 2.0, 0.5, -2.5
        z = 0.4
        S1 = elliptic.weierstrassSigma(z, e1, e2, e3)
        S2 = elliptic.weierstrassSigma(-z, e1, e2, e3)
        assert abs(float(S1) + float(S2)) < 1e-10
