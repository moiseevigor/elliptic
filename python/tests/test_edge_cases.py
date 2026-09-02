"""Edge-case tests anchored on values external to this library.

Every expected value here comes from a source independent of the code under
test:

1. Closed forms from the reference literature — Abramowitz & Stegun (A&S)
   and the NIST Digital Library of Mathematical Functions (DLMF,
   https://dlmf.nist.gov) — evaluated from ``math.gamma`` and ``pi``.
2. Exact functional identities from those references (quasi-period,
   Legendre's relation, Jacobi's theta identity, ...).
3. Independent implementations / methods: ``scipy.special`` and
   ``scipy.integrate.quad``, which share no code with this package.

The cases are the ones this library has historically got wrong: phase beyond
pi/2, exact multiples of pi/2 and their neighbouring floating-point values,
the line Re(u) = pi/2 in the complex plane, and m -> 0 / m -> 1.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import (ellipe, ellipeinc, ellipj as sp_ellipj, ellipk,
                           ellipkinc)

import elliptic


def _s(x):
    """First element of a possibly 0-d / 1-element array result."""
    return complex(np.atleast_1d(x)[0]) if np.iscomplexobj(x) else float(np.atleast_1d(x)[0])


# =====================================================================
# A. Closed-form anchors
# =====================================================================
class TestClosedForms:
    def test_K_and_E_at_m_half(self):
        """K(1/2) = Gamma(1/4)^2/(4 sqrt(pi)); E(1/2) from Legendre's relation.

        Legendre (A&S 17.3.13, DLMF 19.7.1) at k = k' = 1/sqrt(2) reduces to
        2*E*K - K^2 = pi/2, hence E = K/2 + pi/(4K).
        The lemniscate constant omega = 2.6220575542921205 = sqrt(2)*K(1/2).
        """
        K_closed = math.gamma(0.25) ** 2 / (4 * math.sqrt(math.pi))
        E_closed = K_closed / 2 + math.pi / (4 * K_closed)
        assert abs(K_closed - 1.8540746773013723) < 1e-15
        assert abs(E_closed - 1.3506438810476755) < 1e-15
        assert abs(math.sqrt(2) * K_closed - 2.6220575542921205) < 1e-15

        F, E, _ = elliptic.elliptic12(math.pi / 2, 0.5)
        assert abs(_s(F) - K_closed) < 1e-14
        assert abs(_s(E) - E_closed) < 1e-14

    def test_legendre_relation(self):
        """E*K' + E'*K - K*K' = pi/2   [A&S 17.3.13, DLMF 19.7.1]"""
        for m in (0.05, 0.2, 0.5, 0.8, 0.95):
            K, E, _ = (float(np.atleast_1d(v)[0])
                       for v in elliptic.elliptic12(math.pi / 2, m))
            Kp, Ep, _ = (float(np.atleast_1d(v)[0])
                         for v in elliptic.elliptic12(math.pi / 2, 1 - m))
            assert abs(E * Kp + Ep * K - K * Kp - math.pi / 2) < 1e-13, f"m={m}"

    def test_degenerate_parameters(self):
        """m = 0: F = E = phi.  m = 1: F = arctanh(sin phi), E = sin phi [A&S 17.4]"""
        phi = np.array([-2.5, -0.4, 0.0, 0.4, 1.2, 1.5])
        F, E, _ = elliptic.elliptic12(phi, 0.0)
        np.testing.assert_allclose(F, phi, atol=1e-15)
        np.testing.assert_allclose(E, phi, atol=1e-15)

        phi = np.array([-1.2, -0.3, 0.0, 0.3, 1.2])
        F, E, _ = elliptic.elliptic12(phi, 1.0)
        np.testing.assert_allclose(F, np.arctanh(np.sin(phi)), atol=1e-13)
        np.testing.assert_allclose(E, np.sin(phi), atol=1e-13)


# =====================================================================
# B. Quasi-period — the identity issue #35 violated
# =====================================================================
class TestQuasiPeriod:
    def test_F_and_E_quasi_period(self):
        """F(phi + k*pi|m) = F(phi|m) + 2k*K(m); E likewise  [A&S 17.4.3]"""
        for m in (0.1, 0.4, 0.75, 0.99):
            K, Ec = float(ellipk(m)), float(ellipe(m))
            for phi in (-1.3, -0.2, 0.4, 1.1, math.pi / 2):
                F0, E0, _ = (float(np.atleast_1d(v)[0])
                             for v in elliptic.elliptic12(phi, m))
                for k in range(-3, 4):
                    Fk, Ek, _ = (float(np.atleast_1d(v)[0])
                                 for v in elliptic.elliptic12(phi + k * math.pi, m))
                    assert abs(Fk - (F0 + 2 * k * K)) < 1e-12, f"F m={m} phi={phi} k={k}"
                    assert abs(Ek - (E0 + 2 * k * Ec)) < 1e-12, f"E m={m} phi={phi} k={k}"

    def test_half_period_neighbourhood(self):
        """Exact multiples of pi/2 and their float neighbours, vs scipy."""
        m = 0.5
        pts = []
        for j in range(-5, 6):
            b = j * math.pi / 2
            e = math.ulp(max(abs(b), 1.0))
            pts += [b, b + e, b - e, b + 3 * e, b - 3 * e]
        pts = np.array(pts)
        F, E, _ = elliptic.elliptic12(pts, m)
        # scipy's ellipkinc/ellipeinc are an independent implementation
        np.testing.assert_allclose(F, ellipkinc(pts, m), atol=1e-11)
        np.testing.assert_allclose(E, ellipeinc(pts, m), atol=1e-11)

    def test_parity_past_half_period(self):
        phi = np.array([0.3, 1.2, math.pi / 2, 2.4, math.pi, 5.0, 7.7])
        for m in (0.2, 0.6, 0.9):
            Fp, Ep, _ = elliptic.elliptic12(phi, m)
            Fm, Em, _ = elliptic.elliptic12(-phi, m)
            np.testing.assert_allclose(Fp, -np.asarray(Fm), atol=1e-12)
            np.testing.assert_allclose(Ep, -np.asarray(Em), atol=1e-12)


# =====================================================================
# C. Complex phase (issue #35)
# =====================================================================
class TestComplexPhase:
    def test_on_the_line_re_u_equals_half_pi(self):
        """A&S 17.4.11-13 collapse at phi = pi/2, below the branch point at
        pi/2 + i*arccosh(1/sqrt(m)):

            F(pi/2 + i*psi | m) = K(m) + i*F(mu | 1-m),
            sin(mu) = tanh(psi)/sqrt(1-m)
        """
        for m in (0.1, 0.4, 0.75, 0.9):
            psi_c = math.acosh(1 / math.sqrt(m))
            for frac in (0.005, 0.05, 0.3, 0.9):
                psi = frac * psi_c
                mu = math.asin(math.tanh(psi) / math.sqrt(1 - m))
                want = ellipk(m) + 1j * ellipkinc(mu, 1 - m)
                got = _s(elliptic.elliptic12i(math.pi / 2 + 1j * psi, m)[0])
                assert abs(got - want) < 1e-12, f"m={m} psi={psi}: {got} != {want}"
                neg = _s(elliptic.elliptic12i(-(math.pi / 2 + 1j * psi), m)[0])
                assert abs(neg + got) < 1e-12, "F must be odd in u"

    def test_continuity_across_half_periods(self):
        """No jump crossing Re(u) = pi/2 + k*pi while below the branch point."""
        m, psi = 0.4, 0.5          # arccosh(1/sqrt(0.4)) = 1.0317
        phi = np.linspace(-4.5, 4.5, 451)
        Fi, Ei, _ = elliptic.elliptic12i(phi + 1j * psi, m)
        h = phi[1] - phi[0]
        assert np.max(np.abs(np.diff(Fi))) < 5 * h
        assert np.max(np.abs(np.diff(Ei))) < 5 * h

    def test_real_axis_and_sn_roundtrip(self):
        m = 0.4
        phi = np.array([-3.4, -1.2, 0.3, math.pi / 2, 2.0, 3.4])
        Fi, Ei, _ = elliptic.elliptic12i(phi + 0j, m)
        np.testing.assert_allclose(np.real(Fi), ellipkinc(phi, m), atol=1e-9)
        np.testing.assert_allclose(np.real(Ei), ellipeinc(phi, m), atol=1e-9)
        # sn(F(u|m)|m) = sin(u) — definition of the amplitude, DLMF 22.16
        u = np.array([0.4 + 0.3j, 1.2 - 0.8j, 2.9 + 0.2j, -2.0 + 0.6j])
        F, _, _ = elliptic.elliptic12i(u, m)
        sn, _, _ = elliptic.ellipji(F, m)
        np.testing.assert_allclose(sn, np.sin(u), atol=1e-9)


# =====================================================================
# D. Jacobi amplitude — DLMF 22.16
# =====================================================================
class TestAmplitude:
    def test_amplitude_landmarks_and_period(self):
        """am(0)=0, am(K)=pi/2, am(2K)=pi, am(u+2K)=am(u)+pi  [DLMF 22.16]"""
        for m in (0.2, 0.6, 0.95):
            K = float(ellipk(m))
            assert abs(_s(elliptic.ellipj(0.0, m)[3])) < 1e-14
            assert abs(_s(elliptic.ellipj(K, m)[3]) - math.pi / 2) < 1e-12
            assert abs(_s(elliptic.ellipj(2 * K, m)[3]) - math.pi) < 1e-12
            for u in (-3.3, -0.7, 0.5, 2.2, 6.0):
                sn, cn, _, am = (_s(v) for v in elliptic.ellipj(u, m))
                am2 = _s(elliptic.ellipj(u + 2 * K, m)[3])
                assert abs(am2 - (am + math.pi)) < 1e-12, f"u={u} m={m}"
                assert abs(math.sin(am) - sn) < 1e-13
                assert abs(math.cos(am) - cn) < 1e-13

    def test_amplitude_matches_scipy_over_many_periods(self):
        """scipy's ellipj returns the continuous amplitude as its 4th output."""
        for m in (0.3, 0.9, 0.999):
            u = np.linspace(-40, 40, 801)
            am = elliptic.ellipj(u, m)[3]
            np.testing.assert_allclose(am, sp_ellipj(u, m)[3], atol=1e-12)


# =====================================================================
# E. Associate integrals B, D, J — DLMF 19.25
# =====================================================================
class TestAssociates:
    def test_connection_and_quadrature_past_half_period(self):
        n = 0.3
        for m in (0.3, 0.8):
            for phi in (0.7, math.pi / 2, 2.0, math.pi, 4.5, 2 * math.pi, 7.0):
                B, D, J = (float(np.atleast_1d(v)[0])
                           for v in elliptic.ellipticBDJ(phi, m, n))
                F, E = float(ellipkinc(phi, m)), float(ellipeinc(phi, m))
                assert abs(B + D - F) < 1e-12, f"B+D != F  phi={phi} m={m}"
                assert abs(B + (1 - m) * D - E) < 1e-12, f"B+(1-m)D != E  phi={phi} m={m}"
                d = lambda t: math.sqrt(1 - m * math.sin(t) ** 2)
                Bq = quad(lambda t: math.cos(t) ** 2 / d(t), 0, phi, limit=200)[0]
                Dq = quad(lambda t: math.sin(t) ** 2 / d(t), 0, phi, limit=200)[0]
                Jq = quad(lambda t: math.sin(t) ** 2 /
                          ((1 - n * math.sin(t) ** 2) * d(t)), 0, phi, limit=200)[0]
                assert abs(B - Bq) < 1e-10, f"B({phi}|{m}) = {B} vs {Bq}"
                assert abs(D - Dq) < 1e-10, f"D({phi}|{m}) = {D} vs {Dq}"
                assert abs(J - Jq) < 1e-10, f"J({phi},{n}|{m}) = {J} vs {Jq}"

    def test_associate_quasi_period(self):
        """All three integrands are pi-periodic: X(phi+k*pi) = X(phi) + 2k*X(m)."""
        m, n = 0.6, 0.25
        Bc, Dc, Jc = (float(np.atleast_1d(v)[0])
                      for v in elliptic.ellipticBDJ(math.pi / 2, m, n))
        for phi in (-1.1, 0.0, 0.9, math.pi / 2):
            B0, D0, J0 = (float(np.atleast_1d(v)[0])
                          for v in elliptic.ellipticBDJ(phi, m, n))
            for k in range(-2, 3):
                Bk, Dk, Jk = (float(np.atleast_1d(v)[0])
                              for v in elliptic.ellipticBDJ(phi + k * math.pi, m, n))
                assert abs(Bk - (B0 + 2 * k * Bc)) < 1e-12, f"B phi={phi} k={k}"
                assert abs(Dk - (D0 + 2 * k * Dc)) < 1e-12, f"D phi={phi} k={k}"
                assert abs(Jk - (J0 + 2 * k * Jc)) < 1e-12, f"J phi={phi} k={k}"

    def test_jacobi_argument_associates(self):
        """E_u(u|m) = E(am(u|m)|m)  [DLMF 22.16]; D_u against quadrature."""
        for m in (0.3, 0.8):
            K = float(ellipk(m))
            for u in (0.5, K, 1.5 * K, 2 * K, 3 * K, 9.0):
                Eu, Du, _ = elliptic.jacobiEDJ(u, m)
                am = _s(elliptic.ellipj(u, m)[3])
                assert abs(_s(Eu) - ellipeinc(am, m)) < 1e-10, f"u={u} m={m}"
                Dq = quad(lambda v: sp_ellipj(v, m)[0] ** 2, 0, u, limit=400)[0]
                assert abs(_s(Du) - Dq) < 1e-9, f"D_u({u}|{m})"


# =====================================================================
# F. Inverse of the incomplete second kind
# =====================================================================
class TestInverseE:
    def test_roundtrip_over_several_periods(self):
        for m in (0.1, 0.5, 0.9, 0.99):
            E1 = float(ellipe(m))
            z = E1 * np.array([-3, -2, -1, -0.5, 0, 0.01, 0.5, 0.999, 1, 1.5, 2, 3.0])
            phi = elliptic.inverselliptic2(z, m)
            np.testing.assert_allclose(ellipeinc(phi, m), z, atol=1e-11)
            np.testing.assert_allclose(elliptic.inverselliptic2(-z, m), -np.asarray(phi),
                                       atol=1e-11)


# =====================================================================
# G. Complete third kind — DLMF 19.6
# =====================================================================
class TestThirdKind:
    def test_special_values(self):
        """Pi(0,m) = K(m);  Pi(n=m, m) = E(m)/(1-m)  [DLMF 19.6]"""
        for m in (0.1, 0.4, 0.8):
            K, E = float(ellipk(m)), float(ellipe(m))
            assert abs(_s(elliptic.elliptic3(math.pi / 2, m, 0.0)) - K) < 1e-12
            assert abs(_s(elliptic.elliptic3(math.pi / 2, m, m)) - E / (1 - m)) < 1e-11


# =====================================================================
# H. Theta functions — A&S 16.38, DLMF 20.7
# =====================================================================
class TestTheta:
    def test_null_values_and_jacobi_identity(self):
        """theta2(0)^2 = 2 sqrt(m) K/pi, theta3(0)^2 = 2K/pi,
        theta4(0)^2 = 2 sqrt(1-m) K/pi  [A&S 16.38];
        theta3^4 = theta2^4 + theta4^4  [DLMF 20.7.1]."""
        for m in (0.05, 0.3, 0.5, 0.9):
            K = float(ellipk(m))
            t1, t2, t3, t4 = (_s(elliptic.theta(j, 0.0, m)) for j in (1, 2, 3, 4))
            assert abs(t2 ** 2 - 2 * math.sqrt(m) * K / math.pi) < 1e-13, f"m={m}"
            assert abs(t3 ** 2 - 2 * K / math.pi) < 1e-13, f"m={m}"
            assert abs(t4 ** 2 - 2 * math.sqrt(1 - m) * K / math.pi) < 1e-13, f"m={m}"
            assert abs(t3 ** 4 - (t2 ** 4 + t4 ** 4)) < 1e-13, f"m={m}"
            assert abs(t1) < 1e-14

    def test_theta_eta_versus_sn(self):
        """sn(u|m) = H(u)/(m^(1/4) Theta(u))  [A&S 16.31], past the half period."""
        m = 0.5
        K = float(ellipk(m))
        for u in (0.3, K, 1.5 * K, 2 * K, 3 * K, 5.0, -2.2):
            Th, H = (_s(v) for v in elliptic.jacobiThetaEta(u, m))
            sn = float(sp_ellipj(u, m)[0])
            assert abs(H / (m ** 0.25 * Th) - sn) < 1e-12, f"u/K={u / K}"


# =====================================================================
# I. Carlson symmetric forms — DLMF 19.20
# =====================================================================
class TestCarlson:
    def test_special_values(self):
        for x in (0.25, 1.0, 3.7):
            assert abs(_s(elliptic.carlsonRF(x, x, x)) - x ** -0.5) < 1e-14
            assert abs(_s(elliptic.carlsonRD(x, x, x)) - x ** -1.5) < 1e-14
            assert abs(_s(elliptic.carlsonRJ(x, x, x, x)) - x ** -1.5) < 1e-14
            assert abs(_s(elliptic.carlsonRC(x, x)) - x ** -0.5) < 1e-14
        assert abs(_s(elliptic.carlsonRF(0.0, 1.0, 1.0)) - math.pi / 2) < 1e-14
        assert abs(_s(elliptic.carlsonRC(0.0, 1.0)) - math.pi / 2) < 1e-14

    def test_legendre_connection(self):
        """R_F(0,1-m,1) = K(m);  R_D(0,1-m,1)/3 = D(m) = (K-E)/m  [DLMF 19.25]"""
        for m in (0.1, 0.5, 0.9):
            K, E = float(ellipk(m)), float(ellipe(m))
            assert abs(_s(elliptic.carlsonRF(0.0, 1 - m, 1.0)) - K) < 1e-13
            assert abs(_s(elliptic.carlsonRD(0.0, 1 - m, 1.0)) / 3 - (K - E) / m) < 1e-12


# =====================================================================
# K. Weierstrass functions — 30-digit external reference values
# =====================================================================
class TestWeierstrassExternal:
    """Reference: theta-function closed forms (DLMF 23.6.4/8/9/13) evaluated
    with mpmath 1.4.1 at mp.dps = 50 — a formula family independent of this
    library's sn-based / series code.  The reference construction was
    validated to < 1e-44 against P'^2 = 4P^3 - g2 P - g3, zeta' = -P,
    sigma'/sigma = zeta, P(w1) = e1, and the z -> 0 Laurent behaviour.

    Roots are exactly binary-representable with e1+e2+e3 == 0 exactly
    (inexact roots break the depressed cubic by O(eps*P^2)).
    z/w1 = 2.6 is the case the old quadrature Sigma got wrong by -277x.
    """

    # e1, e2, e3, z, P, Pp, Zeta, Sigma  (mpmath, 30 digits)
    ROWS = [
        (2.0, -0.5, -1.5, 0.364678508599203802930916253782,
         7.60991654315483213781676575283, -40.716795203467035857794498067,
         2.73133905217586401231118481172, 0.364322893731329436910487536054),
        (2.0, -0.5, -1.5, 1.09403552579761140879274876135,
         2.31139666938461354617370617507, 3.65334019657835372774485143298,
         0.509456906969716181254227726116, 0.993298793310168551601976056672),
        (2.0, -0.5, -1.5, 2.37041030589482471905095564958,
         3.55973799836232295109911735041, -11.3205842122331758742412322161,
         3.5753186226267378852136602119, -7.35675884743112892037462320321),
        (1.0, 0.25, -1.25, 0.54105576073301731403993176119,
         3.48955873892538680757168974091, -12.3652547699847269163302027234,
         1.83475002192072477467486916472, 0.540061003157461700899655935519),
        (1.0, 0.25, -1.25, 1.62316728219905194211979528357,
         1.13305294818566893955539247187, 1.05828457240220452238075426776,
         0.25696621923989724678288838394, 1.40818598628919804944444268506),
        (1.0, 0.25, -1.25, 3.51686244476461254125955644774,
         1.67783096402393224759843061699, -3.36668198668738980119457645764,
         2.26609152750196941530142216878, -8.29563954795637078192140811817),
    ]

    def test_against_mpmath_reference(self):
        for e1, e2, e3, z, P, Pp, Ze, Si in self.ROWS:
            assert abs(_s(elliptic.weierstrassP(z, e1, e2, e3)) - P) < 1e-12 * max(1, abs(P))
            assert abs(_s(elliptic.weierstrassPPrime(z, e1, e2, e3)) - Pp) < 1e-12 * max(1, abs(Pp))
            assert abs(_s(elliptic.weierstrassZeta(z, e1, e2, e3)) - Ze) < 1e-12 * max(1, abs(Ze))
            assert abs(_s(elliptic.weierstrassSigma(z, e1, e2, e3)) - Si) < 1e-12 * max(1, abs(Si))

    def test_structural_identities(self):
        """DE, zeta quasi-period (eta1 = zeta(w1)), sigma quasi-period
        sigma(z + 2 w1) = -sigma(z) exp(2 eta1 (z + w1))   [DLMF 23.2]"""
        e1, e2, e3 = 1.5, 0.25, -1.75
        g2, g3, _ = (float(np.atleast_1d(v)[0])
                     for v in elliptic.weierstrassInvariants(e1, e2, e3))
        m = (e2 - e3) / (e1 - e3)
        w1 = float(ellipk(m)) / math.sqrt(e1 - e3)
        eta1 = _s(elliptic.weierstrassZeta(w1, e1, e2, e3))
        for zf in (0.23, 0.61, 0.89):
            z = zf * w1
            P = _s(elliptic.weierstrassP(z, e1, e2, e3))
            Pp = _s(elliptic.weierstrassPPrime(z, e1, e2, e3))
            assert abs(Pp ** 2 - (4 * P ** 3 - g2 * P - g3)) < 1e-9 * max(1, Pp ** 2)
            Z0 = _s(elliptic.weierstrassZeta(z, e1, e2, e3))
            for k in (-2, -1, 1, 2):
                Zk = _s(elliptic.weierstrassZeta(z + 2 * k * w1, e1, e2, e3))
                assert abs(Zk - (Z0 + 2 * k * eta1)) < 1e-11, f"zeta q-p z/w1={zf} k={k}"
            Sq = _s(elliptic.weierstrassSigma(z + 2 * w1, e1, e2, e3))
            Sw = -_s(elliptic.weierstrassSigma(z, e1, e2, e3)) * math.exp(2 * eta1 * (z + w1))
            assert abs(Sq - Sw) < 1e-10 * max(1, abs(Sw)), f"sigma q-p z/w1={zf}"


# =====================================================================
# J. Ellipse arclength — spans 2*pi of phase, so exercises the quasi-period
# =====================================================================
class TestArclength:
    def test_circle_and_full_ellipse(self):
        assert abs(_s(elliptic.arclength_ellipse(2.0, 2.0)) - 2 * math.pi * 2) < 1e-12
        for a, b in ((5.0, 10.0), (10.0, 5.0), (1.0, 0.5), (1.0, 2.0)):
            ref = quad(lambda t: math.sqrt(a ** 2 * math.sin(t) ** 2
                                           + b ** 2 * math.cos(t) ** 2),
                       0, 2 * math.pi, limit=200)[0]
            assert abs(_s(elliptic.arclength_ellipse(a, b)) - ref) < 1e-9, f"a={a} b={b}"


# =====================================================================
# L. Third kind past pi/2, complex Jacobi, nome, agm, theta_prime —
#    mpmath 1.4.1 references at mp.dps = 50 (25 significant digits)
# =====================================================================
class TestMpmathAnchors:
    def test_elliptic3_extended_phase(self):
        """mpmath ellippi(c, u, m); exercises the quasi-period reduction."""
        rows = [  # u, c, m, Pi
            (2.0, 0.3, 0.4, 2.906928215899509476566159),
            (3.14159265358979323846, 0.3, 0.4, 4.297590831732913271869402),
            (6.0, 0.3, 0.4, 8.30819635557621288613689),
            (-2.2, 0.25, 0.6, -3.471643753410082592812564),
            (4.5, 0.7, 0.2, 8.475833763417048864937332),
        ]
        for u, c, m, ref in rows:
            got = _s(elliptic.elliptic3(u, m, c))
            assert abs(got - ref) < 1e-12 * max(1, abs(ref)), f"Pi({u}|{m},{c})"
        assert abs(_s(elliptic.elliptic3(-2.0, 0.4, 0.3))
                   + _s(elliptic.elliptic3(2.0, 0.4, 0.3))) < 1e-13

    def test_ellipji_vs_mpmath(self):
        rows = [  # m, u, sn, cn, dn   (mpmath ellipfun)
            (0.3, complex(1.2, -0.8), complex(1.144487529338893062450354, -0.2802393491508895255196492), complex(0.4746455171789828435720278, 0.6757262603879129546403857), complex(0.8030935803341932691444772, 0.1198106104396107763170419)),
            (0.7, complex(2.9, 0.2), complex(0.8945661685904884966950071, -0.0623052446149157017922886), complex(-0.4667874916734775923001321, -0.1194037221486780885477261), complex(0.667799853302213119590502, 0.05842366478197534982879624)),
            (0.7, complex(0.1, 1.9), complex(2.711675412858391500357814, -4.920359160836619695635427), complex(-4.998252993291208495782488, -2.669416089337967774258619), complex(-4.209736402055241155323296, -2.218593037476527723830424)),
            (0.96, complex(4.5, -2.5), complex(0.9833643701168131893646114, -0.09200718087257114311996756), complex(0.3369732007632264400404546, 0.2684978605421883948528415), complex(-0.3680824599803019014518416, -0.2359729940161469769037399)),
        ]
        for m, u, sn_r, cn_r, dn_r in rows:
            sn, cn, dn = (complex(np.atleast_1d(v)[0]) for v in elliptic.ellipji(u, m))
            assert abs(sn - sn_r) < 1e-12 * max(1, abs(sn_r))
            assert abs(cn - cn_r) < 1e-12 * max(1, abs(cn_r))
            assert abs(dn - dn_r) < 1e-12 * max(1, abs(dn_r))

    def test_nome_agm_theta_prime(self):
        """q(1/2) = exp(-pi) exactly (K = K' at m = 1/2); AGM and jtheta rows."""
        assert abs(_s(elliptic.nomeq(0.5)) - math.exp(-math.pi)) < 1e-16
        assert abs(_s(elliptic.inversenomeq(math.exp(-math.pi))) - 0.5) < 1e-12
        assert abs(_s(elliptic.agm(1.0, 0.01)) - 0.2621668872022492366947771) < 1e-15
        th, thp = (_s(v) for v in elliptic.theta_prime(1, 0.4, 0.6))
        assert abs(th - 0.3776251831225481533690215) < 1e-13
        assert abs(thp - 0.8967180488866961938686533) < 1e-13

    def test_theta_prime_at_theta_zeros(self):
        """DLMF 20.4.6: theta1'(0) = theta2(0)*theta3(0)*theta4(0), and the
        derivative must be finite at the zeros of theta itself (the MATLAB
        port's log-derivative form returned NaN there)."""
        for m in (0.1, 0.5, 0.9):
            d1 = _s(elliptic.theta_prime(1, 0.0, m)[1])
            t2, t3, t4 = (_s(elliptic.theta_prime(j, 0.0, m)[0]) for j in (2, 3, 4))
            assert abs(d1 - t2 * t3 * t4) < 1e-13, f"DLMF 20.4.6 at m={m}"
            dpi = _s(elliptic.theta_prime(1, math.pi, m)[1])
            assert not math.isnan(dpi) and abs(dpi + d1) < 1e-13
            assert not math.isnan(_s(elliptic.theta_prime(2, math.pi / 2, m)[1]))


# =====================================================================
# Q. Adversarial-review round (external Codex + mpmath 1.4.1, dps=40).
#    Each test is a counterexample a prior version failed.
# =====================================================================
class TestAdversarialRound:
    def test_elliptic3_negative_amplitude_with_pole(self):
        """0*Inf guard: negative phase never crossing the complete-integral pole."""
        assert abs(_s(elliptic.elliptic3(-1.0, 0.5, 1.0)) - (-1.7319915420235269928)) < 1e-13
        assert abs(_s(elliptic.elliptic3(-1.0, 1.0, 0.2)) - (-1.3115010674599590753)) < 1e-13
        assert not math.isnan(_s(elliptic.elliptic3(-0.4, 1.0, 1.0)))

    def test_complex_FE_small_m_series(self):
        """A&S 17.4.11 path lost sqrt(eps/m) digits; the m^2 series is exact."""
        assert abs(_s(elliptic.elliptic12i(0.2j, 1e-20)[0]) - 0.2j) < 1e-15
        F, E, _ = elliptic.elliptic12i(math.pi / 2 + 0.2j, 1e-14)
        assert abs(_s(F) - (1.5707963267949005462 + 0.20000000000000101344j)) < 1e-13
        assert abs(_s(E) - (1.5707963267948926922 + 0.19999999999999898656j)) < 1e-13
        F, E, _ = elliptic.elliptic12i(math.pi / 2 + 0.2j, 1e-6)
        assert abs(_s(F) - (1.5707967194941992113 + 0.20000010134411776594j)) < 5e-12
        assert abs(_s(E) - (1.5707959340957412894 + 0.19999989865593359446j)) < 5e-12
        # both sides of the series threshold vs mpmath (dps=30)
        Fa = _s(elliptic.elliptic12i(1.1 + 0.3j, 0.99e-4)[0])
        assert abs(Fa - (1.1000153646885162 + 0.3000120622623928j)) < 1e-12
        Fb = _s(elliptic.elliptic12i(1.1 + 0.3j, 1.01e-4)[0])
        assert abs(Fb - (1.1000156750952820 + 0.3000123059589823j)) < 5e-12

    def test_weierstrass_near_origin_finite(self):
        """DLMF 23.9.2: only the exact lattice point is a pole; z = 1e-16 is
        a huge FINITE value (a tolerance here used to return Inf)."""
        assert abs(_s(elliptic.weierstrassP(1e-16, 1.0, 0.0, -1.0)) - 1e32) < 1e19
        assert abs(_s(elliptic.weierstrassPPrime(1e-16, 1.0, 0.0, -1.0)) + 2e48) < 1e36
        assert abs(_s(elliptic.weierstrassZeta(1e-16, 1.0, 0.0, -1.0)) - 1e16) < 1e3
        assert math.isinf(_s(elliptic.weierstrassP(0.0, 1.0, 0.0, -1.0)))

    def test_inverse_nome_all_scales(self):
        """DLMF 20.9.1 closed form: m = (theta2/theta3)^4, exact at every scale
        (the 64-step bisection had a 2^-64 absolute floor: m(1e-30) came back
        2.7e-20)."""
        assert abs(_s(elliptic.inversenomeq(1e-30)) - 1.6e-29) < 1e-41
        assert abs(_s(elliptic.inversenomeq(1e-12)) - 1.5999999999872e-11) < 1e-24
        for mv in (1e-8, 0.3, 0.85, 0.999):
            assert abs(_s(elliptic.inversenomeq(np.asarray(_s(elliptic.nomeq(mv))))) - mv) \
                < 1e-12 * max(mv, 1e-3), f"roundtrip at m={mv}"
        # the computed upper endpoint must be accepted, not rejected
        q_max = _s(elliptic.nomeq(np.nextafter(1.0, 0.0)))
        assert _s(elliptic.inversenomeq(q_max)) > 0.999

    def test_carlson_scale_invariance(self):
        """DLMF 19.20: RF, RC ~ lambda^-1/2; RD, RJ ~ lambda^-3/2.  An absolute
        branch tolerance in RC broke this at small scales (27% at 1e-20)."""
        x, y, z, p = 1.0, 2.0, 3.0, 4.0
        for lam in (1e-20, 1e20):
            for fn, args, power in (
                (elliptic.carlsonRF, (x, y, z), 0.5),
                (elliptic.carlsonRC, (x, y), 0.5),
                (elliptic.carlsonRD, (x, y, z), 1.5),
                (elliptic.carlsonRJ, (x, y, z, p), 1.5),
            ):
                base = _s(fn(*args))
                scaled = _s(fn(*(lam * a for a in args)))
                want = base / lam ** power
                assert abs(scaled - want) < 1e-10 * abs(want), f"{fn.__name__} at {lam}"
        assert abs(_s(elliptic.carlsonRC(1e-20, 2e-20)) - 7853981633.9744830962) < 1e-4

    def test_ellipticBD_nondegenerate_anchors(self):
        """mpmath: B = (E-(1-m)K)/m, D = (K-E)/m at dps=40."""
        rows = [(0.2, 0.8066808960371526438, 0.85294270257337535705),
                (0.7, 0.88437375336868858245, 1.1909893819237805614),
                (0.999, 0.99832798626015502386, 3.8428045742901420065)]
        for m, B_ref, D_ref in rows:
            B, D, _ = elliptic.ellipticBD(m)
            assert abs(_s(B) - B_ref) < 1e-14, f"B({m})"
            assert abs(_s(D) - D_ref) < 1e-13, f"D({m})"

    def test_reversed_arc_intervals_signed(self):
        """Reversal negates the arc for ellipses AND circles alike."""
        assert abs(_s(elliptic.arclength_ellipse(2.0, 3.0, 1.0, 0.1))
                   + _s(elliptic.arclength_ellipse(2.0, 3.0, 0.1, 1.0))) < 1e-13
        assert abs(_s(elliptic.arclength_ellipse(2.0, 2.0, 1.0, 0.1)) - (-1.8)) < 1e-13


# =====================================================================
# R. Second adversarial round — fuzz vs mpmath (dps=40) over parameter
#    endpoints, extreme scales, poles and period multiples.  Every
#    reference was evaluated at the EXACT DOUBLE the library receives:
#    near singularities the decimal input and its double rounding differ
#    enough to move the answer at 1e-9.  scipy reaches machine precision
#    on every case below; each was an implementation loss, now fixed.
# =====================================================================
class TestAdversarialRound2:
    M1 = float(np.nextafter(1.0, 0.0))

    def test_F_near_pole_m_to_1_and_m_equal_1(self):
        assert abs(_s(elliptic.elliptic12(math.pi/2 - 1e-9, self.M1)[0]) - 19.65993026560449767) < 1e-12*20
        assert _s(elliptic.elliptic12(0.0, 1.0)[0]) == 0.0
        assert abs(_s(elliptic.elliptic12(1e-16, 1.0)[0]) - 1e-16) < 1e-31

    def test_third_kind_endpoint_poles(self):
        assert abs(_s(elliptic.elliptic3(math.pi/2 - 1e-6, 0.3, 1.0)) - 1195228.2584444625825) < 1e-12*1.2e6
        assert abs(_s(elliptic.elliptic3(math.pi/2 - 1e-6, 1 - 1e-8, 0.9)) - 88.615055050793590585) < 1e-12*90

    def test_carlson_disparate_scales_tiny_y_near_equal_double_zero(self):
        assert abs(_s(elliptic.carlsonRJ(1e-20, 2e-20, 3e-20, 0.5)) - 43616756114.805842986) < 1e-12*4.4e10
        assert abs(_s(elliptic.carlsonRJ(2.0, 3.0, 4.0, 1e-10)) - 7.179193296087372323) < 1e-13*7.2
        assert abs(_s(elliptic.carlsonRC(3.0, 1e-10)) - 7.3643213780616827229) < 1e-13*7.4
        assert abs(_s(elliptic.carlsonRC(1.0000000000001, 1.0)) - 0.99999999999998334665) < 1e-14
        assert math.isinf(_s(elliptic.carlsonRF(0.0, 0.0, 1.0)))
        assert math.isinf(_s(elliptic.carlsonRD(0.0, 0.0, 1.0)))
        assert math.isinf(_s(elliptic.carlsonRJ(0.0, 0.0, 1.0, 2.0)))

    def test_ellipj_m_to_1(self):
        assert abs(_s(elliptic.ellipj(9.375277798108883, self.M1)[1]) - 0.00016958935096417269446) < 1e-13*1.7e-4
        assert abs(_s(elliptic.ellipj(7.0, 1 - 1e-12)[1]) - 0.0018237622775256289351) < 1e-13*1.8e-3

    def test_inverse_E_tiny_negative_z(self):
        for m in (0.0, 0.5, 1 - 1e-8):
            z = -1e-9 * float(ellipe(m))
            phi = _s(elliptic.inverselliptic2(z, m))
            assert abs(float(ellipeinc(phi, m)) - z) < 1e-13 * abs(z), f"m={m}"

    def test_complex_F_tiny_psi_and_small_m(self):
        assert abs(_s(elliptic.elliptic12i(math.pi/2 + 1e-9j, 0.9)[0]).imag - 3.1622776601683798848e-9) < 1e-12*3.2e-9
        F = _s(elliptic.elliptic12i(math.pi/2 + 1e-9j, self.M1)[0])
        assert abs(F - complex(19.754694640120759063, 0.095049319491958534055)) < 1e-9*20
        F = _s(elliptic.elliptic12i(0.4 + 0.3j, 1e-4)[0])
        assert abs(F - complex(0.39999936996865927219, 0.30000195549059365219)) < 1e-13

    def test_weierstrass_near_m1_lattice(self):
        e1, e2, e3, z = 0.5000001, 0.5, -1.0000001, 13.391953465243201
        assert abs(_s(elliptic.weierstrassP(z, e1, e2, e3)) - 0.51848214450600943279) < 1e-12
        assert abs(_s(elliptic.weierstrassZeta(z, e1, e2, e3)) - -5.4787546526901492279) < 1e-12*5.5
        assert abs(_s(elliptic.weierstrassSigma(z, e1, e2, e3)) - 1.822626274365935705e-13) < 1e-12*1.8e-13


# =====================================================================
# S. Third adversarial round: dense random fuzz + API abuse
# =====================================================================
class TestAdversarialRound3:
    def test_cody_waite_reduction_point(self):
        """u = 5.5*pi + 8e-10, m = 1-1.5e-13 (mpmath at the exact doubles): before
        the Cody-Waite split, k*pi rounding cost eps*|u| in the reduced phase,
        amplified ~1e5x by dZ/dphi near pi/2 at m -> 1."""
        u, m = 17.27875959554386, 0.99999999999985
        F, E, Z = elliptic.elliptic12(u, m)
        assert abs(_s(F) - 177.65640489133312311) < 2e-9
        assert abs(_s(E) - 11.000000000012911122) < 1e-12
        assert abs(_s(Z) - (-0.00012790056416609388015)) < 1e-10
        assert abs(_s(elliptic.elliptic3(u, m, 0.3)) - 248.50046674013002377) < 3e-9

    def test_domain_and_nan_are_honest(self):
        """Out-of-range m raises on numpy; NaN propagates and never becomes a
        placeholder value (ellipj(0.3, 1.5) used to return sn(0.3 | 0.5))."""
        for bad in (1.5, -1e-17, np.nextafter(1.0, 2.0)):
            with pytest.raises(ValueError): elliptic.ellipj(0.3, bad)
            with pytest.raises(ValueError): elliptic.elliptic12(0.3, bad)
            with pytest.raises(ValueError): elliptic.nomeq(bad)
        assert math.isnan(_s(elliptic.ellipj(0.3, float('nan'))[0]))
        assert math.isnan(_s(elliptic.elliptic12(0.3, float('nan'))[0]))
        v = elliptic.elliptic12(np.array([0.3, 0.5, 0.7]), np.array([0.2, np.nan, 0.4]))[0]
        assert np.isnan(v[1]) and not np.isnan(v[[0, 2]]).any()

    def test_exact_zero_complex_and_RJ_extreme_ratio(self):
        # note: the Python literal -0.0+0j already evaluates to 0j; use complex(-0.0, 0.0)
        F0 = _s(elliptic.elliptic12i(complex(-0.0, 0.0), 0.5)[0])
        assert F0 == 0 and math.copysign(1.0, F0.real) < 0      # -0.0 preserved
        assert _s(elliptic.elliptic12i(0j, 0.5)[0]) == 0
        assert abs(_s(elliptic.carlsonRJ(2.798e-18, 5.954e-24, 9.634e-23, 1.134e21)) - 9.9678905686736778972e-12) < 1e-12 * 1e-11


# =====================================================================
# T. Theta at a huge argument; Weierstrass root ordering
# =====================================================================
class TestAdversarialRound4:
    def test_theta_huge_argument(self):
        """mpmath jtheta at the exact double v = 123456789.123: the products
        (2n+1)*v rounded by eps*|k v| before the angle-addition recurrence."""
        v = 123456789.123
        t, tp = elliptic.theta_prime(1, v, 0.4)
        assert abs(_s(t) - (0.84585020823346348431)) < 2e-15 and abs(_s(tp) - (0.015114923736622936955)) < 2e-14
        assert abs(_s(elliptic.theta(1, v, 0.4)) - (0.84585020823346348431)) < 2e-15
        t, tp = elliptic.theta_prime(2, v, 0.4)
        assert abs(_s(t) - (0.014932290326334898348)) < 2e-15 and abs(_s(tp) - (-0.84241862816186020729)) < 2e-14
        assert abs(_s(elliptic.theta(2, v, 0.4)) - (0.014932290326334898348)) < 2e-15
        t, tp = elliptic.theta_prime(3, v, 0.4)
        assert abs(_s(t) - (0.93627542467710214194)) < 2e-15 and abs(_s(tp) - (-0.004519191759268069986)) < 2e-14
        assert abs(_s(elliptic.theta(3, v, 0.4)) - (0.93627542467710214194)) < 2e-15
        t, tp = elliptic.theta_prime(4, v, 0.4)
        assert abs(_s(t) - (1.0637286984176921296)) < 2e-15 and abs(_s(tp) - (0.0045203629452149075654)) < 2e-14
        assert abs(_s(elliptic.theta(4, v, 0.4)) - (1.0637286984176921296)) < 2e-15

    def test_weierstrass_root_order_is_enforced(self):
        with pytest.raises(ValueError): elliptic.weierstrassP(0.5, 0.5, 1.0, -1.5)      # unsorted
        with pytest.raises(ValueError): elliptic.weierstrassZeta(0.5, 1.0, 1.0, 1.0)     # e1 == e3
        # equal neighbours are the legitimate degenerate lattices (m = 1 / m = 0)
        assert math.isfinite(_s(elliptic.weierstrassP(0.5, 1.0, 1.0, -2.0)))
        assert math.isfinite(_s(elliptic.weierstrassP(0.5, 1.0, 0.0, 0.0)))


# =====================================================================
# U. Round 6: cross-port parity sweep at extreme m (anchors: mpmath at the
#    exact double inputs).  The k*pi reduction now uses a 25-bit split of pi
#    (k*float(pi) rounded by eps*|u|), the RJ series E3 coefficient is
#    4 P^3 (DLMF 19.36.2), ellipticBDJ rejects n > 1 beyond the pole and
#    handles n = 1, and the Python-only tiny-m / theta-nome paths are pinned.
# =====================================================================
class TestAdversarialRound6:
    def test_tiny_m_band(self):
        F, E, _ = elliptic.elliptic12(np.array([1.0, 1.0]), np.array([3e-16, 5e-16]))
        assert np.all(np.isfinite(F)) and np.all(np.isfinite(E))
        assert abs(F[0] - 1.0) < 5e-16 and abs(E[0] - 0.99999999999999996) < 5e-16
        assert abs(F[1] - 1.0000000000000001) < 5e-16 and abs(E[1] - 0.99999999999999993) < 5e-16

    def test_theta_nome_from_exact_m(self):
        assert abs(_s(elliptic.theta(1, 34401.9, 1.6e-16)) - 0.00011178415088289534) < 1e-13 * 1.1e-4
        assert abs(_s(elliptic.theta_prime(1, 6577.39, 1.5e-16)[0]) - (-9.8878892558450512e-5)) < 1e-13 * 1e-4

    def test_E_near_m1_and_large_u_pi_split(self):
        F, E, _ = elliptic.elliptic12(-1.65181, 0.99999999999999578)
        assert abs(_s(E) - (-1.0032798131910099)) < 5e-14 and abs(_s(F) - (-32.666065762173088)) < 1e-14 * 33
        u, m = 80101.48788857895, 0.9999533239086507          # 25497*pi + 0.3
        F, E, Z = elliptic.elliptic12(u, m)
        assert abs(_s(Z) - 0.2477141143165845) < 2e-14           # eps*|u| = 1.8e-11 before the split
        assert abs(_s(E) - 51001.284415600044) < 1e-14 * 51001
        assert abs(_s(F) - 324959.38078465716) < 1e-14 * 324959
        _, _, Z = elliptic.elliptic12(1000000.123, 1 - 2**-53)
        assert abs(_s(Z) - (-0.220434859492317)) < 2e-14
        assert abs(_s(elliptic.elliptic3(-2.70143, 1 - 2**-53, 0.9723)) - (-1249.3300419938347)) < 1e-14 * 1249

    def test_RJ_series_and_characteristic_domain(self):
        assert abs(_s(elliptic.carlsonRJ(0.1, 0.2, 1, 3.0)) - 1.1311524759367163) < 5e-15 * 1.13
        assert abs(_s(elliptic.carlsonRJ(0.292, 0.646, 1, 1.354)) - 1.2806109121365949) < 5e-15 * 1.28
        assert abs(_s(elliptic.ellipticBDJ(1.0, 0.5, 1.0)[2]) - 0.64877476917835824) < 5e-15      # was NaN (0 * inf)
        assert abs(_s(elliptic.ellipticBDJ(0.5, 0.5, 1.5)[2]) - 0.052791966372572887) < 5e-16
        with pytest.raises(ValueError, match="principal"):
            elliptic.ellipticBDJ(1.0, 0.5, 1.5)                  # beyond the pole: returned 1.147 silently
        for c, ref in [(-0.5, 0.9560406633267465), (-3.0, 0.66684868942035313), (-100.0, 0.1523863772236308)]:
            assert abs(_s(elliptic.elliptic3(1.0, 0.5, c)) - ref) < 5e-16
        assert abs(_s(elliptic.elliptic3(4.0, 0.9, -100.0)) - 0.4921742710224714) < 5e-15

    def test_bulirsch_cel_is_kc_native(self):
        """Bulirsch's algorithm: the old route through m = 1 - kc**2 lost kc
        below ~1e-8 (cel1(1e-9) returned 2e6; ln(4/kc) = 22.1).  p < 0 is the
        Cauchy principal value = Re Pi(1-p | 1-kc^2) (mpmath)."""
        assert abs(_s(elliptic.cel1(1e-9)) - (22.109560198066302)) < 1e-15 * 22
        assert abs(_s(elliptic.cel1(1e-300)) - (692.1618222593336)) < 1e-15 * 692
        assert abs(_s(elliptic.cel1(2.0)) - (1.0782578237498216)) < 1e-15                      # kc > 1 is m = -3
        assert abs(_s(elliptic.cel(0.5, -0.5, 1.0, 1.0)) - (-1.0782578237498216)) < 1e-15 * 1.1
        assert abs(_s(elliptic.cel(0.7, -5.0, 1.0, 1.0)) - (-0.092277884964284496)) < 1e-15
        assert abs(_s(elliptic.cel(1e-9, 0.3, 1.5, -0.7)) - (-46.045402351061091)) < 1e-15 * 47
        assert _s(elliptic.cel(-0.3, 1.0, 1.0, 1.0)) == _s(elliptic.cel1(0.3))   # depends on kc^2 only
        assert math.isinf(_s(elliptic.cel1(0.0))) and _s(elliptic.cel1(0.0)) > 0
        kc = np.array([1e-12, 0.3, 0.9, 2.5]); p = np.array([-0.4, 0.7, -3.0, 1e-3])
        K = np.asarray(elliptic.cel1(kc))
        assert np.all(np.abs(np.asarray(elliptic.cel(kc, p, 1.0, 0.0)) + p * np.asarray(elliptic.cel(kc, p, 0.0, 1.0)) - K) < 1e-14 * np.maximum(1, K))

    def test_round6c_bdj_delta_jacobiEDJ_complexZ_arclength(self):
        """mpmath at exact doubles: Delta^2 = (1-m) + m cos^2 in ellipticBDJ,
        jacobiEDJ reduces u before the amplitude, complex Z uses the exact
        complete integrals, arclength on mixed arrays."""
        assert abs(_s(elliptic.ellipticBDJ(math.pi/2 - 2e-4, 1 - 1e-8)[1]) - (8.1529993379146678)) < 2e-12 * 9
        Eu, Du, _ = elliptic.jacobiEDJ(1520.3427441800743, 0.999999990458008)
        assert abs(_s(Eu) - 143.00000694616405) < 2e-12 and abs(_s(Du) - 1377.3427503765037) < 2e-12
        Z = _s(elliptic.elliptic12i(complex(-1.053, -0.8215), 1 - 2**-53)[2])
        assert abs(Z - complex(-1.1405612347714637, -0.39945642654606886)) < 1e-14                  # was 1e-11: K taken at double(pi/2)
        v = np.asarray(elliptic.arclength_ellipse(np.array([5, 785.9, 3]), np.array([10, 495.8, 3]), np.array([0, 5.279, 0]), np.array([1, -6.134, 2])))
        assert np.all(np.abs(v - np.array([8.8662512353670695, -7494.1448816975323, 6])) < 1e-14 * np.array([9, 7495, 6]))

    def test_elliptic12i_period_just_below_half_pi(self):
        """pi*ceil(phi/pi - 0.5 + 1e-14) added a period for phi within 3e-14
        below pi/2: Re F came out 3K instead of K (mpmath at the exact double)."""
        F = _s(elliptic.elliptic12i(complex(1.5707963267948961, 0.5), 1/3)[0])
        assert abs(F - complex(1.7339168852579344, 0.62666316872107993)) < 4e-15
        F = _s(elliptic.elliptic12i(complex(math.pi/2 - 1e-12, -2.0), 0.5)[0])
        assert abs(F.real - 0.3901536583) < 1e-9                     # left of the cut, not 2K - ...
