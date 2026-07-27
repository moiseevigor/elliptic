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
