"""Tests for functions added to reach feature parity with the MATLAB package:
agm, nomeq, inversenomeq, inverselliptic2, jacobiThetaEta, theta, theta_prime,
elliptic12i, ellipji, weierstrassInvariants, weierstrassPPrime, arclength_ellipse.
"""
import numpy as np
import pytest
from scipy.special import ellipk, ellipe, ellipkinc, ellipeinc, ellipj as sp_ellipj

import elliptic


# -----------------------------------------------------------------------
# agm
# -----------------------------------------------------------------------
class TestAGM:
    def test_known_value(self):
        """agm(1, 1/√2) = π/(2K(1/2)) — A&S 17.6."""
        from scipy.special import ellipk
        # K(1/2) = agm(1, 1/√2) * π / (2 * agm(1, 1/√2))  — trivially
        # Direct: agm(1, cos(π/4)) = π / (2 K(m=sin²(π/4)=0.5)) is wrong
        # Correct: K(m) = π / (2 * agm(1, sqrt(1-m)))
        m = 0.7
        K_ref = float(ellipk(m))
        agm_val = float(elliptic.agm(1.0, np.sqrt(1.0 - m)))
        assert abs(np.pi / (2.0 * agm_val) - K_ref) < 1e-13

    def test_equal_inputs(self):
        assert abs(float(elliptic.agm(3.0, 3.0)) - 3.0) < 1e-14

    def test_array(self):
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([1.0, 1.0, 1.0])
        result = elliptic.agm(a, b)
        for ai, bi, ri in zip(a, b, np.asarray(result)):
            expected = float(elliptic.agm(ai, bi))
            assert abs(ri - expected) < 1e-14


# -----------------------------------------------------------------------
# nomeq / inversenomeq
# -----------------------------------------------------------------------
class TestNome:
    def test_nomeq_m0(self):
        assert abs(float(elliptic.nomeq(0.0))) < 1e-15

    def test_nomeq_array(self):
        m = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        q = np.asarray(elliptic.nomeq(m))
        # q = exp(-π K'/K), all in (0, 1)
        assert np.all(q > 0) and np.all(q < 1)
        # monotone increasing with m
        assert np.all(np.diff(q) > 0)

    def test_roundtrip(self):
        """inversenomeq(nomeq(m)) ≈ m for m in [0.1, 0.7]."""
        m = np.array([0.1, 0.3, 0.5, 0.7])
        q = np.asarray(elliptic.nomeq(m))
        m2 = np.asarray(elliptic.inversenomeq(q))
        np.testing.assert_allclose(m2, m, atol=1e-10)

    def test_inversenomeq_q0(self):
        assert abs(float(elliptic.inversenomeq(0.0))) < 1e-12


# -----------------------------------------------------------------------
# inverselliptic2
# -----------------------------------------------------------------------
class TestInverselliptic2:
    def test_roundtrip_scalar(self):
        """E(inverselliptic2(E0, m), m) ≈ E0."""
        m = 0.5
        phi0 = np.pi / 3
        _, E0, _ = elliptic.elliptic12(phi0, m)
        phi_inv = float(elliptic.inverselliptic2(float(E0), m))
        assert abs(phi_inv - phi0) < 1e-12

    def test_roundtrip_array(self):
        phi = np.linspace(0.1, np.pi / 2, 10)
        m   = 0.6
        _, E_vals, _ = elliptic.elliptic12(phi, m)
        phi_inv = np.asarray(elliptic.inverselliptic2(E_vals, m))
        np.testing.assert_allclose(phi_inv, phi, atol=1e-10)

    def test_phi0(self):
        assert abs(float(elliptic.inverselliptic2(0.0, 0.5))) < 1e-13


# -----------------------------------------------------------------------
# jacobiThetaEta
# -----------------------------------------------------------------------
class TestJacobiThetaEta:
    def test_Theta_m0(self):
        """At m=0, Θ = 1."""
        Th, H = elliptic.jacobiThetaEta(np.array([0.5, 1.0]), 0.0)
        np.testing.assert_allclose(np.asarray(Th), 1.0, atol=1e-14)

    def test_H_m0(self):
        """At m=0, H = 0 (nome q=0, all series terms vanish)."""
        Th, H = elliptic.jacobiThetaEta(np.array([0.5, 1.0]), 0.0)
        np.testing.assert_allclose(np.asarray(H), 0.0, atol=1e-14)

    def test_Theta_m1_is_nan(self):
        Th, H = elliptic.jacobiThetaEta(np.array([0.5]), 1.0)
        assert np.isnan(float(Th[0]))

    def test_Theta_positive(self):
        """Θ(u,m) > 0 for 0 < m < 1."""
        u = np.linspace(0.1, 2.0, 10)
        Th, _ = elliptic.jacobiThetaEta(u, 0.5)
        assert np.all(np.asarray(Th) > 0)

    def test_sn_from_theta(self):
        """sn(u|m) = (√k)⁻¹ · H(u) / Θ(u)  where k = m^{1/4}   (A&S 16.27)."""
        m = 0.6
        u = np.linspace(0.1, 1.0, 8)
        Th, H = elliptic.jacobiThetaEta(u, m)
        sn_theta = np.asarray(H) / (m**0.25 * np.asarray(Th))
        sn_ref, _, _, _ = elliptic.ellipj(u, m)
        np.testing.assert_allclose(sn_theta, np.asarray(sn_ref), atol=1e-10)


# -----------------------------------------------------------------------
# theta (four types)
# -----------------------------------------------------------------------
class TestTheta:
    def test_th3_m0(self):
        """θ₃(v, q=0) = 1."""
        v = np.array([0.3, 1.0, 2.0])
        Th = elliptic.theta(3, v, 0.0)
        np.testing.assert_allclose(np.asarray(Th), 1.0, atol=1e-14)

    def test_th4_m0(self):
        Th = elliptic.theta(4, np.array([0.5, 1.5]), 0.0)
        np.testing.assert_allclose(np.asarray(Th), 1.0, atol=1e-14)

    def test_th1_odd(self):
        """θ₁(-v) = -θ₁(v)."""
        v = np.array([0.3, 0.8, 1.2])
        m = 0.5
        Tp = np.asarray(elliptic.theta(1,  v, m))
        Tm = np.asarray(elliptic.theta(1, -v, m))
        np.testing.assert_allclose(Tp, -Tm, atol=1e-13)

    def test_th2_even(self):
        """θ₂(-v) = θ₂(v)."""
        v = np.array([0.3, 0.8])
        m = 0.5
        np.testing.assert_allclose(
            np.asarray(elliptic.theta(2, v, m)),
            np.asarray(elliptic.theta(2, -v, m)),
            atol=1e-13)

    def test_th3_th4_relation(self):
        """θ₃(v + π/2) = θ₄(v)  (A&S 16.39)."""
        v = np.array([0.2, 0.7, 1.1])
        m = 0.6
        Th3_shifted = np.asarray(elliptic.theta(3, v + np.pi/2, m))
        Th4         = np.asarray(elliptic.theta(4, v, m))
        np.testing.assert_allclose(Th3_shifted, Th4, atol=1e-12)

    def test_th1_th2_relation(self):
        """θ₁(v + π/2) = θ₂(v)  (A&S 16.39)."""
        v = np.array([0.2, 0.7])
        m = 0.5
        Th1_shifted = np.asarray(elliptic.theta(1, v + np.pi/2, m))
        Th2         = np.asarray(elliptic.theta(2, v, m))
        np.testing.assert_allclose(Th1_shifted, Th2, atol=1e-12)

    def test_jacobi_consistency(self):
        """theta and jacobiThetaEta agree via the u ↔ v mapping."""
        m = 0.5
        K = float(ellipk(m))
        v = np.array([0.3, 0.7, 1.1])
        u = 2.0 * K * v / np.pi
        Th_jac, _ = elliptic.jacobiThetaEta(u, m)
        Th4       = np.asarray(elliptic.theta(4, v, m))
        np.testing.assert_allclose(np.asarray(Th_jac), Th4, atol=1e-12)


# -----------------------------------------------------------------------
# theta_prime
# -----------------------------------------------------------------------
class TestThetaPrime:
    @pytest.mark.parametrize("j", [1, 2, 3, 4])
    def test_derivative_numerical(self, j):
        """theta_prime vs central finite difference for all four types."""
        v0, m, h = 0.6, 0.5, 1e-7
        _, thp = elliptic.theta_prime(j, np.array([v0]), m)
        th_p = float(np.asarray(elliptic.theta(j, np.array([v0 + h]), m))[0])
        th_m = float(np.asarray(elliptic.theta(j, np.array([v0 - h]), m))[0])
        thp_num = (th_p - th_m) / (2 * h)
        assert abs(float(np.asarray(thp)[0]) - thp_num) < 1e-7, \
            f"theta_prime({j}) numerical mismatch: {float(np.asarray(thp)[0])} vs {thp_num}"


# -----------------------------------------------------------------------
# elliptic12i (complex F, E, Z)
# -----------------------------------------------------------------------
class TestElliptic12i:
    def test_real_case(self):
        """For purely real phi, elliptic12i and elliptic12 agree."""
        phi = np.array([0.3, 0.7, 1.2])
        m   = 0.5
        Fi, Ei, Zi = elliptic.elliptic12i(phi.astype(complex), m)
        F,  E,  Z  = elliptic.elliptic12(phi, m)
        np.testing.assert_allclose(np.real(np.asarray(Fi)), np.asarray(F), atol=1e-12)
        np.testing.assert_allclose(np.real(np.asarray(Ei)), np.asarray(E), atol=1e-12)

    def test_purely_imaginary_is_imaginary(self):
        """F(iψ|m) is purely imaginary (real part = 0) for ψ real."""
        psi, m = 0.6, 0.5
        Fi, _, _ = elliptic.elliptic12i(np.array([1j * psi]), m)
        assert abs(np.real(complex(np.asarray(Fi).flat[0]))) < 1e-12

    def test_z_definition(self):
        """Zi = Ei - (E(m)/K(m)) * Fi  [definition of Jacobi zeta]."""
        u  = np.array([0.4 + 0.2j, 0.8 + 0.5j])
        m  = 0.6
        Fi, Ei, Zi = elliptic.elliptic12i(u, m)
        K_m = float(ellipk(m)); E_m = float(ellipe(m))
        Zi_ref = np.asarray(Ei) - (E_m / K_m) * np.asarray(Fi)
        np.testing.assert_allclose(np.asarray(Zi), Zi_ref, atol=1e-12)


# -----------------------------------------------------------------------
# ellipji (complex Jacobi functions)
# -----------------------------------------------------------------------
class TestEllipji:
    def test_real_case(self):
        """For real u, ellipji agrees with ellipj."""
        u = np.array([0.5, 1.0, 1.5])
        m = 0.7
        sni, cni, dni = elliptic.ellipji(u.astype(complex), m)
        sn, cn, dn, _ = elliptic.ellipj(u, m)
        np.testing.assert_allclose(np.real(np.asarray(sni)), np.asarray(sn), atol=1e-12)
        np.testing.assert_allclose(np.real(np.asarray(cni)), np.asarray(cn), atol=1e-12)
        np.testing.assert_allclose(np.real(np.asarray(dni)), np.asarray(dn), atol=1e-12)

    def test_identity_sn2_cn2(self):
        """sn² + cn² = 1 for complex argument."""
        u = np.array([0.3 + 0.4j, 0.8 + 0.2j, 1.0 - 0.3j])
        m = 0.5
        sni, cni, _ = elliptic.ellipji(u, m)
        identity = np.asarray(sni)**2 + np.asarray(cni)**2
        np.testing.assert_allclose(identity, 1.0 + 0j, atol=1e-12)

    def test_identity_dn2(self):
        """dn² + m sn² = 1 for complex argument."""
        u = np.array([0.5 + 0.3j, 1.2 + 0.1j])
        m = 0.6
        sni, _, dni = elliptic.ellipji(u, m)
        identity = np.asarray(dni)**2 + m * np.asarray(sni)**2
        np.testing.assert_allclose(identity, 1.0 + 0j, atol=1e-12)


# -----------------------------------------------------------------------
# weierstrassInvariants
# -----------------------------------------------------------------------
class TestWeierstrassInvariants:
    def test_known_values(self):
        """g2, g3, Delta for e1=2, e2=0, e3=-2 (A&S example)."""
        g2, g3, D = elliptic.weierstrassInvariants(2.0, 0.0, -2.0)
        # g2 = -4*(0 - 4 + 0) = 16
        assert abs(float(g2) - 16.0) < 1e-14
        # g3 = 4*2*0*(-2) = 0
        assert abs(float(g3)) < 1e-14
        # Delta = 16^3 - 0 = 4096
        assert abs(float(D) - 4096.0) < 1e-12

    def test_sum_zero(self):
        """g3 = 0 whenever e2 = 0 and e1 = -e3."""
        e1, e3 = 3.0, -3.0
        _, g3, _ = elliptic.weierstrassInvariants(e1, 0.0, e3)
        assert abs(float(g3)) < 1e-14

    def test_P_differential_equation(self):
        """(℘')² = 4℘³ - g₂℘ - g₃  at a non-singular point."""
        e1, e2, e3 = 2.0, 0.5, -2.5
        z = 0.4
        g2, g3, _ = elliptic.weierstrassInvariants(e1, e2, e3)
        P  = float(elliptic.weierstrassP(z, e1, e2, e3))
        dP = float(elliptic.weierstrassPPrime(z, e1, e2, e3))
        lhs = dP**2
        rhs = 4*P**3 - float(g2)*P - float(g3)
        assert abs(lhs - rhs) < 1e-8, f"Weierstrass ODE failed: {lhs} vs {rhs}"


# -----------------------------------------------------------------------
# weierstrassPPrime
# -----------------------------------------------------------------------
class TestWeierstrassPPrime:
    def test_pole(self):
        """℘'(0) = ∞."""
        dP = elliptic.weierstrassPPrime(0.0, 2.0, 0.0, -2.0)
        assert np.isinf(float(dP))

    def test_numerical_derivative(self):
        """℘'(z) ≈ (℘(z+h) - ℘(z-h)) / (2h)."""
        e1, e2, e3 = 2.0, 0.5, -2.5
        z, h = 0.4, 1e-6
        dP_analytic = float(elliptic.weierstrassPPrime(z, e1, e2, e3))
        Pp = float(elliptic.weierstrassP(z + h, e1, e2, e3))
        Pm = float(elliptic.weierstrassP(z - h, e1, e2, e3))
        dP_num = (Pp - Pm) / (2 * h)
        assert abs(dP_analytic - dP_num) < 1e-7

    def test_odd(self):
        """℘'(-z) = -℘'(z)."""
        e1, e2, e3 = 3.0, 1.0, -4.0
        z = 0.3
        assert abs(float(elliptic.weierstrassPPrime(z, e1, e2, e3)) +
                   float(elliptic.weierstrassPPrime(-z, e1, e2, e3))) < 1e-12


# -----------------------------------------------------------------------
# arclength_ellipse
# -----------------------------------------------------------------------
class TestArclengthEllipse:
    def test_full_perimeter_mathematica(self):
        """Full perimeter a=5,b=10 matches Mathematica to 1e-8."""
        arc = elliptic.arclength_ellipse(5, 10)
        assert abs(arc - 48.442241102738385905) < 1e-8

    def test_arc_segment_mathematica(self):
        """Arc from pi/10 to pi/2, a=5,b=10 matches Mathematica."""
        arc = elliptic.arclength_ellipse(5, 10, np.pi / 10, np.pi / 2)
        assert abs(arc - 9.0073890635296310688) < 1e-8

    def test_circle(self):
        """Circle: arclength = r * (theta1 - theta0)."""
        assert abs(elliptic.arclength_ellipse(3.0, 3.0, 0, np.pi) - 3.0 * np.pi) < 1e-13

    def test_antisymmetry(self):
        """Reversing limits negates arc length."""
        a, b = 2.0, 5.0
        arc1 = elliptic.arclength_ellipse(a, b, 0.2, 0.8)
        arc2 = elliptic.arclength_ellipse(a, b, 0.8, 0.2)
        assert abs(arc1 + arc2) < 1e-12

    def test_a_gt_b(self):
        """Prolate case (a > b) returns positive arc length."""
        arc = elliptic.arclength_ellipse(10, 5, 0, np.pi / 2)
        assert arc > 0
