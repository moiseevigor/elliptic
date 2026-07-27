"""Special / boundary value tests not covered by the existing suite.

New coverage
------------
elliptic12  – m=1 closed forms: F=atanh(sin φ), E=sin φ
              φ=π periodicity: F(π|m)=2K(m), E(π|m)=2E(m)
              φ<0 odd symmetry at m=0 and m=1
ellipj      – quarter-period K: sn=1, cn=0, dn=√(1−m), am=π/2
              half-period 2K: sn=0, cn=−1, dn=1
              full-period 4K: sn=0, cn=1, dn=1
              parity: sn odd, cn and dn even
elliptic3   – n=0 degenerate: Π(φ|m,0)=F(φ|m)
              φ↦−φ antisymmetry
              φ=π periodicity: Π(π|m,n)=2Π(π/2|m,n)
ellipticBD  – m→1 limit: B(m)→1, D(m)→∞
              m=0 degenerate: B(0)=D(0)=π/4, K(0)=π/2
ellipticBDJ – φ=π/2 equals complete B/D
              F=B+D and E=B+(1−m)D relations
              odd symmetry in φ
carlsonRC   – RC(0,9)=π/6, RC(0,1/9)=3π/2, RC(1,1)=1
carlsonRF   – RF(0,1,1)=π/2, RF(1,1,1)=1, RF(0,4,4)=π/4
carlsonRD   – RD(1,1,1)=1, RD(0,1,1)=3π/4
carlsonRJ   – RJ(1,1,1,1)=1
agm         – agm(a,a)=a; agm(a,b)=agm(b,a); agm→0 for zero input
nomeq       – nomeq(0)=0, nomeq(1)=1, monotone increasing
cel1        – cel1(1)=π/2
theta       – θ₂(π/2,m)=0; θ₁(π,m)=0; Jacobi derivative product
"""
import numpy as np
import pytest
from scipy.special import ellipk, ellipe

import elliptic


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _f(v):
    """Collapse any array/scalar to a plain Python float."""
    return float(np.asarray(v).flat[0])


# ============================================================================
# elliptic12 — special inputs
# ============================================================================

class TestElliptic12Special:
    @pytest.mark.parametrize("phi", [np.pi / 6, np.pi / 4, np.pi / 3])
    def test_F_m1_atanh(self, phi):
        """F(φ|1) = atanh(sin φ) for |φ| < π/2."""
        F, _, _ = elliptic.elliptic12(phi, 1.0)
        np.testing.assert_allclose(_f(F), np.arctanh(np.sin(phi)), atol=1e-12,
                                   err_msg=f"phi={phi}")

    @pytest.mark.parametrize("phi", [0.0, np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2])
    def test_E_m1_sin(self, phi):
        """E(φ|1) = sin φ for φ ∈ [0, π/2]."""
        _, E, _ = elliptic.elliptic12(phi, 1.0)
        np.testing.assert_allclose(_f(E), np.sin(phi), atol=1e-12,
                                   err_msg=f"phi={phi}")

    @pytest.mark.parametrize("m", [0.2, 0.5, 0.7, 0.9])
    def test_F_phi_pi_is_2K(self, m):
        """F(π|m) = 2K(m)."""
        F, _, _ = elliptic.elliptic12(np.pi, m)
        np.testing.assert_allclose(_f(F), 2.0 * float(ellipk(m)), atol=1e-11,
                                   err_msg=f"m={m}")

    @pytest.mark.parametrize("m", [0.2, 0.5, 0.7, 0.9])
    def test_E_phi_pi_is_2E(self, m):
        """E(π|m) = 2E(m)."""
        _, E, _ = elliptic.elliptic12(np.pi, m)
        np.testing.assert_allclose(_f(E), 2.0 * float(ellipe(m)), atol=1e-11,
                                   err_msg=f"m={m}")

    def test_F_odd_m0(self):
        """F(−φ|0) = −F(φ|0) = −φ  (odd symmetry at m=0)."""
        phis = np.array([0.5, 1.0, 1.5])
        Fp, _, _ = elliptic.elliptic12(phis, 0.0)
        Fm, _, _ = elliptic.elliptic12(-phis, 0.0)
        np.testing.assert_allclose(np.asarray(Fp), -np.asarray(Fm), atol=1e-14)

    def test_F_odd_m1(self):
        """F(−φ|1) = −atanh(sin φ) = −F(φ|1)  (odd symmetry at m=1)."""
        phis = np.array([np.pi / 6, np.pi / 4, np.pi / 3])
        for phi in phis:
            Fp, _, _ = elliptic.elliptic12(phi, 1.0)
            Fm, _, _ = elliptic.elliptic12(-phi, 1.0)
            np.testing.assert_allclose(_f(Fp), -_f(Fm), atol=1e-12,
                                       err_msg=f"phi={phi}")


# ============================================================================
# ellipj — period special values and parity
# ============================================================================

class TestEllipjSpecial:
    @pytest.mark.parametrize("m", [0.1, 0.3, 0.5, 0.7, 0.9])
    def test_quarter_period(self, m):
        """At u=K(m): sn=1, cn=0, dn=√(1−m), am=π/2."""
        K = float(ellipk(m))
        sn, cn, dn, am = elliptic.ellipj(K, m)
        np.testing.assert_allclose(_f(sn), 1.0,            atol=1e-12)
        np.testing.assert_allclose(_f(cn), 0.0,            atol=1e-11)
        np.testing.assert_allclose(_f(dn), np.sqrt(1 - m), atol=1e-12)
        np.testing.assert_allclose(_f(am), np.pi / 2,      atol=1e-11)

    @pytest.mark.parametrize("m", [0.1, 0.3, 0.5, 0.7, 0.9])
    def test_half_period(self, m):
        """At u=2K(m): sn=0, cn=−1, dn=1."""
        K = float(ellipk(m))
        sn, cn, dn, _ = elliptic.ellipj(2 * K, m)
        np.testing.assert_allclose(_f(sn), 0.0,  atol=1e-11)
        np.testing.assert_allclose(_f(cn), -1.0, atol=1e-12)
        np.testing.assert_allclose(_f(dn), 1.0,  atol=1e-12)

    @pytest.mark.parametrize("m", [0.1, 0.5, 0.9])
    def test_full_period(self, m):
        """At u=4K(m): sn=0, cn=1, dn=1  (full period)."""
        K = float(ellipk(m))
        sn, cn, dn, _ = elliptic.ellipj(4 * K, m)
        np.testing.assert_allclose(_f(sn), 0.0, atol=1e-11)
        np.testing.assert_allclose(_f(cn), 1.0, atol=1e-12)
        np.testing.assert_allclose(_f(dn), 1.0, atol=1e-12)

    def test_sn_odd(self):
        """sn(−u|m) = −sn(u|m)."""
        u = np.array([0.3, 0.8, 1.4])
        m = 0.6
        sn_p, _, _, _ = elliptic.ellipj(u, m)
        sn_m, _, _, _ = elliptic.ellipj(-u, m)
        np.testing.assert_allclose(np.asarray(sn_p), -np.asarray(sn_m), atol=1e-14)

    def test_cn_dn_even(self):
        """cn(−u|m) = cn(u|m) and dn(−u|m) = dn(u|m)."""
        u = np.array([0.5, 1.0, 1.5])
        m = 0.7
        _, cn_p, dn_p, _ = elliptic.ellipj(u, m)
        _, cn_m, dn_m, _ = elliptic.ellipj(-u, m)
        np.testing.assert_allclose(np.asarray(cn_p), np.asarray(cn_m), atol=1e-14)
        np.testing.assert_allclose(np.asarray(dn_p), np.asarray(dn_m), atol=1e-14)


# ============================================================================
# elliptic3 — degenerate, antisymmetry, period
# ============================================================================

class TestElliptic3Special:
    @pytest.mark.parametrize("m", [0.2, 0.5, 0.8])
    def test_n0_equals_F(self, m):
        """Π(φ|m,0) = F(φ|m)  (degenerate: no characteristic)."""
        phi = np.array([np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2])
        Pi = np.asarray(elliptic.elliptic3(phi, m, 0.0))
        F, _, _ = elliptic.elliptic12(phi, m)
        np.testing.assert_allclose(Pi, np.asarray(F), atol=1e-12,
                                   err_msg=f"m={m}")

    @pytest.mark.parametrize("n", [0.2, 0.5, 0.8])
    def test_odd_symmetry(self, n):
        """Π(−φ|m,n) = −Π(φ|m,n)."""
        m = 0.5
        phi = np.array([np.pi / 6, np.pi / 4, np.pi / 3])
        for p in phi:
            Pi_p = _f(elliptic.elliptic3(p,  m, n))
            Pi_m = _f(elliptic.elliptic3(-p, m, n))
            np.testing.assert_allclose(Pi_m, -Pi_p, atol=1e-12,
                                       err_msg=f"phi={p} n={n}")

    @pytest.mark.parametrize("m,n", [(0.2, 0.2), (0.5, 0.5), (0.7, 0.2)])
    def test_phi_pi_periodicity(self, m, n):
        """Π(π|m,n) = 2 Π(π/2|m,n)."""
        Pi_pi  = _f(elliptic.elliptic3(np.pi,     m, n))
        Pi_hpi = _f(elliptic.elliptic3(np.pi / 2, m, n))
        np.testing.assert_allclose(Pi_pi, 2.0 * Pi_hpi, atol=1e-7,
                                   err_msg=f"m={m} n={n}")


# ============================================================================
# ellipticBD — complete integrals at boundary moduli
# ============================================================================

class TestEllipticBDSpecial:
    def test_m0_degenerate(self):
        """B(0)=D(0)=π/4;  K(0)=B+D=π/2."""
        B, D, _ = elliptic.ellipticBD(0.0)
        np.testing.assert_allclose(_f(B), np.pi / 4, atol=1e-13)
        np.testing.assert_allclose(_f(D), np.pi / 4, atol=1e-13)

    def test_B_approaches_1_near_m1(self):
        """B(m) → 1 as m → 1  [B=(K−D), E→1, (1−m)K→0 ⟹ D→0 and K→∞]."""
        m_vals = np.array([0.9, 0.99, 0.999, 0.9999])
        B_vals = np.array([_f(elliptic.ellipticBD(m)[0]) for m in m_vals])
        np.testing.assert_allclose(B_vals[-1], 1.0, atol=5e-4)
        assert np.all(np.diff(B_vals) > 0), "B(m) should be increasing toward m=1"

    def test_D_monotone_near_m1(self):
        """D(m) = (K−E)/m grows monotonically as m → 1."""
        m_vals = np.array([0.9, 0.99, 0.999, 0.9999])
        D_vals = np.array([_f(elliptic.ellipticBD(m)[1]) for m in m_vals])
        assert np.all(np.diff(D_vals) > 0), f"D not monotone: {D_vals}"
        assert D_vals[-1] > D_vals[0], f"D(0.9999)={D_vals[-1]} not > D(0.9)={D_vals[0]}"


# ============================================================================
# ellipticBDJ — incomplete integrals: relations and symmetry
# ============================================================================

class TestEllipticBDJSpecial:
    @pytest.mark.parametrize("m", [0.2, 0.5, 0.7, 0.9])
    def test_phi_pi2_equals_complete(self, m):
        """BDJ(π/2|m) equals the complete B(m), D(m)."""
        B_j, D_j, _ = elliptic.ellipticBDJ(np.pi / 2, m)
        B_c, D_c, _ = elliptic.ellipticBD(m)
        np.testing.assert_allclose(_f(B_j), _f(B_c), atol=1e-12, err_msg=f"B m={m}")
        np.testing.assert_allclose(_f(D_j), _f(D_c), atol=1e-12, err_msg=f"D m={m}")

    @pytest.mark.parametrize("m", [0.2, 0.5, 0.7])
    def test_F_equals_B_plus_D(self, m):
        """F(φ|m) = B(φ|m) + D(φ|m) for φ ∈ [0, π/2]."""
        phis = np.array([0.3, 0.6, 0.9, 1.2])
        for phi in phis:
            B, D, _ = elliptic.ellipticBDJ(phi, m)
            F, _, _ = elliptic.elliptic12(phi, m)
            np.testing.assert_allclose(_f(B) + _f(D), _f(F), atol=1e-12,
                                       err_msg=f"phi={phi} m={m}")

    @pytest.mark.parametrize("m", [0.2, 0.5, 0.7])
    def test_E_equals_B_plus_1mm_D(self, m):
        """E(φ|m) = B(φ|m) + (1−m)·D(φ|m) for φ ∈ [0, π/2]."""
        phis = np.array([0.3, 0.6, 0.9, 1.2])
        for phi in phis:
            B, D, _ = elliptic.ellipticBDJ(phi, m)
            _, E, _ = elliptic.elliptic12(phi, m)
            np.testing.assert_allclose(_f(B) + (1 - m) * _f(D), _f(E), atol=1e-12,
                                       err_msg=f"phi={phi} m={m}")

    def test_B_odd_symmetry(self):
        """B(−φ|m) = −B(φ|m) and D(−φ|m) = −D(φ|m)."""
        phis = [0.3, 0.6, 0.9]
        m = 0.5
        for phi in phis:
            B_p, D_p, _ = elliptic.ellipticBDJ(phi, m)
            B_m, D_m, _ = elliptic.ellipticBDJ(-phi, m)
            np.testing.assert_allclose(_f(B_m), -_f(B_p), atol=1e-12)
            np.testing.assert_allclose(_f(D_m), -_f(D_p), atol=1e-12)


# ============================================================================
# carlsonRC — isolated special values
# ============================================================================

class TestCarlsonRCSpecial:
    def test_RC_0_9(self):
        """RC(0, 9) = π/(2·√9) = π/6."""
        RC = elliptic.carlsonRC(0.0, 9.0)
        np.testing.assert_allclose(_f(RC), np.pi / 6, atol=1e-14)

    def test_RC_0_1over9(self):
        """RC(0, 1/9) = π/(2·√(1/9)) = 3π/2."""
        RC = elliptic.carlsonRC(0.0, 1.0 / 9.0)
        np.testing.assert_allclose(_f(RC), 3.0 * np.pi / 2, atol=1e-13)

    def test_RC_1_1(self):
        """RC(1, 1) = 1  [RC(x,x)=1/√x at x=1]."""
        RC = elliptic.carlsonRC(1.0, 1.0)
        np.testing.assert_allclose(_f(RC), 1.0, atol=1e-14)

    @pytest.mark.parametrize("x", [0.25, 1.0, 4.0])
    def test_RC_x_x(self, x):
        """RC(x, x) = 1/√x for any x > 0."""
        RC = elliptic.carlsonRC(x, x)
        np.testing.assert_allclose(_f(RC), 1.0 / np.sqrt(x), atol=1e-13,
                                   err_msg=f"x={x}")


# ============================================================================
# carlsonRF — isolated special values
# ============================================================================

class TestCarlsonRFSpecial:
    def test_RF_0_1_1(self):
        """RF(0, 1, 1) = π/2  [= K(0)]."""
        RF = elliptic.carlsonRF(0.0, 1.0, 1.0)
        np.testing.assert_allclose(_f(RF), np.pi / 2, atol=1e-13)

    def test_RF_1_1_1(self):
        """RF(1, 1, 1) = 1  [RF(x,x,x)=x^{−1/2} at x=1]."""
        RF = elliptic.carlsonRF(1.0, 1.0, 1.0)
        np.testing.assert_allclose(_f(RF), 1.0, atol=1e-13)

    def test_RF_0_4_4(self):
        """RF(0, 4, 4) = π/4  [RC(0,4)=π/(2·2)=π/4]."""
        RF = elliptic.carlsonRF(0.0, 4.0, 4.0)
        np.testing.assert_allclose(_f(RF), np.pi / 4, atol=1e-13)


# ============================================================================
# carlsonRD — isolated special values
# ============================================================================

class TestCarlsonRDSpecial:
    def test_RD_1_1_1(self):
        """RD(1, 1, 1) = 1  [RD(x,x,x)=x^{−3/2} at x=1]."""
        RD = elliptic.carlsonRD(1.0, 1.0, 1.0)
        np.testing.assert_allclose(_f(RD), 1.0, atol=1e-13)

    def test_RD_0_1_1(self):
        """RD(0, 1, 1) = 3π/4  [RD(0,y,y)=3RF(0,y,y)/(2y)=3π/(4√y·y) at y=1]."""
        RD = elliptic.carlsonRD(0.0, 1.0, 1.0)
        np.testing.assert_allclose(_f(RD), 3.0 * np.pi / 4, atol=1e-13)


# ============================================================================
# carlsonRJ — isolated special values
# ============================================================================

class TestCarlsonRJSpecial:
    def test_RJ_1_1_1_1(self):
        """RJ(1, 1, 1, 1) = 1  [RJ(x,x,x,x)=x^{−3/2} at x=1]."""
        RJ = elliptic.carlsonRJ(1.0, 1.0, 1.0, 1.0)
        np.testing.assert_allclose(_f(RJ), 1.0, atol=1e-13)


# ============================================================================
# agm — boundary and symmetry
# ============================================================================

class TestAGMSpecial:
    @pytest.mark.parametrize("a", [0.5, 1.0, 2.0, 5.0])
    def test_equal_inputs(self, a):
        """agm(a, a) = a."""
        result = elliptic.agm(a, a)
        np.testing.assert_allclose(float(np.asarray(result)), a, atol=1e-14)

    @pytest.mark.parametrize("a,b", [(1.0, 0.5), (3.0, 1.0), (0.5, 0.25)])
    def test_symmetry(self, a, b):
        """agm(a, b) = agm(b, a)."""
        r1 = float(np.asarray(elliptic.agm(a, b)))
        r2 = float(np.asarray(elliptic.agm(b, a)))
        np.testing.assert_allclose(r1, r2, atol=1e-14)

    def test_zero_input_small(self):
        """agm(a, 0) ≈ 0  (converges to 0 when one argument is zero)."""
        result = float(np.asarray(elliptic.agm(1.0, 0.0)))
        assert abs(result) < 1e-6, f"agm(1,0)={result}, expected ~0"

    def test_homogeneity(self):
        """agm(λa, λb) = λ agm(a, b)."""
        a, b, lam = 2.0, 1.0, 3.0
        r1 = float(np.asarray(elliptic.agm(lam * a, lam * b)))
        r2 = lam * float(np.asarray(elliptic.agm(a, b)))
        np.testing.assert_allclose(r1, r2, atol=1e-13)


# ============================================================================
# nomeq — boundary values and monotonicity
# ============================================================================

class TestNomeqSpecial:
    def test_m0_is_zero(self):
        """nomeq(0) = 0."""
        q = float(np.asarray(elliptic.nomeq(0.0)))
        np.testing.assert_allclose(q, 0.0, atol=1e-14)

    def test_m1_is_one(self):
        """nomeq(1) = 1  (nome reaches 1 as m→1)."""
        q = float(np.asarray(elliptic.nomeq(1.0)))
        np.testing.assert_allclose(q, 1.0, atol=1e-14)

    def test_monotone_increasing(self):
        """nomeq(m) is strictly increasing in m."""
        m = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        q = np.asarray(elliptic.nomeq(m))
        assert np.all(np.diff(q) > 0), "nomeq must be monotone increasing"

    def test_range(self):
        """nomeq(m) ∈ (0, 1) for m ∈ (0, 1)."""
        m = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        q = np.asarray(elliptic.nomeq(m))
        assert np.all(q > 0) and np.all(q < 1)


# ============================================================================
# cel1 — degenerate case
# ============================================================================

class TestCel1Special:
    def test_kc_1_is_pi_over_2(self):
        """cel1(1) = π/2  [cel1(kc)=K(m) at kc=1 corresponds to K(0)=π/2]."""
        c = float(np.asarray(elliptic.cel1(1.0)))
        np.testing.assert_allclose(c, np.pi / 2, atol=1e-13)

    def test_array(self):
        """cel1 on an array matches individual scalar calls."""
        kc_vals = np.array([0.5, 0.8, 1.0])
        c_arr = np.asarray(elliptic.cel1(kc_vals))
        for i, kc in enumerate(kc_vals):
            c_s = float(np.asarray(elliptic.cel1(kc)))
            np.testing.assert_allclose(c_arr[i], c_s, atol=1e-14)


# ============================================================================
# theta — zeros and Jacobi derivative product
# ============================================================================

class TestThetaSpecial:
    @pytest.mark.parametrize("m", [0.1, 0.3, 0.5, 0.7, 0.9])
    def test_theta2_zero_at_pi_over_2(self, m):
        """θ₂(π/2, m) = 0  (all terms contain cos((2n+1)·π/2) = 0)."""
        th2 = float(np.asarray(elliptic.theta(2, np.pi / 2, m)))
        np.testing.assert_allclose(th2, 0.0, atol=1e-13, err_msg=f"m={m}")

    @pytest.mark.parametrize("m", [0.1, 0.3, 0.5, 0.7, 0.9])
    def test_theta1_zero_at_pi(self, m):
        """θ₁(π, m) = 0  (quasi-period: θ₁(z+π) = −θ₁(z), so θ₁(π)=−θ₁(0)=0)."""
        th1 = float(np.asarray(elliptic.theta(1, np.pi, m)))
        np.testing.assert_allclose(th1, 0.0, atol=1e-13, err_msg=f"m={m}")

    @pytest.mark.parametrize("m", [0.1, 0.3, 0.5, 0.7, 0.9])
    def test_jacobi_derivative_product(self, m):
        """θ₁′(0, m) = θ₂(0, m) · θ₃(0, m) · θ₄(0, m)."""
        _, dth1 = elliptic.theta_prime(1, 0.0, m)
        th2, _  = elliptic.theta_prime(2, 0.0, m)
        th3, _  = elliptic.theta_prime(3, 0.0, m)
        th4, _  = elliptic.theta_prime(4, 0.0, m)
        product = _f(th2) * _f(th3) * _f(th4)
        np.testing.assert_allclose(_f(dth1), product, rtol=1e-10, err_msg=f"m={m}")
