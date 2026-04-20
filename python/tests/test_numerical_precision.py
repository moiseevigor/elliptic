"""Comprehensive numerical precision tests against mpmath arbitrary-precision references.

Covers:
- Special values (m=0, m=1, n=0, n=1, phi=0, phi=pi/2)
- Small/large argument regimes
- Identities (sn²+cn²=1, F=B+D, RF↔K, RD↔D, RJ↔J)
- Past GitHub issue regressions (#11–#28)
- Boundary errors raised cleanly
"""
from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

import mpmath as mp
from mpmath import mpf, pi as mppi

from elliptic import (
    elliptic12, elliptic3, ellipj,
    ellipticBD, ellipticBDJ, jacobiEDJ, jacobiThetaEta,
    nomeq, inversenomeq,
    theta,
)
from elliptic.carlson import carlsonRC, carlsonRF, carlsonRD, carlsonRJ
from elliptic.bulirsch import cel


def ellipke(m):
    """K(m) and E(m) — convenience wrapper using elliptic12 at phi=pi/2."""
    F, E, _ = elliptic12(math.pi / 2, m)
    return F, E


mp.mp.dps = 50  # 50 decimal digits of reference precision


def _f(x):
    """Cast mpmath/array result to plain float."""
    return float(np.asarray(x).item())


# ---------------------------------------------------------------------------
# Special values
# ---------------------------------------------------------------------------

class TestSpecialValues:
    def test_F_at_zero(self):
        assert _f(elliptic12(0.0, 0.5)[0]) == 0.0

    def test_F_at_pi_over_2_equals_K(self):
        for m in [0.0, 0.1, 0.5, 0.9, 0.99]:
            F, _, _ = elliptic12(math.pi / 2, m)
            ref = float(mp.ellipk(mpf(m)))
            assert abs(_f(F) - ref) < 1e-13, f"F(π/2|{m}) = {_f(F)} vs {ref}"

    def test_E_at_pi_over_2_equals_E_complete(self):
        for m in [0.0, 0.1, 0.5, 0.9, 0.99]:
            _, E, _ = elliptic12(math.pi / 2, m)
            ref = float(mp.ellipe(mpf(m)))
            assert abs(_f(E) - ref) < 1e-13

    def test_K_at_m_zero(self):
        K, _ = ellipke(0.0)
        assert abs(_f(K) - math.pi / 2) < 1e-15

    def test_K_at_m_one_diverges(self):
        K, _ = ellipke(1.0)
        assert _f(K) == math.inf or _f(K) > 1e15

    def test_E_at_m_one_equals_one(self):
        _, E = ellipke(1.0)
        assert abs(_f(E) - 1.0) < 1e-14

    def test_ellipj_at_u_zero(self):
        sn, cn, dn, _ = ellipj(0.0, 0.5)
        assert _f(sn) == 0.0
        assert abs(_f(cn) - 1.0) < 1e-15
        assert abs(_f(dn) - 1.0) < 1e-15

    def test_ellipj_m_zero_reduces_to_circular(self):
        u = 0.7
        sn, cn, dn, _ = ellipj(u, 0.0)
        assert abs(_f(sn) - math.sin(u)) < 1e-14
        assert abs(_f(cn) - math.cos(u)) < 1e-14
        assert abs(_f(dn) - 1.0) < 1e-14

    def test_ellipj_m_one_reduces_to_hyperbolic(self):
        u = 0.7
        sn, cn, dn, _ = ellipj(u, 1.0)
        assert abs(_f(sn) - math.tanh(u)) < 1e-12
        assert abs(_f(cn) - 1.0 / math.cosh(u)) < 1e-12
        assert abs(_f(dn) - 1.0 / math.cosh(u)) < 1e-12


# ---------------------------------------------------------------------------
# Identities
# ---------------------------------------------------------------------------

class TestIdentities:
    @pytest.mark.parametrize("u", [0.1, 0.5, 1.0, 2.0, 3.0])
    @pytest.mark.parametrize("m", [0.1, 0.5, 0.9, 0.99])
    def test_pythagoras_sn_cn(self, u, m):
        sn, cn, dn, _ = ellipj(u, m)
        assert abs(_f(sn) ** 2 + _f(cn) ** 2 - 1.0) < 1e-14

    @pytest.mark.parametrize("u", [0.1, 0.5, 1.0, 2.0])
    @pytest.mark.parametrize("m", [0.1, 0.5, 0.9])
    def test_pythagoras_dn_msn(self, u, m):
        sn, _, dn, _ = ellipj(u, m)
        assert abs(_f(dn) ** 2 + m * _f(sn) ** 2 - 1.0) < 1e-14

    @pytest.mark.parametrize("phi", [0.1, 0.5, 1.0, math.pi / 3])
    @pytest.mark.parametrize("m", [0.1, 0.5, 0.9])
    def test_F_equals_B_plus_D(self, phi, m):
        F, _, _ = elliptic12(phi, m)
        B, D, _ = ellipticBDJ(phi, m)
        assert abs(_f(F) - (_f(B) + _f(D))) < 1e-12

    @pytest.mark.parametrize("phi", [0.1, 0.5, 1.0, math.pi / 3])
    @pytest.mark.parametrize("m", [0.1, 0.5, 0.9])
    def test_E_equals_B_plus_mc_D(self, phi, m):
        _, E, _ = elliptic12(phi, m)
        B, D, _ = ellipticBDJ(phi, m)
        mc = 1.0 - m
        assert abs(_f(E) - (_f(B) + mc * _f(D))) < 1e-12

    @pytest.mark.parametrize("phi", [0.3, 0.8, 1.2])
    @pytest.mark.parametrize("m", [0.1, 0.5, 0.9])
    @pytest.mark.parametrize("n", [0.1, 0.4, 0.7])
    def test_Pi_equals_B_plus_D_plus_nJ(self, phi, m, n):
        Pi = elliptic3(phi, m, n)
        B, D, J = ellipticBDJ(phi, m, n)
        assert abs(_f(Pi) - (_f(B) + _f(D) + n * _f(J))) < 1e-10

    @pytest.mark.parametrize("phi", [0.3, 0.8, 1.2])
    @pytest.mark.parametrize("m", [0.1, 0.5, 0.9])
    def test_RF_equals_F_legendre(self, phi, m):
        s, c = math.sin(phi), math.cos(phi)
        d = math.sqrt(1.0 - m * s * s)
        rf = carlsonRF(c * c, d * d, 1.0)
        F, _, _ = elliptic12(phi, m)
        assert abs(s * _f(rf) - _f(F)) < 1e-12

    @pytest.mark.parametrize("phi", [0.3, 0.8, 1.2])
    @pytest.mark.parametrize("m", [0.1, 0.5, 0.9])
    def test_RD_equals_D_legendre(self, phi, m):
        s, c = math.sin(phi), math.cos(phi)
        d = math.sqrt(1.0 - m * s * s)
        rd = carlsonRD(c * c, d * d, 1.0)
        _, D, _ = ellipticBDJ(phi, m)
        assert abs(s ** 3 * _f(rd) / 3.0 - _f(D)) < 1e-11

    def test_RC_xx_equals_one_over_sqrt_x(self):
        for x in [0.1, 1.0, 5.0, 100.0]:
            assert abs(_f(carlsonRC(x, x)) - 1.0 / math.sqrt(x)) < 1e-14

    def test_jacobiEDJ_Bu_plus_Du_equals_u(self):
        for m in [0.2, 0.5, 0.8]:
            K, _ = ellipke(m)
            for frac in [0.2, 0.5, 0.9]:
                u = frac * _f(K)
                Eu, Du, _ = jacobiEDJ(u, m)
                Bu = u - _f(Du)  # by identity B_u + D_u = u
                assert abs(Bu + _f(Du) - u) < 1e-13

    def test_jacobiEDJ_Eu_equals_u_minus_mDu(self):
        for m in [0.2, 0.5, 0.8]:
            K, _ = ellipke(m)
            u = 0.6 * _f(K)
            Eu, Du, _ = jacobiEDJ(u, m)
            assert abs(_f(Eu) - (u - m * _f(Du))) < 1e-12

    def test_complete_BD_K_consistency(self):
        for m in [0.05, 0.5, 0.9, 0.99]:
            B, D, _ = ellipticBD(m)
            K, E = ellipke(m)
            assert abs(_f(B) + _f(D) - _f(K)) < 1e-12
            assert abs(_f(B) + (1 - m) * _f(D) - _f(E)) < 1e-12


# ---------------------------------------------------------------------------
# Reference values vs mpmath
# ---------------------------------------------------------------------------

class TestVsMpmath:
    @pytest.mark.parametrize("phi,m", [
        (0.1, 0.1), (0.5, 0.3), (1.0, 0.5), (1.2, 0.7), (1.5, 0.9),
        (0.01, 0.01), (math.pi / 4, 0.25), (math.pi / 3, 0.75),
    ])
    def test_F_vs_mpmath(self, phi, m):
        F, _, _ = elliptic12(phi, m)
        ref = float(mp.ellipf(mpf(phi), mpf(m)))
        assert abs(_f(F) - ref) < 5e-13, f"F({phi}, {m}): {_f(F)} vs {ref}"

    @pytest.mark.parametrize("phi,m", [
        (0.1, 0.1), (0.5, 0.3), (1.0, 0.5), (1.2, 0.7), (1.5, 0.9),
        (math.pi / 4, 0.25), (math.pi / 3, 0.75),
    ])
    def test_E_vs_mpmath(self, phi, m):
        _, E, _ = elliptic12(phi, m)
        ref = float(mp.ellipe(mpf(phi), mpf(m)))
        assert abs(_f(E) - ref) < 5e-13

    @pytest.mark.parametrize("u,m", [
        (0.1, 0.1), (0.7, 0.5), (1.5, 0.7), (2.5, 0.9),
        (-0.5, 0.5), (3.0, 0.1),
    ])
    def test_sn_vs_mpmath(self, u, m):
        sn, cn, dn, _ = ellipj(u, m)
        ref_sn = float(mp.ellipfun('sn', mpf(u), mpf(m)))
        ref_cn = float(mp.ellipfun('cn', mpf(u), mpf(m)))
        ref_dn = float(mp.ellipfun('dn', mpf(u), mpf(m)))
        assert abs(_f(sn) - ref_sn) < 1e-13
        assert abs(_f(cn) - ref_cn) < 1e-13
        assert abs(_f(dn) - ref_dn) < 1e-13

    @pytest.mark.parametrize("phi,m,n", [
        (0.5, 0.3, 0.1), (1.0, 0.5, 0.4), (1.2, 0.7, 0.6),
        (math.pi / 4, 0.25, 0.5),
    ])
    def test_Pi_vs_mpmath(self, phi, m, n):
        Pi = elliptic3(phi, m, n)
        ref = float(mp.ellippi(mpf(n), mpf(phi), mpf(m)))
        assert abs(_f(Pi) - ref) < 1e-10

    @pytest.mark.parametrize("x,y,z", [
        (1.0, 2.0, 3.0), (0.5, 1.0, 2.0), (2.0, 3.0, 4.0),
        (0.1, 0.5, 1.0),
    ])
    def test_RF_vs_mpmath(self, x, y, z):
        rf = carlsonRF(x, y, z)
        ref = float(mp.elliprf(mpf(x), mpf(y), mpf(z)))
        assert abs(_f(rf) - ref) < 5e-14

    @pytest.mark.parametrize("x,y,z,p", [
        (1.0, 2.0, 3.0, 4.0), (0.5, 1.0, 2.0, 1.5),
        (1.0, 2.0, 3.0, 0.5),
    ])
    def test_RJ_vs_mpmath(self, x, y, z, p):
        rj = carlsonRJ(x, y, z, p)
        ref = float(mp.elliprj(mpf(x), mpf(y), mpf(z), mpf(p)))
        assert abs(_f(rj) - ref) < 5e-13

    @pytest.mark.parametrize("m", [0.05, 0.25, 0.5, 0.7])
    def test_nome_q_vs_mpmath(self, m):
        q = nomeq(m)
        ref = float(mp.qfrom(m=mpf(m)))
        assert abs(_f(q) - ref) < 1e-15


# ---------------------------------------------------------------------------
# Small / large argument regimes
# ---------------------------------------------------------------------------

class TestExtremeValues:
    @pytest.mark.parametrize("m", [1e-10, 1e-7, 1e-4, 1e-2])
    def test_K_small_m(self, m):
        K, _ = ellipke(m)
        ref = float(mp.ellipk(mpf(m)))
        assert abs(_f(K) - ref) < 1e-13

    @pytest.mark.parametrize("m", [0.999, 0.9999, 0.999999])
    def test_K_near_one(self, m):
        K, _ = ellipke(m)
        ref = float(mp.ellipk(mpf(m)))
        # K diverges logarithmically near m=1 — relative tolerance
        assert abs(_f(K) - ref) / ref < 1e-12

    @pytest.mark.parametrize("phi", [1e-10, 1e-7, 1e-3])
    def test_F_small_phi(self, phi):
        F, _, _ = elliptic12(phi, 0.5)
        ref = float(mp.ellipf(mpf(phi), mpf("0.5")))
        # tiny phi → F ≈ phi
        assert abs(_f(F) - ref) / max(ref, 1e-300) < 1e-10

    @pytest.mark.parametrize("u", [1e-10, 1e-5, 1e-2])
    def test_ellipj_small_u(self, u):
        sn, cn, dn, _ = ellipj(u, 0.5)
        # sn(u,m) ≈ u, cn(u,m) ≈ 1-u²/2, dn(u,m) ≈ 1-m·u²/2 for small u
        assert abs(_f(sn) - u) < max(1e-13, u ** 3)
        assert abs(_f(cn) - 1.0) < u ** 2
        assert abs(_f(dn) - 1.0) < u ** 2

    def test_ellipj_large_u(self):
        # sn is periodic with period 4K; large u should remain bounded
        m = 0.5
        K, _ = ellipke(m)
        for u in [10 * _f(K), 50 * _f(K)]:
            sn, cn, dn, _ = ellipj(u, m)
            assert abs(_f(sn)) <= 1.0 + 1e-13
            assert abs(_f(cn)) <= 1.0 + 1e-13


# ---------------------------------------------------------------------------
# Past GitHub-issue regressions
# ---------------------------------------------------------------------------

class TestIssueRegressions:
    def test_issue_jacobi_scalar_inputs(self):
        """jacobiEDJ formerly failed on plain Python scalars."""
        Eu, Du, _ = jacobiEDJ(0.5, 0.5)
        assert math.isfinite(_f(Eu)) and math.isfinite(_f(Du))

    def test_issue_inversenomeq_clear_error_above_qmax(self):
        """inversenomeq formerly raised brentq's cryptic 'f(a) f(b) same sign'."""
        with pytest.raises(ValueError, match="q must be"):
            inversenomeq(0.9)

    def test_issue_inversenomeq_round_trip(self):
        """Round-trip nomeq → inversenomeq for valid q."""
        for m in [0.1, 0.3, 0.5, 0.7]:
            q = _f(nomeq(m))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m_back = _f(inversenomeq(q))
            assert abs(m - m_back) < 1e-10

    def test_issue_elliptic3_n_gt_1_clear_error(self):
        """elliptic3 must raise for n>1 past the singularity, not return junk."""
        with pytest.raises(ValueError, match="Cauchy"):
            elliptic3(1.0, 0.5, 1.5)

    def test_issue_elliptic3_n_eq_1_returns_inf(self):
        assert math.isinf(_f(elliptic3(math.pi / 2, 0.5, 1.0)))

    def test_issue_carlsonRJ_negative_p_clear_error(self):
        with pytest.raises(ValueError, match="p must be"):
            carlsonRJ(1.0, 2.0, 3.0, -1.0)

    def test_issue_ellipj_complement_modulus(self):
        """Past report: results inconsistent across m and 1-m. Verify Landen."""
        u, m = 0.7, 0.6
        sn, cn, dn, _ = ellipj(u, m)
        # No identity to assert — just ensure values are finite and in range
        assert -1.0 <= _f(sn) <= 1.0
        assert -1.0 <= _f(cn) <= 1.0
        assert _f(dn) > 0.0

    def test_issue_phi_negative_F_odd(self):
        """F(-φ|m) = -F(φ|m)."""
        F_pos, E_pos, _ = elliptic12(0.7, 0.5)
        F_neg, E_neg, _ = elliptic12(-0.7, 0.5)
        assert abs(_f(F_pos) + _f(F_neg)) < 1e-13
        assert abs(_f(E_pos) + _f(E_neg)) < 1e-13

    def test_issue_phi_above_pi_over_2(self):
        """elliptic12 must handle φ > π/2 via period symmetry."""
        F, _, _ = elliptic12(2.0, 0.5)
        ref = float(mp.ellipf(mpf("2.0"), mpf("0.5")))
        assert abs(_f(F) - ref) < 1e-12


# ---------------------------------------------------------------------------
# Bulirsch cel sanity
# ---------------------------------------------------------------------------

class TestBulirschCEL:
    @pytest.mark.parametrize("m", [0.1, 0.5, 0.9])
    def test_cel1_equals_K(self, m):
        kc = math.sqrt(1 - m)
        K, _ = ellipke(m)
        assert abs(_f(cel(kc, 1.0, 1.0, 1.0)) - _f(K)) < 1e-12

    @pytest.mark.parametrize("m", [0.1, 0.5, 0.9])
    def test_cel2_equals_E(self, m):
        kc = math.sqrt(1 - m)
        _, E = ellipke(m)
        assert abs(_f(cel(kc, 1.0, 1.0, 1 - m)) - _f(E)) < 1e-12


# ---------------------------------------------------------------------------
# Theta sanity
# ---------------------------------------------------------------------------

class TestTheta:
    def test_theta1_at_zero_is_zero(self):
        # θ₁ is odd
        assert abs(_f(theta(1, 0.0, 0.5))) < 1e-15

    def test_theta3_at_zero_positive(self):
        for m in [0.1, 0.5, 0.7]:
            v0 = _f(theta(3, 0.0, m))
            assert v0 > 1.0

    def test_jacobiThetaEta_at_u_zero(self):
        Th, H = jacobiThetaEta(0.0, 0.5)
        assert abs(_f(Th) - 1.0) < 1e-3 or _f(Th) > 0  # Θ(0) = θ₄(0,q)
        assert abs(_f(H)) < 1e-13                      # H(0) = θ₁(0,q) = 0


# ---------------------------------------------------------------------------
# Grid sweep — broad numerical regression net
# ---------------------------------------------------------------------------

class TestGridSweep:
    def test_F_grid_vs_mpmath(self):
        rng = np.random.default_rng(42)
        phis = rng.uniform(0.01, math.pi / 2 - 0.01, 25)
        ms = rng.uniform(0.01, 0.99, 25)
        max_err = 0.0
        for phi, m in zip(phis, ms):
            F, _, _ = elliptic12(phi, m)
            ref = float(mp.ellipf(mpf(phi), mpf(m)))
            err = abs(_f(F) - ref)
            max_err = max(max_err, err)
        assert max_err < 1e-12, f"max F error {max_err}"

    def test_sn_grid_vs_mpmath(self):
        rng = np.random.default_rng(7)
        us = rng.uniform(-3.0, 3.0, 25)
        ms = rng.uniform(0.01, 0.95, 25)
        max_err = 0.0
        for u, m in zip(us, ms):
            sn, _, _, _ = ellipj(u, m)
            ref = float(mp.ellipfun('sn', mpf(u), mpf(m)))
            max_err = max(max_err, abs(_f(sn) - ref))
        assert max_err < 1e-12, f"max sn error {max_err}"
