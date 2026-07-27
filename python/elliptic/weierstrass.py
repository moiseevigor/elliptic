"""Weierstrass elliptic functions P, Zeta, Sigma.

All three parameterised by roots (e1, e2, e3) with e1 > e2 > e3.

References
----------
Abramowitz & Stegun §18.9–18.10; NIST DLMF §23.
"""
from __future__ import annotations

import math
import numpy as np
from ._xputils import get_xp

def _broadcast4(z, e1, e2, e3):
    xp = get_xp(z, e1, e2, e3)
    z  = xp.asarray(z,  dtype=xp.float64)
    e1 = xp.asarray(e1, dtype=xp.float64)
    e2 = xp.asarray(e2, dtype=xp.float64)
    e3 = xp.asarray(e3, dtype=xp.float64)
    z, e1, e2, e3 = xp.broadcast_arrays(z, e1, e2, e3)
    return xp, z, e1, e2, e3


# ---------------------------------------------------------------------------
# Weierstrass P
# ---------------------------------------------------------------------------

def weierstrassP(z, e1, e2, e3):
    """Weierstrass P-function.  e1 > e2 > e3."""
    xp, z, e1, e2, e3 = _broadcast4(z, e1, e2, e3)
    return _weierP_xp(xp, z, e1, e2, e3)


def _weierP_xp(xp, z, e1, e2, e3):
    from .ellipj import _ellipj_xp
    m  = (e2 - e3) / (e1 - e3)
    w  = z * xp.sqrt(e1 - e3)
    sn, _, _, _ = _ellipj_xp(xp, w, m)
    sn2 = sn * sn
    pole = xp.abs(sn) < 1e-10
    P = e3 + (e1 - e3) / xp.where(pole, xp.ones_like(sn2), sn2)
    return xp.where(pole, xp.full_like(P, math.inf), P)


def _weierP_numpy(z, e1, e2, e3):
    return _weierP_xp(np, np.asarray(z, dtype=np.float64),
                          np.asarray(e1, dtype=np.float64),
                          np.asarray(e2, dtype=np.float64),
                          np.asarray(e3, dtype=np.float64))


# ---------------------------------------------------------------------------
# Weierstrass Zeta
# ---------------------------------------------------------------------------

def weierstrassZeta(z, e1, e2, e3):
    """Weierstrass zeta function (NOT Riemann zeta)."""
    _, z, e1, e2, e3 = _broadcast4(z, e1, e2, e3)
    orig_shape = z.shape
    Z = _weierZ_numpy(np.asarray(z).ravel(), np.asarray(e1).ravel(),
                      np.asarray(e2).ravel(), np.asarray(e3).ravel())
    return np.asarray(Z.reshape(orig_shape))


def _lattice_theta_numpy(z, e1, e2, e3):
    """omega1, eta1 and the theta1 series values needed by zeta/sigma.

    Closed theta forms (DLMF 23.6.8/9/13): no quadrature.  Returns
    (omega1, eta1, th1, th1p) with the common factor 2 of A&S 16.27
    dropped from every series -- only ratios are used downstream, except
    th1/th1p0 where both drop the same factor.
    """
    from .elliptic12 import _elliptic12_xp
    m_param = (e2 - e3) / (e1 - e3)
    phi_half = np.full_like(m_param, math.pi / 2)
    K,  _, _ = _elliptic12_xp(np, phi_half, m_param)
    Kp, _, _ = _elliptic12_xp(np, phi_half, 1.0 - m_param)
    omega1 = K / np.sqrt(e1 - e3)
    q = np.exp(-math.pi * Kp / K)
    v = math.pi * z / (2.0 * omega1)

    qmax = float(np.max(q)) if q.size else 0.0
    nT = min(30, max(2, math.ceil(math.sqrt(abs(math.log(np.finfo(float).eps)
                                                / math.log(qmax)))))) if qmax > 0 else 1

    th1 = np.zeros_like(v); th1p = np.zeros_like(v)
    th1p0 = np.zeros_like(v); th1ppp0 = np.zeros_like(v)
    for n in range(nT + 1):
        qq = (-1.0) ** n * q ** ((n + 0.5) ** 2)
        k = 2 * n + 1
        th1     += qq * np.sin(k * v)
        th1p    += qq * k * np.cos(k * v)
        th1p0   += qq * k
        th1ppp0 -= qq * k ** 3
    eta1 = -math.pi ** 2 / (12.0 * omega1) * th1ppp0 / th1p0
    return omega1, eta1, th1, th1p, th1p0


def _weierZ_numpy(z, e1, e2, e3):
    # zeta(z) = eta1*z/omega1 + pi/(2*omega1) * theta1'(v)/theta1(v)  [DLMF 23.6.13]
    # Quasi-periodicity zeta(z + 2k*omega1) = zeta(z) + 2k*eta1 is carried
    # exactly by the formula; no period reduction needed.
    omega1, eta1, th1, th1p, _ = _lattice_theta_numpy(z, e1, e2, e3)
    with np.errstate(divide='ignore', invalid='ignore'):
        Z = eta1 * z / omega1 + math.pi / (2.0 * omega1) * th1p / th1
    Z[th1 == 0.0] = np.inf          # lattice points z = 2k*omega1
    return Z


# ---------------------------------------------------------------------------
# Weierstrass Sigma
# ---------------------------------------------------------------------------

def weierstrassSigma(z, e1, e2, e3):
    """Weierstrass sigma function (entire, odd, sigma'(z)/sigma(z) = zeta(z))."""
    _, z, e1, e2, e3 = _broadcast4(z, e1, e2, e3)
    orig_shape = z.shape
    S = _weierS_numpy(np.asarray(z).ravel(), np.asarray(e1).ravel(),
                      np.asarray(e2).ravel(), np.asarray(e3).ravel())
    return np.asarray(S.reshape(orig_shape))


def _weierS_numpy(z, e1, e2, e3):
    # sigma(z) = 2*omega1/pi * exp(eta1*z^2/(2*omega1)) * theta1(v)/theta1'(0)
    # [DLMF 23.6.9].  Entire function: every lattice zero and sign change
    # comes out of theta1 itself.  The previous form integrated log(sigma)
    # through the zeta pole at 2*omega1 and was wrong beyond it.
    omega1, eta1, th1, _, th1p0 = _lattice_theta_numpy(z, e1, e2, e3)
    return 2.0 * omega1 / math.pi * np.exp(eta1 * z * z / (2.0 * omega1)) * th1 / th1p0


# -----------------------------------------------------------------------
# Weierstrass invariants and P-derivative
# -----------------------------------------------------------------------

def weierstrassInvariants(e1, e2, e3):
    """Lattice invariants g₂, g₃, Δ from the three roots."""
    xp = get_xp(e1, e2, e3)
    e1 = xp.asarray(e1, dtype=xp.float64)
    e2 = xp.asarray(e2, dtype=xp.float64)
    e3 = xp.asarray(e3, dtype=xp.float64)
    e1, e2, e3 = xp.broadcast_arrays(e1, e2, e3)
    g2    = -4.0 * (e1*e2 + e1*e3 + e2*e3)
    g3    =  4.0 * e1 * e2 * e3
    Delta = g2**3 - 27.0 * g3**2
    return g2, g3, Delta


def weierstrassPPrime(z, e1, e2, e3):
    """℘'(z) = -2(e1-e3)^{3/2} · cn·dn / sn³."""
    xp, z, e1, e2, e3 = _broadcast4(z, e1, e2, e3)
    from .ellipj import _ellipj_xp
    m     = (e2 - e3) / (e1 - e3)
    w     = z * xp.sqrt(e1 - e3)
    sn, cn, dn, _ = _ellipj_xp(xp, w, m)
    scale = -2.0 * (e1 - e3) ** 1.5
    pole  = xp.abs(sn) < 1e-10
    dP    = scale * cn * dn / xp.where(pole, xp.ones_like(sn), sn * sn * sn)
    return xp.where(pole, xp.full_like(dP, math.inf), dP)
