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
from .carlson import _rf_xp


def _reject_complex_inputs(*values):
    for value in values:
        dtype = getattr(value, "dtype", None)
        if dtype is not None and "complex" in str(dtype).lower():
            raise ValueError("Weierstrass functions currently support real inputs only")
        if dtype is None and isinstance(value, complex):
            raise ValueError("Weierstrass functions currently support real inputs only")

def _broadcast4(z, e1, e2, e3):
    _reject_complex_inputs(z, e1, e2, e3)
    xp = get_xp(z, e1, e2, e3)
    # Python-scalar roots must be materialised on the SAME DEVICE as the
    # array inputs: xp.asarray(1.5) is a CPU tensor in torch, and
    # broadcast_arrays refuses to mix it with CUDA tensors (found on L4).
    args = [z, e1, e2, e3]
    ref = next((a for a in args if hasattr(a, 'device') or hasattr(a, 'shape')), None)
    def dev(a):
        if ref is not None and not hasattr(a, 'shape'):
            return xp.full_like(xp.asarray(ref, dtype=xp.float64), float(a))
        return xp.asarray(a, dtype=xp.float64)
    z, e1, e2, e3 = (dev(a) for a in args)
    # Root ordering e1 >= e2 >= e3 with e1 > e3 (equal neighbours are the
    # legitimate m = 0 / m = 1 degenerate lattices).  Unsorted roots used to
    # fall through as m > 1 and return NaN silently; MATLAB errors there.
    from ._xputils import is_numpy
    ok = (e1 >= e2) & (e2 >= e3) & (e1 > e3)
    if is_numpy(xp) and not bool(xp.all(ok | xp.isnan(e1 + e2 + e3))):
        raise ValueError("Weierstrass roots must satisfy e1 >= e2 >= e3 with e1 > e3")
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
    mp = (e1 - e2) / (e1 - e3)          # 1-m without cancellation
    scale = xp.sqrt(e1 - e3)
    K = _rf_xp(xp, xp.zeros_like(m), mp, xp.ones_like(m))
    omega1 = K / scale
    period = xp.round(z / (2.0 * omega1))
    z_reduced = z - 2.0 * period * omega1
    w = z_reduced * scale
    sn, _, _, _ = _ellipj_xp(xp, w, m)
    sn2 = sn * sn
    # Pole only at the exact lattice point (DLMF 23.9.2: the Laurent
    # expansion makes every nearby representable z a huge FINITE value --
    # P(1e-16) ~ 1e32, not Inf; a tolerance here destroys that data).
    pole = z_reduced == 0.0
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
    xp, z, e1, e2, e3 = _broadcast4(z, e1, e2, e3)
    return _weierZ_xp(xp, z, e1, e2, e3)


def _lattice_theta_xp(xp, z, e1, e2, e3):
    """omega1, eta1 and the theta1 series values needed by zeta/sigma.

    Closed theta forms (DLMF 23.6.8/9/13): no quadrature.  Returns
    (omega1, eta1, th1, th1p) with the common factor 2 of A&S 16.27
    dropped from every series -- only ratios are used downstream, except
    th1/th1p0 where both drop the same factor.
    """
    m_param = (e2 - e3) / (e1 - e3)
    mp_param = (e1 - e2) / (e1 - e3)    # 1-m without cancellation: near
    # m -> 1 lattices, 1.0 - m loses the digits omega1 and eta1 depend on
    zero = xp.zeros_like(m_param)
    one = xp.ones_like(m_param)
    K = _rf_xp(xp, zero, mp_param, one)
    Kp = _rf_xp(xp, zero, m_param, one)
    omega1 = K / xp.sqrt(e1 - e3)
    q = xp.exp(-math.pi * Kp / K)
    v = math.pi * z / (2.0 * omega1)

    th1 = xp.zeros_like(v)
    th1p = xp.zeros_like(v)
    th1p0 = xp.zeros_like(v)
    th1ppp0 = xp.zeros_like(v)
    # sin/cos of (2n+1)v by angle-addition from sin v, cos v: the products
    # k*v round by eps*|k v| (1e-12 at v ~ 1e8) -- see theta._trig_start
    sk, ck = xp.sin(v), xp.cos(v)
    s2, c2 = 2.0 * sk * ck, 1.0 - 2.0 * sk * sk
    for n in range(31):
        qq = (-1.0) ** n * q ** ((n + 0.5) ** 2)
        k = 2 * n + 1
        th1 = th1 + qq * sk
        th1p = th1p + qq * k * ck
        th1p0 = th1p0 + qq * k
        th1ppp0 = th1ppp0 - qq * k ** 3
        sk, ck = sk * c2 + ck * s2, ck * c2 - sk * s2
    eta1 = -math.pi ** 2 / (12.0 * omega1) * th1ppp0 / th1p0
    return omega1, eta1, th1, th1p, th1p0


def _weierZ_xp(xp, z, e1, e2, e3):
    # zeta(z) = eta1*z/omega1 + pi/(2*omega1) * theta1'(v)/theta1(v)  [DLMF 23.6.13]
    # Quasi-periodicity zeta(z + 2k*omega1) = zeta(z) + 2k*eta1 is carried
    # exactly by the formula; no period reduction needed.
    omega1, eta1, th1, th1p, _ = _lattice_theta_xp(xp, z, e1, e2, e3)
    period = xp.round(z / (2.0 * omega1))
    z_reduced = z - 2.0 * period * omega1
    # Pole only at the exact lattice point (DLMF 23.9.2: the Laurent
    # expansion makes every nearby representable z a huge FINITE value --
    # P(1e-16) ~ 1e32, not Inf; a tolerance here destroys that data).
    pole = z_reduced == 0.0
    ratio = th1p / xp.where(pole, xp.ones_like(th1), th1)
    Z = eta1 * z / omega1 + math.pi / (2.0 * omega1) * ratio
    return xp.where(pole, xp.full_like(Z, math.inf), Z)


def _weierZ_numpy(z, e1, e2, e3):
    return _weierZ_xp(np, z, e1, e2, e3)


# ---------------------------------------------------------------------------
# Weierstrass Sigma
# ---------------------------------------------------------------------------

def weierstrassSigma(z, e1, e2, e3):
    """Weierstrass sigma function (entire, odd, sigma'(z)/sigma(z) = zeta(z))."""
    xp, z, e1, e2, e3 = _broadcast4(z, e1, e2, e3)
    return _weierS_xp(xp, z, e1, e2, e3)


def _weierS_xp(xp, z, e1, e2, e3):
    # sigma(z) = 2*omega1/pi * exp(eta1*z^2/(2*omega1)) * theta1(v)/theta1'(0)
    # [DLMF 23.6.9].  Entire function: every lattice zero and sign change
    # comes out of theta1 itself.  The previous form integrated log(sigma)
    # through the zeta pole at 2*omega1 and was wrong beyond it.
    omega1, eta1, th1, _, th1p0 = _lattice_theta_xp(xp, z, e1, e2, e3)
    return (
        2.0
        * omega1
        / math.pi
        * xp.exp(eta1 * z * z / (2.0 * omega1))
        * th1
        / th1p0
    )


def _weierS_numpy(z, e1, e2, e3):
    return _weierS_xp(np, z, e1, e2, e3)


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
    mp    = (e1 - e2) / (e1 - e3)
    root_scale = xp.sqrt(e1 - e3)
    K = _rf_xp(xp, xp.zeros_like(m), mp, xp.ones_like(m))
    omega1 = K / root_scale
    period = xp.round(z / (2.0 * omega1))
    z_reduced = z - 2.0 * period * omega1
    w = z_reduced * root_scale
    sn, cn, dn, _ = _ellipj_xp(xp, w, m)
    scale = -2.0 * (e1 - e3) ** 1.5
    # Pole only at the exact lattice point (DLMF 23.9.2: the Laurent
    # expansion makes every nearby representable z a huge FINITE value --
    # P(1e-16) ~ 1e32, not Inf; a tolerance here destroys that data).
    pole = z_reduced == 0.0
    dP    = scale * cn * dn / xp.where(pole, xp.ones_like(sn), sn * sn * sn)
    return xp.where(pole, xp.full_like(dP, math.inf), dP)
