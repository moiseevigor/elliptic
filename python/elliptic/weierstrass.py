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

# 10-point GL nodes/weights as plain Python floats
_GL_T = [
    0.9931285991850949, 0.9639719272779138,
    0.9122344282513259, 0.8391169718222188,
    0.7463319064601508, 0.6360536807265150,
    0.5108670019508271, 0.3737060887154195,
    0.2277858511416451, 0.07652652113349734,
]
_GL_W = [
    0.01761400713915212, 0.04060142980038694,
    0.06267204833410907, 0.08327674157670475,
    0.10193011981724040, 0.11819453196151840,
    0.13168863844917660, 0.14209610931838200,
    0.14917298647260370, 0.15275338713072580,
]


def _broadcast4(z, e1, e2, e3):
    xp = get_xp(z, e1, e2, e3)
    z  = xp.asarray(z,  dtype=xp.float64)
    e1 = xp.asarray(e1, dtype=xp.float64)
    e2 = xp.asarray(e2, dtype=xp.float64)
    e3 = xp.asarray(e3, dtype=xp.float64)
    z, e1, e2, e3 = xp.broadcast_arrays(z, e1, e2, e3)
    return xp, z, e1, e2, e3


def _gl_integral_numpy(f, a, b, e1, e2, e3):
    """10-pt GL quadrature of f(t, e1, e2, e3) on [a, b]."""
    mid = 0.5 * (a + b)
    half = 0.5 * (b - a)
    I = np.zeros_like(a)
    for ti, wi in zip(_GL_T, _GL_W):
        t_p = mid + half * ti
        t_m = mid - half * ti
        I += wi * (f(t_p, e1, e2, e3) + f(t_m, e1, e2, e3))
    return half * I


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


def _weierZ_numpy(z, e1, e2, e3):
    from .elliptic12 import _elliptic12_xp
    m_param = (e2 - e3) / (e1 - e3)
    phi_half = np.full_like(m_param, math.pi / 2)
    K, _, _ = _elliptic12_xp(np, phi_half, m_param)
    omega1 = K / np.sqrt(e1 - e3)

    g2 = -4.0 * (e1 * e2 + e1 * e3 + e2 * e3)
    g3 =  4.0 * e1 * e2 * e3

    eps0 = 0.1 * omega1
    anal = -g2 * eps0 ** 3 / 60.0 - g3 * eps0 ** 5 / 140.0

    def feta(t, e1, e2, e3):
        P = _weierP_numpy(t, e1, e2, e3)
        v = -P + 1.0 / t ** 2
        g2_ = -4.0 * (e1 * e2 + e1 * e3 + e2 * e3)
        g3_ =  4.0 * e1 * e2 * e3
        small = np.abs(t) < 1e-4
        v[small] = -g2_[small] * t[small] ** 2 / 20.0 - g3_[small] * t[small] ** 4 / 28.0
        return v

    gl_eta = _gl_integral_numpy(feta, eps0, omega1, e1, e2, e3)
    eta1 = 1.0 / omega1 + anal + gl_eta

    k = np.floor((z + omega1 * (1.0 - 1e-14)) / (2.0 * omega1))
    z_red = z - k * 2.0 * omega1
    neg = z_red < 0.0
    az_red = np.abs(z_red)

    def fzeta(t, e1, e2, e3):
        v = -_weierP_numpy(t, e1, e2, e3)
        v[~np.isfinite(v)] = 0.0
        return v

    int_val = _gl_integral_numpy(fzeta, omega1, az_red, e1, e2, e3)
    Z_red = eta1 + int_val
    Z_red[neg] = -Z_red[neg]
    Z = Z_red + k * 2.0 * eta1
    Z[az_red < 1e-5] = np.inf
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
    from .elliptic12 import _elliptic12_xp
    m_param = (e2 - e3) / (e1 - e3)
    phi_half = np.full_like(m_param, math.pi / 2)
    K, _, _ = _elliptic12_xp(np, phi_half, m_param)
    omega1 = K / np.sqrt(e1 - e3)

    g2 = -4.0 * (e1 * e2 + e1 * e3 + e2 * e3)
    g3 =  4.0 * e1 * e2 * e3

    t_split = np.minimum(0.25 * omega1, np.abs(z))

    # Analytic part of integral on [0, t_split]: integral of zeta(t)-1/t
    # Laurent: zeta(t)-1/t = -g2*t^3/60 - g3*t^5/140 + ...
    # integral = -g2*t_split^4/240 - g3*t_split^6/840
    anal = -g2 * t_split ** 4 / 240.0 - g3 * t_split ** 6 / 840.0

    def flog(t, e1, e2, e3):
        """zeta(t) - 1/t (integrand for log sigma)."""
        Z = _weierZ_numpy(t, e1, e2, e3)
        v = Z - 1.0 / t
        v[~np.isfinite(v)] = 0.0
        return v

    gl_part = _gl_integral_numpy(flog, t_split, np.abs(z), e1, e2, e3)
    log_sigma = np.log(np.abs(z) + (z == 0.0)) + anal + gl_part
    S = np.sign(z) * np.exp(log_sigma)
    S[z == 0.0] = 0.0
    return S


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
