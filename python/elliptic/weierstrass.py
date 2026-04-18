"""Weierstrass elliptic functions P, Zeta, Sigma.

All three parameterised by roots (e1, e2, e3) with e1 > e2 > e3.

References
----------
Abramowitz & Stegun §18.9–18.10; NIST DLMF §23.
"""
from __future__ import annotations

import math
import numpy as np
from array_api_compat import array_namespace

# 10-point GL nodes/weights on [-1, 1]
_GL_T = np.array([
    0.9931285991850949, 0.9639719272779138,
    0.9122344282513259, 0.8391169718222188,
    0.7463319064601508, 0.6360536807265150,
    0.5108670019508271, 0.3737060887154195,
    0.2277858511416451, 0.07652652113349734,
])
_GL_W = np.array([
    0.01761400713915212, 0.04060142980038694,
    0.06267204833410907, 0.08327674157670475,
    0.10193011981724040, 0.11819453196151840,
    0.13168863844917660, 0.14209610931838200,
    0.14917298647260370, 0.15275338713072580,
])


def _broadcast4(xp, z, e1, e2, e3):
    z  = np.asarray(z,  dtype=np.float64)
    e1 = np.asarray(e1, dtype=np.float64)
    e2 = np.asarray(e2, dtype=np.float64)
    e3 = np.asarray(e3, dtype=np.float64)
    xp2 = array_namespace(z, e1, e2, e3)
    z, e1, e2, e3 = xp2.broadcast_arrays(z, e1, e2, e3)
    return z, e1, e2, e3


def _gl_integral_numpy(f, a, b, e1, e2, e3):
    """10-pt GL quadrature of f(t, e1, e2, e3) on [a, b] (vectorised over N)."""
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
    """Weierstrass P-function.

    Parameters
    ----------
    z, e1, e2, e3 : array_like
        Argument and roots (e1 > e2 > e3).

    Returns
    -------
    P : array
    """
    xp = None  # determined inside _broadcast4
    z, e1, e2, e3 = _broadcast4(xp, z, e1, e2, e3)
    xp = array_namespace(z, e1, e2, e3)
    orig_shape = z.shape
    P = _weierP_numpy(np.asarray(z).ravel(), np.asarray(e1).ravel(),
                      np.asarray(e2).ravel(), np.asarray(e3).ravel())
    return xp.asarray(P.reshape(orig_shape))


def _weierP_numpy(z, e1, e2, e3):
    from .ellipj import _ellipj_numpy
    m   = (e2 - e3) / (e1 - e3)
    w   = z * np.sqrt(e1 - e3)
    sn, _, _, _ = _ellipj_numpy(w, m)
    P = e3 + (e1 - e3) / sn ** 2
    P[np.abs(sn) < np.finfo(np.float64).eps ** (1 / 3)] = np.inf
    return P


# ---------------------------------------------------------------------------
# Weierstrass Zeta
# ---------------------------------------------------------------------------

def weierstrassZeta(z, e1, e2, e3):
    """Weierstrass zeta function (NOT Riemann zeta).

    Satisfies: zeta'(z) = -P(z),  zeta(-z) = -zeta(z).
    """
    z, e1, e2, e3 = _broadcast4(None, z, e1, e2, e3)
    xp = array_namespace(z, e1, e2, e3)
    orig_shape = z.shape
    Z = _weierZ_numpy(np.asarray(z).ravel(), np.asarray(e1).ravel(),
                      np.asarray(e2).ravel(), np.asarray(e3).ravel())
    return xp.asarray(Z.reshape(orig_shape))


def _weierZ_numpy(z, e1, e2, e3):
    from .elliptic12 import _elliptic12_numpy
    eps = np.finfo(np.float64).eps
    m_param = (e2 - e3) / (e1 - e3)
    phi_half = np.full_like(m_param, math.pi / 2)
    K, _, _ = _elliptic12_numpy(phi_half, m_param, eps)
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

    k = np.floor((z + omega1 * (1.0 - eps)) / (2.0 * omega1))
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
    Z[az_red < eps ** (1 / 3)] = np.inf
    return Z


# ---------------------------------------------------------------------------
# Weierstrass Sigma
# ---------------------------------------------------------------------------

def weierstrassSigma(z, e1, e2, e3):
    """Weierstrass sigma function (entire, odd, sigma'(z)/sigma(z) = zeta(z))."""
    z, e1, e2, e3 = _broadcast4(None, z, e1, e2, e3)
    xp = array_namespace(z, e1, e2, e3)
    orig_shape = z.shape
    S = _weierS_numpy(np.asarray(z).ravel(), np.asarray(e1).ravel(),
                      np.asarray(e2).ravel(), np.asarray(e3).ravel())
    return xp.asarray(S.reshape(orig_shape))


def _weierS_numpy(z, e1, e2, e3):
    from .elliptic12 import _elliptic12_numpy
    eps = np.finfo(np.float64).eps
    m_param = (e2 - e3) / (e1 - e3)
    phi_half = np.full_like(m_param, math.pi / 2)
    K, _, _ = _elliptic12_numpy(phi_half, m_param, eps)
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
    """Lattice invariants g₂, g₃ and discriminant Δ from the three roots.

    The roots e1, e2, e3 are the values of ℘ at the three half-periods,
    satisfying e1 > e2 > e3 and e1 + e2 + e3 = 0 (A&S 18.1.2).

    g₂ = -4(e1 e2 + e1 e3 + e2 e3)  [A&S 18.1.3]
    g₃ =  4 e1 e2 e3
    Δ  = g₂³ - 27 g₃²

    Parameters
    ----------
    e1, e2, e3 : array_like (real, e1 > e2 > e3)

    Returns
    -------
    g2, g3, Delta : arrays
    """
    e1 = np.asarray(e1, dtype=np.float64)
    e2 = np.asarray(e2, dtype=np.float64)
    e3 = np.asarray(e3, dtype=np.float64)
    xp = array_namespace(e1, e2, e3)
    e1, e2, e3 = (np.asarray(x) for x in (e1, e2, e3))
    e1, e2, e3 = np.broadcast_arrays(e1, e2, e3)
    g2    = -4.0 * (e1*e2 + e1*e3 + e2*e3)
    g3    =  4.0 * e1 * e2 * e3
    Delta = g2**3 - 27.0 * g3**2
    return xp.asarray(g2), xp.asarray(g3), xp.asarray(Delta)


def weierstrassPPrime(z, e1, e2, e3):
    """Derivative ℘'(z; e1, e2, e3) of the Weierstrass P-function.

    ℘'(z) = -2(e1-e3)^{3/2} · cn(w,m) · dn(w,m) / sn³(w,m)

    where  w = z√(e1-e3),  m = (e2-e3)/(e1-e3)  (A&S 18.9.8).

    Verification identity: (℘')² = 4℘³ - g₂℘ - g₃.

    Parameters
    ----------
    z : array_like
    e1, e2, e3 : array_like  (e1 > e2 > e3)

    Returns
    -------
    dP : array
    """
    z  = np.asarray(z,  dtype=np.float64)
    e1 = np.asarray(e1, dtype=np.float64)
    e2 = np.asarray(e2, dtype=np.float64)
    e3 = np.asarray(e3, dtype=np.float64)
    xp = array_namespace(z, e1, e2, e3)
    z, e1, e2, e3 = np.broadcast_arrays(
        np.asarray(z), np.asarray(e1), np.asarray(e2), np.asarray(e3))
    orig_shape = z.shape

    z_f  = z.ravel();  e1_f = e1.ravel(); e2_f = e2.ravel(); e3_f = e3.ravel()
    m    = (e2_f - e3_f) / (e1_f - e3_f)
    w    = z_f * np.sqrt(e1_f - e3_f)
    sn, cn, dn, _ = _ellipj_numpy(w, m)
    scale = -2.0 * (e1_f - e3_f)**1.5
    dP = scale * cn * dn / np.where(np.abs(sn) < np.finfo(float).eps**(1/3), np.nan, sn**3)
    # poles at sn = 0
    dP[np.abs(sn) < np.finfo(float).eps**(1/3)] = np.inf
    return xp.asarray(dP.reshape(orig_shape))


def _ellipj_numpy(w, m):
    """Thin shim — import avoids circular dependency."""
    from .ellipj import _ellipj_numpy as _ej
    return _ej(w, m)
