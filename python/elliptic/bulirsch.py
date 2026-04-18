"""Bulirsch's generalised complete elliptic integral and special cases.

    cel(kc, p, a, b) = integral_0^inf (a + b t^2) / ((1 + p t^2) sqrt((1+t^2)(1+kc^2 t^2))) dt

Special cases:
    K(m)   = cel(sqrt(1-m), 1, 1, 1)
    E(m)   = cel(sqrt(1-m), 1, 1, 1-m)
    B(m)   = cel(sqrt(1-m), 1, 1, 0)
    D(m)   = cel(sqrt(1-m), 1, 0, 1)
    Pi(n|m) = cel(sqrt(1-m), 1-n, 1, 1)
"""
from __future__ import annotations

import numpy as np
from array_api_compat import array_namespace
from .ellipticBD import _bd_numpy
from .elliptic12 import _elliptic12_numpy
from .carlson import _rj_numpy


def cel(kc, p, a, b):
    """Bulirsch generalised complete elliptic integral.

    Parameters
    ----------
    kc : array_like
        Complementary modulus, kc >= 0.
    p : array_like
        Characteristic parameter, p > 0.
    a, b : array_like
        Numerator coefficients.

    Returns
    -------
    C : array
    """
    kc = np.asarray(kc, dtype=np.float64)
    p  = np.asarray(p,  dtype=np.float64)
    a  = np.asarray(a,  dtype=np.float64)
    b  = np.asarray(b,  dtype=np.float64)
    xp = array_namespace(kc, p, a, b)
    kc, p, a, b = xp.broadcast_arrays(kc, p, a, b)
    orig_shape = kc.shape

    C = _cel_numpy(np.asarray(kc).ravel(), np.asarray(p).ravel(),
                   np.asarray(a).ravel(), np.asarray(b).ravel())
    return xp.asarray(C.reshape(orig_shape))


def _cel_numpy(kc, p, a, b):
    import math
    N = kc.size
    C = np.zeros(N)

    bad  = kc < 0.0
    pole = p <= 0.0
    C[bad]  = np.nan
    C[pole] = np.inf

    ok = ~bad & ~pole
    if not np.any(ok):
        return C

    m  = 1.0 - kc[ok] ** 2
    eps = np.finfo(np.float64).eps

    phi_half = np.full(m.shape, math.pi / 2)
    K, _, _  = _elliptic12_numpy(phi_half, m, eps)  # K(m)
    Bv, Dv, _ = _bd_numpy(m)

    pp = p[ok]; aa = a[ok]; bb = b[ok]

    p1 = np.abs(pp - 1.0) < 1e-12
    pn = ~p1
    Cv = np.zeros(np.sum(ok))

    # p ≈ 1: cel = a*B + b*D
    Cv[p1] = aa[p1] * Bv[p1] + bb[p1] * Dv[p1]

    # p ≠ 1: cel = a*K + (b - a*p)*(Pi - K)/(1-p)
    if np.any(pn):
        m_pn  = m[pn];  p_pn  = pp[pn]
        a_pn  = aa[pn]; b_pn  = bb[pn]; K_pn = K[pn]
        n_val = 1.0 - p_pn   # n for Pi(n|m)
        kc_pn = np.sqrt(1.0 - m_pn)
        # J_complete = (1/3) * RJ(0, kc^2, 1, p)  at phi=pi/2 → s=1, c=0
        RJ_val = _rj_numpy(np.zeros_like(m_pn), kc_pn ** 2, np.ones_like(m_pn), p_pn)
        J_n    = RJ_val / 3.0
        Pi_n   = K_pn + n_val * J_n
        Cv[pn] = a_pn * K_pn + (b_pn - a_pn * p_pn) * (Pi_n - K_pn) / n_val

    C[ok] = Cv
    return C


def cel1(kc):
    """K(m) via Bulirsch: cel(kc, 1, 1, 1) where m = 1 - kc^2."""
    return cel(kc, 1.0, 1.0, 1.0)


def cel2(kc, a, b):
    """Bulirsch cel2(kc, a, b) = cel(kc, 1, a, b)."""
    return cel(kc, 1.0, a, b)


def cel3(kc, p):
    """Pi(n|m) via Bulirsch: cel(kc, p, 1, 1) where n = 1-p, m = 1-kc^2."""
    return cel(kc, p, 1.0, 1.0)
