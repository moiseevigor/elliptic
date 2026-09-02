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

import math
import numpy as np

from ._xputils import get_xp


def cel(kc, p, a, b):
    """Bulirsch generalised complete elliptic integral."""
    xp = get_xp(kc, p, a, b)
    kc = xp.asarray(kc, dtype=xp.float64)
    p  = xp.asarray(p,  dtype=xp.float64)
    a  = xp.asarray(a,  dtype=xp.float64)
    b  = xp.asarray(b,  dtype=xp.float64)
    kc, p, a, b = xp.broadcast_arrays(kc, p, a, b)
    return _cel_xp(xp, kc, p, a, b)


def _cel_xp(xp, kc, p, a, b):
    """Bulirsch's algorithm (Numer. Math. 13 (1969) 305), backend-native.

    Works directly with kc: the previous route through m = 1 - kc**2 lost kc
    entirely below ~1e-8 (cel1(1e-9) returned 2e6; ln(4/kc) = 22.1).  Any
    real kc (kc > 1 is m < 0); p < 0 is the Cauchy principal value; p = 0
    gives inf.  Quadratic Landen ascent run for a fixed number of steps with
    converged elements frozen (no data-dependent branching).
    """
    CA = 1e-9
    k = xp.abs(kc)
    zero_kc = k == 0.0
    # kc = 0: divergent unless b = 0, where the kc -> 0 limit is finite and
    # the ascent evaluates it from the smallest normal number.
    k = xp.where(zero_kc, xp.full_like(k, 2.2250738585072014e-308), k)
    e = k
    em = xp.ones_like(k)
    pos = p > 0.0
    p_safe = xp.where(pos, p, xp.ones_like(p))
    # p > 0
    sp = xp.sqrt(p_safe)
    # p <= 0: transform to the p > 0 case (principal value)
    g0 = 1.0 - p
    g0_safe = xp.where(pos, xp.ones_like(g0), g0)
    f0 = k * k - p
    q0 = (1.0 - k * k) * (b - a * p)
    pn = xp.sqrt(xp.where(pos, xp.ones_like(f0), f0 / g0_safe))
    an = (a - b) / g0_safe
    bn = -q0 / (g0_safe * g0_safe * pn) + an * pn
    p = xp.where(pos, sp, pn)
    b = xp.where(pos, b / sp, bn)
    a = xp.where(pos, a, an)
    active = k == k                      # all-True boolean of the right backend/shape
    for _ in range(40):
        f = a
        a = xp.where(active, a + b / p, a)
        g = e / p
        b = xp.where(active, 2.0 * (b + f * g), b)
        p = xp.where(active, p + g, p)
        g = em
        em = xp.where(active, em + k, em)
        conv = xp.abs(g - k) <= g * CA
        step = active & ~conv
        kk = 2.0 * xp.sqrt(e)
        e = xp.where(step, kk * em, e)
        k = xp.where(step, kk, k)
        active = step
    C = math.pi / 2.0 * (b + a * em) / (em * (em + p))
    C = xp.where(zero_kc & (b != 0.0), xp.sign(b / xp.where(p == 0, xp.ones_like(p), p)) * math.inf, C)
    # kc = NaN never became active above (NaN == NaN is False) and returned
    # the untouched pi/2; propagate it like every other input
    return xp.where(k != k, xp.full_like(C, math.nan), C)


def cel1(kc):
    """K(m) via Bulirsch: cel(kc, 1, 1, 1) where m = 1 - kc^2."""
    return cel(kc, 1.0, 1.0, 1.0)


def cel2(kc, a, b):
    """Bulirsch cel2(kc, a, b) = cel(kc, 1, a, b)."""
    return cel(kc, 1.0, a, b)


def cel3(kc, p):
    """Pi(n|m) via Bulirsch: cel(kc, p, 1, 1) where n = 1-p, m = 1-kc^2."""
    return cel(kc, p, 1.0, 1.0)
