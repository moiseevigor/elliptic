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
from .ellipticBD import _bd_xp
from .elliptic12 import _elliptic12_xp
from .carlson import _rj_xp


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
    m    = 1.0 - kc * kc
    phi  = xp.full_like(m, math.pi * 0.5)
    K, _, _ = _elliptic12_xp(xp, phi, m)
    B, D, _ = _bd_xp(xp, m)

    # p ≈ 1 branch: C = a*B + b*D
    C_p1 = a * B + b * D

    # p ≠ 1 branch: C = a*K + (b - a*p)*(Pi - K)/(1-p)
    n_val  = 1.0 - p
    mc     = 1.0 - m
    n_safe = xp.where(xp.abs(n_val) < 1e-14, xp.ones_like(n_val), n_val)
    RJ     = _rj_xp(xp, xp.zeros_like(m), mc, xp.ones_like(m), p)
    J_n    = RJ / 3.0
    Pi_n   = K + n_val * J_n
    C_pn   = a * K + (b - a * p) * (Pi_n - K) / n_safe

    C = xp.where(xp.abs(p - 1.0) < 1e-12, C_p1, C_pn)
    C = xp.where(kc < 0.0, xp.full_like(C, math.nan), C)
    C = xp.where(p <= 0.0, xp.full_like(C, math.inf), C)
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
