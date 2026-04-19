"""Jacobi elliptic functions sn, cn, dn, am — native on any array backend.

Algorithm: Arithmetic-Geometric Mean + descending Landen back-substitution
(Abramowitz & Stegun §16.4).  Fixed 25 AGM iterations, no per-element
convergence tracking → fully data-parallel on CUDA / JAX.
"""
from __future__ import annotations

import numpy as np

from ._xputils import get_xp

_AGM_ITERS = 25


def ellipj(u, m):
    """Jacobi elliptic functions sn(u|m), cn(u|m), dn(u|m), am(u|m).

    Parameters
    ----------
    u : array_like  Argument.
    m : array_like  Parameter, 0 <= m <= 1.

    Returns
    -------
    sn, cn, dn, am : arrays with broadcast shape of *u* and *m*.
    """
    xp = get_xp(u, m)
    u = xp.asarray(u, dtype=xp.float64)
    m = xp.asarray(m, dtype=xp.float64)
    u, m = xp.broadcast_arrays(u, m)
    return _ellipj_xp(xp, u, m)


def _ellipj_xp(xp, u, m):
    # Clamp m away from exact 0/1 so sqrt is defined; edge cases handled below.
    m_safe = xp.clip(m, 1e-15, 1.0 - 1e-15)

    a = xp.ones_like(m_safe)
    b = xp.sqrt(1.0 - m_safe)

    # Forward AGM: store ratio = (a-b)/(a+b) = c_new/a_new for back-sub
    ratios = []
    for _ in range(_AGM_ITERS):
        ab_sum = a + b
        ratios.append((a - b) / ab_sum)
        b = xp.sqrt(a * b)
        a = ab_sum * 0.5

    # Starting amplitude:  phi_N = 2^N * a_N * u
    phin = (2.0 ** _AGM_ITERS) * a * u

    # Descending Landen back-substitution (all elements, fixed 25 steps)
    for i in range(_AGM_ITERS - 1, -1, -1):
        arg  = xp.clip(ratios[i] * xp.sin(phin), -1.0, 1.0)
        phin = 0.5 * (xp.arcsin(arg) + phin)

    sn_g = xp.sin(phin)
    cn_g = xp.cos(phin)
    dn_g = xp.sqrt(xp.clip(1.0 - m_safe * sn_g * sn_g, 0.0, None))

    # Blend exact m=0 and m=1 results
    sn = xp.where(m == 0.0, xp.sin(u),
         xp.where(m == 1.0, xp.tanh(u), sn_g))
    cn = xp.where(m == 0.0, xp.cos(u),
         xp.where(m == 1.0, 1.0 / xp.cosh(u), cn_g))
    dn = xp.where(m == 0.0, xp.ones_like(u),
         xp.where(m == 1.0, 1.0 / xp.cosh(u), dn_g))
    am = xp.arcsin(xp.clip(sn, -1.0, 1.0))
    return sn, cn, dn, am


def _ellipj_numpy(u, m):
    """Legacy alias for internal callers."""
    u = np.asarray(u, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64)
    return _ellipj_xp(np, u, m)
