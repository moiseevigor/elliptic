"""Jacobi elliptic functions sn, cn, dn, am — native on any array backend.

Accuracy limit for large arguments: the phase is reduced modulo 2K in
double precision, so the residual carries an absolute uncertainty ~|u|*eps.
Full precision holds for |u| up to ~1e12; by |u| ~ 1e16 the phase is lost
entirely (a bound shared by every double implementation, scipy included).

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
    # Keep every representable interior parameter unchanged.  The old global
    # clip to [1e-15, 1-1e-15] silently changed valid inputs near both
    # endpoints (most visibly m=nextafter(1, 0)).  Exact endpoints use a
    # harmless interior placeholder here and are replaced below.
    interior = (m > 0.0) & (m < 1.0)
    m_safe = xp.where(interior, m, xp.full_like(m, 0.5))

    a = xp.ones_like(m_safe)
    b = xp.sqrt(1.0 - m_safe)

    # Forward AGM: store ratio = (a-b)/(a+b) = c_new/a_new for back-sub
    ratios = []
    for _ in range(_AGM_ITERS):
        ab_sum = a + b
        ratios.append((a - b) / ab_sum)
        b = xp.sqrt(a * b)
        a = ab_sum * 0.5

    # Reduce u before multiplying by 2^N.  Without this, the large argument
    # enters the Landen recursion directly and loses phase bits; near m=1 the
    # error can become O(1) after only a few dozen periods.
    K = np.pi / (2.0 * a)
    period = xp.floor((u + K) / (2.0 * K))
    u_reduced = u - 2.0 * period * K

    # Starting amplitude on [-K, K]: phi_N = 2^N * a_N * u_reduced
    phin = (2.0 ** _AGM_ITERS) * a * u_reduced

    # Descending Landen back-substitution (all elements, fixed 25 steps)
    for i in range(_AGM_ITERS - 1, -1, -1):
        arg  = xp.clip(ratios[i] * xp.sin(phin), -1.0, 1.0)
        phin = 0.5 * (xp.arcsin(arg) + phin)

    period_mod2 = period - 2.0 * xp.floor(period * 0.5)
    quasi_sign = 1.0 - 2.0 * period_mod2
    sn_g = quasi_sign * xp.sin(phin)
    cn_g = quasi_sign * xp.cos(phin)
    # The cn form avoids subtracting two nearly equal numbers when m and
    # |sn| are both close to one.
    dn_g = xp.sqrt(xp.clip((1.0 - m_safe) + m_safe * cn_g * cn_g, 0.0, None))
    am_g = phin + period * np.pi

    # Stable sech avoids overflow in cosh for large non-m=1 elements.  Array
    # backends evaluate both sides of where, so a nominally unselected cosh
    # still emitted warnings/overflowed during ordinary calls.
    exp_neg = xp.exp(-xp.abs(u))
    sech_u = 2.0 * exp_neg / (1.0 + exp_neg * exp_neg)

    # Blend exact m=0 and m=1 results
    sn = xp.where(m == 0.0, xp.sin(u),
         xp.where(m == 1.0, xp.tanh(u), sn_g))
    cn = xp.where(m == 0.0, xp.cos(u),
         xp.where(m == 1.0, sech_u, cn_g))
    dn = xp.where(m == 0.0, xp.ones_like(u),
         xp.where(m == 1.0, sech_u, dn_g))
    # am is the *continuous* amplitude from the Landen recursion, not
    # arcsin(sn): the latter folds it into [-pi/2, pi/2] and so loses the
    # period count (DLMF 22.16.1: am(u + 2K) = am(u) + pi).
    am = xp.where(m == 0.0, u,
         xp.where(m == 1.0, xp.arcsin(xp.clip(xp.tanh(u), -1.0, 1.0)), am_g))
    return sn, cn, dn, am


def _ellipj_numpy(u, m):
    """Legacy alias for internal callers."""
    u = np.asarray(u, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64)
    return _ellipj_xp(np, u, m)
