"""Complete associate elliptic integrals B(m), D(m), S(m).

    B(m) = B(pi/2 | m)
    D(m) = D(pi/2 | m) = (K(m) - E(m)) / m
    S(m) = (D - B) / m

Relations:
    K(m) = B(m) + D(m)
    E(m) = B(m) + (1-m)*D(m)
"""
from __future__ import annotations

import math
import numpy as np

from ._xputils import get_xp
from .carlson import _rf_xp, _rd_xp


def ellipticBD(m):
    """Complete associate elliptic integrals B(m), D(m), S(m).

    Parameters
    ----------
    m : array_like   Parameter, 0 <= m < 1.

    Returns
    -------
    B, D, S : arrays
    """
    xp = get_xp(m)
    m = xp.asarray(m, dtype=xp.float64)
    return _bd_xp(xp, m)


def _bd_xp(xp, m):
    zero = xp.zeros_like(m)
    one = xp.ones_like(m)
    K = _rf_xp(xp, zero, 1.0 - m, one)
    D = _rd_xp(xp, zero, 1.0 - m, one) / 3.0
    B = K - D

    # S=(D-B)/m is catastrophically cancelling near m=0.  Evaluate its
    # convergent binomial/integral series there:
    #   1/sqrt(1-m sin²t) = Σ C_k m^k sin^(2k)t.
    S_series = xp.zeros_like(m)
    for k in range(1, 9):
        ck = math.comb(2 * k, k) / (4.0 ** k)
        ik = math.pi * math.comb(2 * k, k) / (2.0 * 4.0 ** k)
        ik1 = math.pi * math.comb(2 * k + 2, k + 1) / (2.0 * 4.0 ** (k + 1))
        S_series = S_series + ck * (2.0 * ik1 - ik) * m ** (k - 1)
    m_safe = xp.where(m == 0.0, xp.ones_like(m), m)
    S_direct = (D - B) / m_safe
    S = xp.where(xp.abs(m) < 1e-2, S_series, S_direct)
    return B, D, S


def _bd_numpy(m: np.ndarray):
    """Legacy numpy entry point used by bulirsch.py."""
    m = np.asarray(m, dtype=np.float64)
    return _bd_xp(np, m)
