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
from .elliptic12 import _elliptic12_xp


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
    phi = xp.full_like(m, math.pi * 0.5)
    K, E, _ = _elliptic12_xp(xp, phi, m)

    # D = (K - E) / m,  limit at m=0: π/4
    D = xp.where(m == 0.0,
                 xp.full_like(m, math.pi * 0.25),
                 (K - E) / xp.where(m == 0.0, xp.ones_like(m), m))
    B = K - D
    S = xp.where(m == 0.0,
                 xp.full_like(m, math.pi / 16.0),
                 (D - B) / xp.where(m == 0.0, xp.ones_like(m), m))
    return B, D, S


def _bd_numpy(m: np.ndarray):
    """Legacy numpy entry point used by bulirsch.py."""
    m = np.asarray(m, dtype=np.float64)
    return _bd_xp(np, m)
