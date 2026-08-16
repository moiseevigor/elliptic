"""Incomplete elliptic integral of the third kind Π(u, m, n).

    Pi(u, m, n) = integral_0^u  1 / ((1 - n sin^2 t) sqrt(1 - m sin^2 t))  dt

Algorithm: Carlson symmetric forms (DLMF 19.25.1). Pure array-namespace
operations run natively on NumPy, PyTorch CUDA, and JAX.
"""
from __future__ import annotations

import math
import numpy as np

from ._xputils import get_xp
from .carlson import _rf_xp, _rj_xp


def elliptic3(u, m, n):
    """Incomplete elliptic integral of the third kind.

    Parameters
    ----------
    u : array_like   Phase in radians.
    m : array_like   Parameter, 0 <= m <= 1.
    n : array_like   Characteristic, n <= 1. For n > 1 the integral is a
                     Cauchy principal value (circular case, DLMF 19.7.3)
                     which this real-valued implementation does not resolve;
                     NumPy calls raise ValueError when the integration path
                     crosses that pole.

    Returns
    -------
    Pi : array
    """
    xp = get_xp(u, m, n)
    u = xp.asarray(u, dtype=xp.float64)
    m = xp.asarray(m, dtype=xp.float64)
    n = xp.asarray(n, dtype=xp.float64)
    u, m, n = xp.broadcast_arrays(u, m, n)

    # Eager NumPy calls can provide a precise domain error.  Traced backends
    # cannot branch on array values; their invalid elements naturally become
    # non-finite through the Carlson expression instead.
    if xp is np and np.any(n > 1.0):
        n_np = np.asarray(n)
        u_np = np.asarray(u)
        # Check whether the singularity sin²θ = 1/n lies in [0, u]
        with np.errstate(invalid="ignore", divide="ignore"):
            sing = np.where(n_np > 1.0, np.arcsin(np.sqrt(1.0 / n_np)), np.inf)
        if np.any((n_np > 1.0) & (np.abs(u_np) >= sing)):
            raise ValueError(
                "elliptic3: n > 1 with phase beyond the pole at arcsin(1/sqrt(n)) "
                "is a Cauchy principal-value integral (DLMF 19.7.3); not supported. "
                "Use a transformation (DLMF 19.7.4) "
                "or compute via Carlson R_J with complex arguments."
            )

    # Reduce the phase to [0, pi/2] using oddness and the quasi-period
    # (the integrand is pi-periodic and even about every multiple of pi/2):
    #   Pi(-u)      = -Pi(u)
    #   Pi(u+k*pi)  =  Pi(u) + 2k*Pi(pi/2)
    #   Pi(pi-u)    =  2*Pi(pi/2) - Pi(u)
    sign_u = xp.where(u < 0, -xp.ones_like(u), xp.ones_like(u))
    ua     = xp.abs(u)
    k_per  = xp.floor(ua / math.pi)
    r      = ua - k_per * math.pi                       # in [0, pi)
    refl   = r > math.pi * 0.5
    u_red  = xp.where(refl, math.pi - r, r)             # in [0, pi/2]
    s = xp.sin(u_red)
    c = xp.cos(u_red)
    s2 = s * s
    d2 = 1.0 - m * s2
    p = 1.0 - n * s2
    one = xp.ones_like(s)

    RF = _rf_xp(xp, c * c, d2, one)
    RJ = _rj_xp(xp, c * c, d2, one, p)
    Pi_red = s * RF + n * s * s2 * RJ / 3.0
    Pi_red = xp.where(s == 0.0, xp.zeros_like(Pi_red), Pi_red)

    # Complete Π is needed to restore reflected/full periods.  Replace
    # singular complete parameters while evaluating the eager expression so
    # unselected branches stay finite on autodiff backends.
    complete_singular = (m == 1.0) | (n >= 1.0)
    m_complete = xp.where(complete_singular, xp.zeros_like(m), m)
    n_complete = xp.where(complete_singular, xp.zeros_like(n), n)
    zero = xp.zeros_like(s)
    RF_complete = _rf_xp(xp, zero, 1.0 - m_complete, one)
    RJ_complete = _rj_xp(xp, zero, 1.0 - m_complete, one, 1.0 - n_complete)
    Pi_complete = RF_complete + n_complete * RJ_complete / 3.0

    Pi_abs = (
        2.0 * k_per * Pi_complete
        + xp.where(refl, 2.0 * Pi_complete - Pi_red, Pi_red)
    )
    Pi = sign_u * Pi_abs

    # At m=1 or n=1 the path diverges once it reaches the first π/2 pole.
    crosses_endpoint_pole = complete_singular & (ua >= math.pi * 0.5)
    signed_inf = sign_u * xp.full_like(Pi, math.inf)
    return xp.where(crosses_endpoint_pole, signed_inf, Pi)
