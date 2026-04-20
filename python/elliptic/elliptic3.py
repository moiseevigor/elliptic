"""Incomplete elliptic integral of the third kind Π(u, m, n).

    Pi(u, m, n) = integral_0^u  1 / ((1 - n sin^2 t) sqrt(1 - m sin^2 t))  dt

Algorithm: 10-point Gauss-Legendre quadrature. Pure array-namespace ops —
runs natively on NumPy, PyTorch CUDA, and JAX without any np.asarray conversion.
"""
from __future__ import annotations

import math

from ._xputils import get_xp

# 10-point Gauss-Legendre nodes and weights on [0, 1] as plain Python floats
# (so they broadcast correctly against any backend tensor)
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


def elliptic3(u, m, n):
    """Incomplete elliptic integral of the third kind.

    Parameters
    ----------
    u : array_like   Phase in radians, 0 <= u <= pi/2.
    m : array_like   Parameter, 0 <= m <= 1.
    n : array_like   Characteristic, n <= 1. For n > 1 the integral is a
                     Cauchy principal value (circular case, DLMF 19.7.3)
                     which 10-point Gauss–Legendre cannot resolve; this
                     function raises ValueError in that regime.

    Returns
    -------
    Pi : array
    """
    xp = get_xp(u, m, n)
    u = xp.asarray(u, dtype=xp.float64)
    m = xp.asarray(m, dtype=xp.float64)
    n = xp.asarray(n, dtype=xp.float64)
    u, m, n = xp.broadcast_arrays(u, m, n)

    import numpy as _np
    n_np = _np.asarray(n)
    u_np = _np.asarray(u)
    if _np.any(n_np > 1.0):
        # Check whether the singularity sin²θ = 1/n lies in [0, u]
        with _np.errstate(invalid='ignore', divide='ignore'):
            sing = _np.where(n_np > 1.0, _np.arcsin(_np.sqrt(1.0 / n_np)), _np.inf)
        if _np.any((n_np > 1.0) & (u_np >= sing)):
            raise ValueError(
                "elliptic3: n > 1 with phase beyond the pole at arcsin(1/sqrt(n)) "
                "is a Cauchy principal-value integral (DLMF 19.7.3); not supported "
                "by 10-point Gauss–Legendre. Use a transformation (DLMF 19.7.4) "
                "or compute via Carlson R_J with complex arguments."
            )

    half_u = u * 0.5
    P = xp.zeros_like(u)
    for ti, wi in zip(_GL_T, _GL_W):
        c0  = half_u * ti
        tp  = half_u + c0
        tm  = half_u - c0
        s2p = xp.sin(tp) ** 2
        s2m = xp.sin(tm) ** 2
        P = P + wi * (
            1.0 / ((1.0 - n * s2p) * xp.sqrt(xp.clip(1.0 - m * s2p, 0.0, None))) +
            1.0 / ((1.0 - n * s2m) * xp.sqrt(xp.clip(1.0 - m * s2m, 0.0, None)))
        )
    Pi = half_u * P

    # u == pi/2 and (m == 1 or n == 1) → inf
    inf_mask = ((u == math.pi * 0.5) & (m == 1.0)) | ((u == math.pi * 0.5) & (n == 1.0))
    return xp.where(inf_mask, xp.full_like(Pi, math.inf), Pi)
