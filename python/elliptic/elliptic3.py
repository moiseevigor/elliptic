"""Incomplete elliptic integral of the third kind Π(u, m, n).

    Pi(u, m, n) = integral_0^u  1 / ((1 - n sin^2 t) sqrt(1 - m sin^2 t))  dt

Algorithm: 10-point Gauss-Legendre quadrature (Zhang & Jin 1996).
Fills the gap left by scipy.special (issue #4452 open since 2015).
"""
from __future__ import annotations

import numpy as np
from array_api_compat import array_namespace

# 10-point Gauss-Legendre nodes and weights on [0, 1]
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


def elliptic3(u, m, n):
    """Incomplete elliptic integral of the third kind.

    Parameters
    ----------
    u : array_like
        Phase in radians, 0 <= u <= pi/2.
    m : array_like
        Parameter, 0 <= m <= 1.
    n : array_like
        Characteristic, 0 <= n < 1.

    Returns
    -------
    Pi : array
        Value of Pi(u|m, n).
    """
    u = np.asarray(u, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    xp = array_namespace(u, m, n)
    u, m, n = xp.broadcast_arrays(u, m, n)
    orig_shape = u.shape

    Pi = _elliptic3_numpy(
        np.asarray(u).ravel(),
        np.asarray(m).ravel(),
        np.asarray(n).ravel(),
    )
    return xp.asarray(Pi.reshape(orig_shape))


def _integrand(t, m, n):
    sn2 = np.sin(t) ** 2
    return 1.0 / ((1.0 - n * sn2) * np.sqrt(np.clip(1.0 - m * sn2, 0.0, None)))


def _elliptic3_numpy(u: np.ndarray, m: np.ndarray, n: np.ndarray):
    N = u.size
    Pi = np.zeros(N)

    # Special: u == pi/2 and (m == 1 or n == 1) → inf
    inf_mask = ((u == np.pi / 2) & (m == 1.0)) | ((u == np.pi / 2) & (n == 1.0))

    t_nodes = _GL_T  # shape (10,)
    w_nodes = _GL_W

    P = np.zeros(N)
    half_u = u / 2.0
    for i in range(10):
        c0 = half_u * t_nodes[i]
        P += w_nodes[i] * (
            _integrand(half_u + c0, m, n) +
            _integrand(half_u - c0, m, n)
        )
    Pi = half_u * P
    Pi[inf_mask] = np.inf
    return Pi
