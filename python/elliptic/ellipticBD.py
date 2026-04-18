"""Complete associate elliptic integrals B(m), D(m), S(m).

    B(m) = B(pi/2 | m)
    D(m) = D(pi/2 | m) = (K(m) - E(m)) / m
    S(m) = (D - B) / m

Relations:
    K(m) = B(m) + D(m)
    E(m) = B(m) + (1-m)*D(m)
"""
from __future__ import annotations

import numpy as np
from array_api_compat import array_namespace
from .elliptic12 import elliptic12


def ellipticBD(m):
    """Complete associate elliptic integrals B(m), D(m), S(m).

    Parameters
    ----------
    m : array_like
        Parameter, 0 <= m < 1.

    Returns
    -------
    B, D, S : arrays
    """
    m = np.asarray(m, dtype=np.float64)
    xp = array_namespace(m)
    orig_shape = m.shape
    m_np = np.asarray(m).ravel()

    B_np, D_np, S_np = _bd_numpy(m_np)
    return (xp.asarray(B_np.reshape(orig_shape)),
            xp.asarray(D_np.reshape(orig_shape)),
            xp.asarray(S_np.reshape(orig_shape)))


def _bd_numpy(m: np.ndarray):
    """Use K and E at phi=pi/2 computed via elliptic12."""
    import math
    phi = np.full_like(m, math.pi / 2)
    # elliptic12 returns (F, E, Z) — at phi=pi/2 F=K, E=E(m)
    F_np, E_np, _ = _elliptic12_numpy_direct(phi, m)
    K = F_np
    E = E_np

    D = np.zeros_like(m)
    nz = m != 0.0
    D[nz]  = (K[nz] - E[nz]) / m[nz]
    D[~nz] = np.pi / 4.0

    B = K - D

    S = np.zeros_like(m)
    S[nz]  = (D[nz] - B[nz]) / m[nz]
    S[~nz] = np.pi / 16.0

    return B, D, S


def _elliptic12_numpy_direct(phi, m):
    """Call the numpy core of elliptic12 directly (avoids circular import)."""
    from .elliptic12 import _elliptic12_numpy
    eps = np.finfo(np.float64).eps
    return _elliptic12_numpy(phi, m, eps)
