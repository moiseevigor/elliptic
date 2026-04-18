"""Incomplete elliptic integrals F(phi, m), E(phi, m), Z(phi, m).

    F(phi, m) = integral_0^phi  1 / sqrt(1 - m sin^2 t)  dt
    E(phi, m) = integral_0^phi  sqrt(1 - m sin^2 t)      dt
    Z(phi, m) = E(phi, m) - E(m)/K(m) * F(phi, m)   [Jacobi Zeta]

Algorithm: Carlson symmetric forms (DLMF 19.25.5-6):
    F = sin(phi) * RF(cos^2, 1-m sin^2, 1)
    E = F - m * sin^3(phi)/3 * RD(cos^2, 1-m sin^2, 1)
    Z = E - E(m)/K(m) * F   where K=RF(0,1-m,1), E(m)=K - m/3*RD(0,1-m,1)

This formulation is numerically stable, avoids the Landen descent, and is
consistent with ellipticBDJ.
"""
from __future__ import annotations

import math
import numpy as np
from array_api_compat import array_namespace

from .carlson import _rf_numpy, _rd_numpy


def elliptic12(u, m):
    """Incomplete elliptic integrals of the first and second kind and Jacobi Zeta.

    Parameters
    ----------
    u : array_like
        Phase in radians.
    m : array_like
        Parameter, 0 <= m <= 1.  Scalar or same shape as *u*.

    Returns
    -------
    F, E, Z : arrays
        Values of F(u|m), E(u|m), Z(u|m) with the same shape as the
        broadcast of *u* and *m*.
    """
    u = np.asarray(u, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64)
    xp = array_namespace(u, m)
    u, m = xp.broadcast_arrays(u, m)
    orig_shape = u.shape

    F, E, Z = _elliptic12_numpy(np.asarray(u).ravel(), np.asarray(m).ravel(),
                                np.finfo(np.float64).eps)

    return (xp.asarray(F.reshape(orig_shape)),
            xp.asarray(E.reshape(orig_shape)),
            xp.asarray(Z.reshape(orig_shape)))


def _elliptic12_numpy(u: np.ndarray, m: np.ndarray, eps: float):
    """Pure-NumPy implementation (flat 1-D arrays in, flat 1-D arrays out)."""
    N = u.size
    F = np.zeros(N)
    E = np.zeros(N)
    Z = np.zeros(N)

    # --- special cases ---
    mask0 = m == 0.0
    mask1 = m == 1.0
    maskN = ~mask0 & ~mask1 & (m >= 0.0) & (m <= 1.0)

    # m == 0: F = E = u, Z = 0
    F[mask0] = u[mask0]
    E[mask0] = u[mask0]

    # m == 1
    if np.any(mask1):
        um1 = np.abs(u[mask1])
        Nf = np.floor((um1 + math.pi / 2) / math.pi)
        idx1 = np.where(mask1)[0]
        good = um1 < math.pi / 2
        F[idx1[good]] = np.log(np.tan(math.pi / 4 + u[idx1[good]] / 2))
        F[idx1[~good]] = np.inf * np.sign(u[idx1[~good]])
        E[mask1] = ((-1.0) ** Nf * np.sin(um1) + 2.0 * Nf) * np.sign(u[mask1])
        Z[mask1] = (-1.0) ** Nf * np.sin(u[mask1])

    if not np.any(maskN):
        return F, E, Z

    # --- general case via Carlson forms ---
    u_g = u[maskN]
    m_g = m[maskN]

    s  = np.sin(u_g)
    c  = np.cos(u_g)
    d2 = 1.0 - m_g * s ** 2         # Δ²

    # K(m) and E(m) for Zeta computation
    K_m = _rf_numpy(np.zeros_like(m_g), 1.0 - m_g, np.ones_like(m_g))
    Em  = K_m - m_g / 3.0 * _rd_numpy(np.zeros_like(m_g), 1.0 - m_g, np.ones_like(m_g))

    # Incomplete integrals
    RF = _rf_numpy(c ** 2, d2, np.ones_like(u_g))
    RD = _rd_numpy(c ** 2, d2, np.ones_like(u_g))

    F_g = s * RF
    E_g = F_g - m_g * s ** 3 / 3.0 * RD
    Z_g = E_g - (Em / K_m) * F_g

    # Handle s == 0 (phi == 0 or multiple of pi)
    zero = s == 0.0
    F_g[zero] = 0.0
    E_g[zero] = 0.0
    Z_g[zero] = 0.0

    # sin(phi) already carries the correct sign for odd symmetry.
    F[maskN] = F_g
    E[maskN] = E_g
    Z[maskN] = Z_g

    return F, E, Z
