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

from ._xputils import get_xp
from .carlson import _rf_xp, _rd_xp, _rf_numpy, _rd_numpy


def elliptic12(u, m):
    """Incomplete elliptic integrals F(u|m), E(u|m) and Jacobi Zeta Z(u|m).

    Parameters
    ----------
    u : array_like   Phase in radians.
    m : array_like   Parameter, 0 <= m <= 1.

    Returns
    -------
    F, E, Z : arrays with broadcast shape of *u* and *m*.
    """
    xp = get_xp(u, m)
    u = xp.asarray(u, dtype=xp.float64)
    m = xp.asarray(m, dtype=xp.float64)
    u, m = xp.broadcast_arrays(u, m)
    return _elliptic12_xp(xp, u, m)


def _elliptic12_xp(xp, u, m):
    """Backend-native F, E, Z via Carlson forms.  u and m are 1-D xp arrays."""
    # Period reduction: F(u+kπ|m) = F(u|m) + 2k·K(m), Z period π
    k   = xp.round(u / math.pi)
    u_r = u - k * math.pi          # reduced to (-π/2, π/2]

    # Complete integrals K(m), E(m) via Carlson
    z0  = xp.zeros_like(m)
    o1  = xp.ones_like(m)
    mc  = 1.0 - m
    K_m = _rf_xp(xp, z0, mc, o1)
    Em  = K_m - m / 3.0 * _rd_xp(xp, z0, mc, o1)

    s   = xp.sin(u_r)
    c   = xp.cos(u_r)
    d2  = 1.0 - m * s * s

    RF  = _rf_xp(xp, c * c, d2, xp.ones_like(u_r))
    RD  = _rd_xp(xp, c * c, d2, xp.ones_like(u_r))

    F_r = s * RF
    E_r = F_r - m * s * s * s / 3.0 * RD
    Z_r = E_r - (Em / K_m) * F_r

    # s == 0 → zero (handles u = multiple of π)
    F_r = xp.where(s == 0.0, xp.zeros_like(F_r), F_r)
    E_r = xp.where(s == 0.0, xp.zeros_like(E_r), E_r)
    Z_r = xp.where(s == 0.0, xp.zeros_like(Z_r), Z_r)

    F = F_r + 2.0 * k * K_m
    E = E_r + 2.0 * k * Em
    Z = Z_r

    # m == 0: F = E = u, Z = 0
    F = xp.where(m == 0.0, u, F)
    E = xp.where(m == 0.0, u, E)
    Z = xp.where(m == 0.0, xp.zeros_like(Z), Z)

    # m == 1: F = log(tan(π/4 + u_r/2)), E via sin, Z = sin(u_r)
    F_m1 = xp.log(xp.tan(math.pi / 4 + u_r * 0.5))
    um1  = xp.abs(u_r)
    Nf   = xp.floor((um1 + math.pi * 0.5) / math.pi)
    sgn  = xp.where(u >= 0.0, xp.ones_like(u), -xp.ones_like(u))
    E_m1 = ((-1.0) ** Nf * xp.sin(um1) + 2.0 * Nf) * sgn
    Z_m1 = xp.sin(u_r)           # (-1)^Nf * sin(u), Nf=0 for |u_r|<π/2

    near_pole_m1 = xp.abs(u_r) >= math.pi * 0.5 - 1e-14
    F_m1 = xp.where(near_pole_m1, xp.full_like(F_m1, math.inf) * sgn, F_m1)
    F = xp.where(m == 1.0, F_m1, F)
    E = xp.where(m == 1.0, E_m1, E)
    Z = xp.where(m == 1.0, Z_m1, Z)

    return F, E, Z


def _elliptic12_numpy(u: np.ndarray, m: np.ndarray, eps: float):
    """Legacy numpy entry point (flat 1-D in, flat 1-D out)."""
    u = np.asarray(u, dtype=np.float64).ravel()
    m = np.asarray(m, dtype=np.float64).ravel()
    return _elliptic12_xp(np, u, m)
