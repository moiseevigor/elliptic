"""Associate incomplete elliptic integrals B(phi|m), D(phi|m), J(phi,n|m).

    B(phi|m) = integral_0^phi  cos^2 t / sqrt(1 - m sin^2 t)  dt
    D(phi|m) = integral_0^phi  sin^2 t / sqrt(1 - m sin^2 t)  dt
    J(phi,n|m) = integral_0^phi  sin^2 t / ((1-n sin^2 t) sqrt(1-m sin^2 t))  dt

Relations to standard forms:
    F(phi|m)   = B + D
    E(phi|m)   = B + (1-m)*D
    Pi(n,phi|m) = B + D + n*J

Algorithm: Carlson symmetric forms (DLMF §19.25):
    F = sin(phi) * RF(cos^2, Delta^2, 1)
    D = (sin^3/3) * RD(cos^2, Delta^2, 1)
    B = F - D
    J = (sin^3/3) * RJ(cos^2, Delta^2, 1, 1 - n*sin^2)
"""
from __future__ import annotations

from ._xputils import get_xp
from .carlson import _rf_xp, _rd_xp, _rj_xp


def ellipticBDJ(phi, m, n=None):
    """Associate incomplete elliptic integrals B(phi|m), D(phi|m), J(phi,n|m).

    Parameters
    ----------
    phi : array_like   Amplitude in radians.
    m   : array_like   Parameter, 0 <= m < 1.
    n   : array_like, optional  Characteristic for J.

    Returns
    -------
    B, D : arrays
    J    : array or None
    """
    compute_J = n is not None
    args = (phi, m, n) if compute_J else (phi, m)
    xp = get_xp(*args)

    phi = xp.asarray(phi, dtype=xp.float64)
    m   = xp.asarray(m,   dtype=xp.float64)
    if compute_J:
        n = xp.asarray(n, dtype=xp.float64)
        phi, m, n = xp.broadcast_arrays(phi, m, n)
    else:
        phi, m = xp.broadcast_arrays(phi, m)

    s    = xp.sin(phi)
    c    = xp.cos(phi)
    d2   = 1.0 - m * s * s
    s3o3 = s * s * s / 3.0

    RF = _rf_xp(xp, c * c, d2, xp.ones_like(phi))
    RD = _rd_xp(xp, c * c, d2, xp.ones_like(phi))

    F_val = s * RF
    D_val = s3o3 * RD
    B_val = F_val - D_val

    zero  = s == 0.0
    B_val = xp.where(zero, xp.zeros_like(B_val), B_val)
    D_val = xp.where(zero, xp.zeros_like(D_val), D_val)

    if compute_J:
        p     = 1.0 - n * s * s
        RJ    = _rj_xp(xp, c * c, d2, xp.ones_like(phi), p)
        J_val = s3o3 * RJ
        J_val = xp.where(zero, xp.zeros_like(J_val), J_val)
    else:
        J_val = None

    return B_val, D_val, J_val
