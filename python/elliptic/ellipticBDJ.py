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

import math

from ._xputils import get_xp, is_numpy, sub_kpi
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

    # The Carlson forms below are only valid for |phi| <= pi/2.  All three
    # integrands are pi-periodic, so reduce phi and restore the periods:
    #   B(phi+k*pi|m)   = B(phi|m)   + 2k*B(m)
    #   D(phi+k*pi|m)   = D(phi|m)   + 2k*D(m)
    #   J(phi+k*pi,n|m) = J(phi,n|m) + 2k*J(n|m)
    k   = xp.ceil(phi / math.pi - 0.5)
    # Cody-Waite split of pi: (u - k*PI_HI) - k*PI_LO keeps the reduction error at eps*|u_r| instead of eps*|u|
    phi = sub_kpi(phi, k)                                     # now in (-pi/2, pi/2], error eps*|phi|

    s    = xp.sin(phi)
    c    = xp.cos(phi)
    d2   = 1.0 - m * s * s
    s3o3 = s * s * s / 3.0

    one  = xp.ones_like(phi)
    zed  = xp.zeros_like(phi)

    RF = _rf_xp(xp, c * c, d2, one)
    RD = _rd_xp(xp, c * c, d2, one)

    F_val = s * RF
    D_val = s3o3 * RD
    B_val = F_val - D_val

    zero  = s == 0.0
    B_val = xp.where(zero, xp.zeros_like(B_val), B_val)
    D_val = xp.where(zero, xp.zeros_like(D_val), D_val)

    # Complete associate integrals — same Carlson forms at phi = pi/2 (s=1, c=0)
    D_cpl = _rd_xp(xp, zed, 1.0 - m, one) / 3.0          # D(m)
    B_cpl = _rf_xp(xp, zed, 1.0 - m, one) - D_cpl        # B(m) = K(m) - D(m)

    B_val = B_val + 2.0 * k * B_cpl
    D_val = D_val + 2.0 * k * D_cpl

    if compute_J:
        p     = 1.0 - n * s * s
        # n > 1 with the phase beyond the pole at arcsin(1/sqrt(n)) is a
        # Cauchy principal-value integral (DLMF 19.7.3): R_J needs p > 0.
        # The private _rj_xp silently returned garbage there (1.147 for
        # J(1, 1.5|0.5); the principal value is 0.859).  The complete J(n|m)
        # is only needed when a period was removed (k != 0); at n = 1 it is
        # a pole, and 0 * inf turned J(phi, 1|m) into NaN for |phi| <= pi/2.
        bad_inc = p <= 0.0
        bad_cpl = (k != 0.0) & (n >= 1.0)
        if is_numpy(xp) and bool(xp.any(bad_inc | bad_cpl)):
            raise ValueError(
                "ellipticBDJ: n > 1 with phase beyond the pole at arcsin(1/sqrt(n)) "
                "(or n >= 1 with |phi| > pi/2) is a Cauchy principal-value integral "
                "(DLMF 19.7.3); not supported.")
        p_safe = xp.where(bad_inc, one, p)
        RJ    = _rj_xp(xp, c * c, d2, one, p_safe)
        J_val = s3o3 * RJ
        J_val = xp.where(zero, xp.zeros_like(J_val), J_val)
        n_safe = xp.where(bad_cpl | (n == 1.0), xp.zeros_like(n), n)
        J_cpl = _rj_xp(xp, zed, 1.0 - m, one, 1.0 - n_safe) / 3.0   # J(n|m)
        J_val = J_val + xp.where(k == 0.0, xp.zeros_like(J_val), 2.0 * k * J_cpl)
        J_val = xp.where(bad_inc | bad_cpl, xp.full_like(J_val, math.nan), J_val)
    else:
        J_val = None

    return B_val, D_val, J_val
