"""Jacobi-argument associate integrals Eu, Du, Ju.

    E_u(u|m)   = integral_0^u dn^2(v|m) dv  =  E(am(u|m) | m)
    D_u(u|m)   = integral_0^u sn^2(v|m) dv  =  D(am(u|m) | m)
    J_u(u,n|m) = integral_0^u sn^2(v)/(1-n sn^2(v)) dv

Relations:
    B_u + D_u = u
    E_u = u - m * D_u
"""
from __future__ import annotations

from ._xputils import get_xp
from .ellipj import ellipj
from .carlson import _rf_xp, _rd_xp, _rj_xp
from .ellipticBDJ import ellipticBDJ


def jacobiEDJ(u, m, n=None):
    """Jacobi-argument associate integrals.

    Parameters
    ----------
    u : array_like
        Jacobi argument.
    m : array_like
        Parameter, 0 <= m < 1.
    n : array_like, optional
        Characteristic for J_u.

    Returns
    -------
    Eu, Du : arrays
    Ju : array or None
    """
    args = (u, m, n) if n is not None else (u, m)
    xp = get_xp(*args)
    u_arr = xp.asarray(u, dtype=xp.float64)
    m_arr = xp.asarray(m, dtype=xp.float64)
    if n is None:
        u_arr, m_arr = xp.broadcast_arrays(u_arr, m_arr)
    else:
        n = xp.asarray(n, dtype=xp.float64)
        u_arr, m_arr, n = xp.broadcast_arrays(u_arr, m_arr, n)

    # Reduce u by the period 2K BEFORE taking the amplitude.  am(u) at
    # |u| ~ 1e3 carries eps*|am| ~ 5e-14 of rounding, and near phi = pi/2
    # with m -> 1 the map phi -> D(phi|m) is steep (dD/dphi = sin^2/Delta
    # ~ 1/sqrt(1-m)): D_u(1520|1-1e-8) was off by 3e-10.  As functions of
    # u the integrals are perfectly conditioned (dD_u/du = sn^2 <= 1).
    zed = xp.zeros_like(m_arr)
    one = xp.ones_like(m_arr)
    K = _rf_xp(xp, zed, 1.0 - m_arr, one)
    k = xp.floor((u_arr + K) / (2.0 * K))
    u_r = u_arr - 2.0 * k * K
    _, _, _, phi_r = ellipj(u_r, m_arr)               # |phi_r| <= pi/2
    B_r, D_r, J_r = ellipticBDJ(phi_r, m_arr, n)
    D_cpl = _rd_xp(xp, zed, 1.0 - m_arr, one) / 3.0   # D(m)
    Du = D_r + 2.0 * k * D_cpl
    Eu = u_arr - m_arr * Du                            # error eps*|u|: the conditioning floor
    if n is not None:
        # J(n|m) only where a period was removed (n = 1 is a pole there)
        n_safe = xp.where(k == 0.0, zed, n)
        J_cpl = _rj_xp(xp, zed, 1.0 - m_arr, one, 1.0 - n_safe) / 3.0
        Ju = J_r + xp.where(k == 0.0, zed, 2.0 * k * J_cpl)
    else:
        Ju = None
    return Eu, Du, Ju
