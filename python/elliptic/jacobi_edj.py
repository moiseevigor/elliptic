"""Jacobi-argument associate integrals Eu, Du, Ju.

    E_u(u|m)   = integral_0^u dn^2(v|m) dv  =  E(am(u|m) | m)
    D_u(u|m)   = integral_0^u sn^2(v|m) dv  =  D(am(u|m) | m)
    J_u(u,n|m) = integral_0^u sn^2(v)/(1-n sn^2(v)) dv

Relations:
    B_u + D_u = u
    E_u = u - m * D_u
"""
from __future__ import annotations

from array_api_compat import array_namespace
from .ellipj import ellipj
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
    compute_J = n is not None
    xp = array_namespace(u, m) if not compute_J else array_namespace(u, m, n)

    _, _, _, phi = ellipj(u, m)
    B, D, J = ellipticBDJ(phi, m, n)

    m_arr = xp.asarray(m, dtype=xp.float64)
    u_arr = xp.asarray(u, dtype=xp.float64)
    m_arr, u_arr = xp.broadcast_arrays(m_arr, u_arr)

    Du = D
    Eu = u_arr - m_arr * D
    Ju = J
    return Eu, Du, Ju
