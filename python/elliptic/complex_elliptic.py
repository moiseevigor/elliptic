"""Elliptic functions for complex arguments.

elliptic12i(u, m)  — F(u|m), E(u|m), Z(u|m) for complex u
ellipji(u, m)      — sn(u|m), cn(u|m), dn(u|m) for complex u

Algorithms follow A&S §17.4 (Milne-Thomson) — decompose the complex
argument into two real-valued calls, then combine analytically.
"""
from __future__ import annotations
import numpy as np
from array_api_compat import array_namespace

from .elliptic12 import _elliptic12_numpy
from .ellipj import _ellipj_numpy


def elliptic12i(u, m):
    """Incomplete elliptic integrals F, E, Z for complex phase u.

    Parameters
    ----------
    u : array_like (complex)
        Complex phase in radians.
    m : array_like (real)
        Parameter 0 <= m <= 1 (scalar or same shape as u).

    Returns
    -------
    Fi, Ei, Zi : complex arrays
    """
    u = np.asarray(u, dtype=np.complex128)
    m = np.asarray(m, dtype=np.float64)
    xp = array_namespace(np.real(u), m)

    u_bc, m_bc = np.broadcast_arrays(u, m)
    orig_shape = u_bc.shape
    u_f = u_bc.ravel()
    m_f = m_bc.ravel().astype(np.float64)

    if np.any(m_f < 0) or np.any(m_f > 1):
        raise ValueError("m must be in [0, 1]")

    phi = np.real(u_f)
    psi = np.imag(u_f)

    # Avoid cot(phi) singularity at phi = 0
    phi_s = np.where(np.abs(phi) < np.finfo(float).eps, np.finfo(float).eps, phi)

    # Roots of   X² - b*X - c = 0   (A&S 17.4.11)
    cot2   = (np.cos(phi_s) / np.sin(phi_s))**2
    sinh2  = np.sinh(psi)**2
    csc2   = 1.0 / np.sin(phi_s)**2
    b      = -(cot2 + m_f * sinh2 * csc2 - 1.0 + m_f)
    c      = -(1.0 - m_f) * cot2

    disc  = np.sqrt(np.maximum(b**2 / 4.0 - c, 0.0))
    X1    = -b / 2.0 + disc
    X2    = -b / 2.0 - disc

    X     = np.where(X1 >= 0, X1, X2)

    lam  = np.arctan(1.0 / np.sqrt(np.maximum(X, 0.0) + 1e-300))
    mu   = np.arctan(np.sqrt(np.maximum(
               1.0 / m_f * (np.tan(phi_s)**2 * (np.cos(lam) / np.sin(lam))**2 - 1.0),
               0.0)))

    # Account for periodicity
    lam = (-1.0)**np.floor(phi / np.pi * 2) * lam + np.pi * np.ceil(phi / np.pi - 0.5 + 1e-14)
    mu  = np.sign(psi) * np.real(mu)

    F1, E1, _ = _elliptic12_numpy(lam, m_f, np.finfo(np.float64).eps)
    F2, E2, _ = _elliptic12_numpy(mu,  1.0 - m_f, np.finfo(np.float64).eps)

    Fi = F1 + 1j * F2

    # E addition formula (A&S 17.4.16)
    sl = np.sin(lam);  cl = np.cos(lam)
    sm = np.sin(mu);   cm = np.cos(mu)
    d2l = 1.0 - m_f * sl**2
    d2m = 1.0 - (1.0 - m_f) * sm**2
    den = cm**2 + m_f * sl**2 * sm**2
    b1  = m_f * sl * cl * sm**2 * np.sqrt(d2l)
    b2  = sm * cm * d2l * np.sqrt(d2m)
    Ei  = (b1 + 1j * b2) / den + E1 + 1j * (-E2 + F2)

    # Z = E - (E_complete / K) * F
    from scipy.special import ellipk, ellipe as _ellipe
    K_m  = ellipk(m_f)
    E_m  = _ellipe(m_f)
    Zi   = Ei - (E_m / K_m) * Fi

    Fi = Fi.reshape(orig_shape)
    Ei = Ei.reshape(orig_shape)
    Zi = Zi.reshape(orig_shape)
    return xp.asarray(Fi), xp.asarray(Ei), xp.asarray(Zi)


def ellipji(u, m):
    """Jacobi elliptic functions sn, cn, dn for complex argument u.

    Uses the addition formulae (A&S 17.4.14-16):
        sn(x+iy|m) via sn(x|m) and sn(y|1-m)

    Parameters
    ----------
    u : array_like (complex)
        Complex argument.
    m : array_like (real)
        Parameter 0 <= m <= 1.

    Returns
    -------
    sn, cn, dn : complex arrays
    """
    u = np.asarray(u, dtype=np.complex128)
    m = np.asarray(m, dtype=np.float64)
    xp = array_namespace(np.real(u), m)

    u_bc, m_bc = np.broadcast_arrays(u, m)
    orig_shape = u_bc.shape
    u_f = u_bc.ravel()
    m_f = m_bc.ravel().astype(np.float64)

    if np.any(m_f < 0) or np.any(m_f > 1):
        raise ValueError("m must be in [0, 1]")

    phi = np.real(u_f)
    psi = np.imag(u_f)

    s,  c,  d,  _ = _ellipj_numpy(phi, m_f)
    s1, c1, d1, _ = _ellipj_numpy(psi, 1.0 - m_f)

    delta = c1**2 + m_f * s**2 * s1**2

    sni = (s * d1 + 1j * c * d * s1 * c1) / delta
    cni = (c * c1 - 1j * s * d * s1 * d1) / delta
    dni = (d * c1 * d1 - 1j * m_f * s * c * s1) / delta

    sni = sni.reshape(orig_shape)
    cni = cni.reshape(orig_shape)
    dni = dni.reshape(orig_shape)
    return xp.asarray(sni), xp.asarray(cni), xp.asarray(dni)
