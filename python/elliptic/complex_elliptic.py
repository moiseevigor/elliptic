"""Elliptic functions for complex arguments.

elliptic12i(u, m)  — F(u|m), E(u|m), Z(u|m) for complex u
ellipji(u, m)      — sn(u|m), cn(u|m), dn(u|m) for complex u

Algorithms follow A&S §17.4 (Milne-Thomson) — decompose the complex
argument into two real-valued calls, then combine analytically.
"""
from __future__ import annotations
import numpy as np

from ._xputils import get_xp
from .elliptic12 import _elliptic12_xp
from .ellipj import _ellipj_xp


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
    xp = get_xp(u, m)
    u = xp.asarray(u, dtype=xp.complex128)
    m = xp.asarray(m, dtype=xp.float64)
    u_f, m_f = xp.broadcast_arrays(u, m)

    if xp is np and np.any((m_f < 0.0) | (m_f > 1.0)):
        raise ValueError("m must be in [0, 1]")

    phi = xp.real(u_f)
    psi = xp.imag(u_f)

    # Avoid cot(phi) singularity at phi = 0
    eps = np.finfo(np.float64).eps
    phi_s = xp.where(xp.abs(phi) < eps, xp.full_like(phi, eps), phi)

    # Roots of   X² - b*X - c = 0   (A&S 17.4.11)
    cot2   = (xp.cos(phi_s) / xp.sin(phi_s))**2
    sinh2  = xp.sinh(psi)**2
    csc2   = 1.0 / xp.sin(phi_s)**2
    b      = -(cot2 + m_f * sinh2 * csc2 - 1.0 + m_f)
    c      = -(1.0 - m_f) * cot2

    disc = xp.sqrt(xp.maximum(b**2 / 4.0 - c, xp.zeros_like(c)))

    # c <= 0, so the roots straddle zero and -b/2 + disc is the non-negative
    # one.  Near phi = pi/2 that form cancels catastrophically (both terms
    # ~ |b|/2 while the root ~ 0), so use the equal -c/(b/2 + disc) when b > 0.
    den = xp.where(b > 0, b / 2.0 + disc, xp.ones_like(b))
    X = xp.where(b > 0, -c / den, -b / 2.0 + disc)
    ratio = xp.where(
        b > 0,
        (1.0 - m_f) / den,
        (-b / 2.0 + disc) / cot2,
    )

    lam = xp.arctan(1.0 / xp.sqrt(xp.maximum(X, xp.zeros_like(X)) + 1e-300))
    # tan(mu)² = (tan(phi)²·cot(lam)² - 1)/m, taken from *ratio* rather than
    # from lam: at phi = pi/2 the root X underflows, lam rounds to exactly
    # pi/2 and cot(lam) loses every digit of it, collapsing Im to zero.
    m_calc = xp.where(m_f == 0.0, xp.ones_like(m_f), m_f)
    mu = xp.arctan(
        xp.sqrt(xp.maximum((ratio - 1.0) / m_calc, xp.zeros_like(ratio)))
    )

    # Account for periodicity
    lam = (
        (-1.0) ** xp.floor(phi / np.pi * 2.0) * lam
        + np.pi * xp.ceil(phi / np.pi - 0.5 + 1e-14)
    )
    mu = xp.sign(psi) * xp.real(mu)

    F1, E1, _ = _elliptic12_xp(xp, lam, m_f)
    F2, E2, _ = _elliptic12_xp(xp, mu, 1.0 - m_f)

    Fi = F1 + 1j * F2

    # E addition formula (A&S 17.4.16)
    sl = xp.sin(lam);  cl = xp.cos(lam)
    sm = xp.sin(mu);   cm = xp.cos(mu)
    d2l = 1.0 - m_f * sl**2
    d2m = 1.0 - (1.0 - m_f) * sm**2
    den = cm**2 + m_f * sl**2 * sm**2
    b1  = m_f * sl * cl * sm**2 * xp.sqrt(d2l)
    b2  = sm * cm * d2l * xp.sqrt(d2m)
    Ei  = (b1 + 1j * b2) / den + E1 + 1j * (-E2 + F2)

    # Z = E - (E_complete / K) * F
    K_m, E_m, _ = _elliptic12_xp(xp, xp.full_like(m_f, np.pi * 0.5), m_f)
    Zi   = Ei - (E_m / K_m) * Fi

    # Small-m Maclaurin series (through m^2).  The A&S 17.4.11 decomposition
    # loses ~sqrt(eps/m) digits as m -> 0 (0.2 absolute at m = 1e-16); the
    # series is exact there and covers m = 0 itself:
    #   F = u + m(u/4 - sin2u/8) + m^2(9u/64 - 3sin2u/32 + 3sin4u/256) + O(m^3)
    #   E = u - m(u/4 - sin2u/8) - m^2(3u/64 -  sin2u/32 +  sin4u/256) + O(m^3)
    # Valid while |m sin^2 u| is small: switch on m*max(1, e^(2|psi|)) < 1e-4,
    # where the crossover error is ~2e-12 (measured against 40-digit mpmath).
    m_eff = m_f * xp.maximum(xp.ones_like(m_f), xp.exp(2.0 * xp.abs(psi)))
    small = m_eff < 1e-4
    s2 = xp.sin(2.0 * u_f)
    s4 = xp.sin(4.0 * u_f)
    F_ser = (u_f + m_f * (u_f / 4.0 - s2 / 8.0)
             + m_f**2 * (9.0 * u_f / 64.0 - 3.0 * s2 / 32.0 + 3.0 * s4 / 256.0))
    E_ser = (u_f - m_f * (u_f / 4.0 - s2 / 8.0)
             - m_f**2 * (3.0 * u_f / 64.0 - s2 / 32.0 + s4 / 256.0))
    Z_ser = E_ser - (E_m / K_m) * F_ser
    Fi = xp.where(small, F_ser, Fi)
    Ei = xp.where(small, E_ser, Ei)
    Zi = xp.where(small, Z_ser, Zi)
    return Fi, Ei, Zi


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
    xp = get_xp(u, m)
    u_f = xp.asarray(u, dtype=xp.complex128)
    m_f = xp.asarray(m, dtype=xp.float64)
    u_f, m_f = xp.broadcast_arrays(u_f, m_f)

    if xp is np and np.any((m_f < 0.0) | (m_f > 1.0)):
        raise ValueError("m must be in [0, 1]")

    phi = xp.real(u_f)
    psi = xp.imag(u_f)

    s,  c,  d,  _ = _ellipj_xp(xp, phi, m_f)
    s1, c1, d1, _ = _ellipj_xp(xp, psi, 1.0 - m_f)

    delta = c1**2 + m_f * s**2 * s1**2

    sni = (s * d1 + 1j * c * d * s1 * c1) / delta
    cni = (c * c1 - 1j * s * d * s1 * d1) / delta
    dni = (d * c1 * d1 - 1j * m_f * s * c * s1) / delta

    return sni, cni, dni
