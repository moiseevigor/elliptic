"""Elliptic functions for complex arguments.

elliptic12i(u, m)  — F(u|m), E(u|m), Z(u|m) for complex u
ellipji(u, m)      — sn(u|m), cn(u|m), dn(u|m) for complex u

Algorithms follow A&S §17.4 (Milne-Thomson) — decompose the complex
argument into two real-valued calls, then combine analytically.
"""
from __future__ import annotations
import numpy as np

from ._xputils import get_xp, is_numpy
from .carlson import _rf_xp, _rd_xp
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

    if is_numpy(xp) and np.any((m_f < 0.0) | (m_f > 1.0)):
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

    # Positive root X1 = cot^2(lambda) of X^2 + bX + c = 0 and tan^2(mu),
    # both without cancellation.  Writing X1 = cot^2(phi) + Y, Y solves
    #     Y^2 + B'Y - C' = 0,  B' = cot^2 + (1-m) - m sinh^2 csc^2,
    #                          C' = cot^2 * m sinh^2 csc^2 >= 0,
    # and A&S 17.4.11's tan^2(mu) = (tan^2 phi cot^2 lambda - 1)/m = Y/(m cot^2)
    # collapses to
    #     tan^2(mu) = 2 sinh^2 csc^2 / (B' + sqrt(B'^2 + 4C'))        (B' >= 0)
    #               = (|B'| + sqrt(B'^2 + 4C')) / (2 m cot^2)         (B' <  0)
    # -- m cancels analytically in the first form, so m -> 0 (and m = 0
    # exactly) is handled to full precision.  The old (ratio - 1)/m lost
    # sqrt(eps/m) digits and returned Im F = 0 for psi = 1e-9.
    s2c2 = sinh2 * csc2
    Bp   = cot2 + (1.0 - m_f) - m_f * s2c2
    Cp   = cot2 * m_f * s2c2
    root = xp.sqrt(Bp * Bp + 4.0 * Cp)
    pos  = Bp >= 0.0
    Y    = xp.where(pos, 2.0 * Cp / xp.where(pos, Bp + root, xp.ones_like(Bp)),
                         0.5 * (-Bp + root))
    X    = cot2 + Y
    m_cot = xp.where(pos, xp.ones_like(cot2), m_f * cot2)
    tan2mu = xp.where(pos,
                      2.0 * s2c2 / xp.where(pos, Bp + root, xp.ones_like(Bp)),
                      0.5 * (-Bp + root) / xp.where(pos, xp.ones_like(m_cot), m_cot))

    lam = xp.arctan(1.0 / xp.sqrt(X + 1e-300))
    mu  = xp.arctan(xp.sqrt(tan2mu))

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
    # Complete integrals from the exact Carlson forms, not F(double(pi/2)|m):
    # cos(double(pi/2)) = 6e-17 is not 0 and K was 5.8e-9 relative short at
    # m = 1 - eps/2, which Z inherited (1.8e-11).  E = RF - (m/3) RD (DLMF 19.25.1).
    zed = xp.zeros_like(m_f);  one = xp.ones_like(m_f)
    K_m = _rf_xp(xp, zed, 1.0 - m_f, one)
    E_m = K_m - m_f / 3.0 * _rd_xp(xp, zed, 1.0 - m_f, one)
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
    # u == 0 exactly (incl. -0.0): the cot(phi) nudge above would return eps
    z0 = u_f == 0
    Fi = xp.where(z0, u_f, Fi); Ei = xp.where(z0, u_f, Ei); Zi = xp.where(z0, xp.zeros_like(Zi), Zi)
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

    if is_numpy(xp) and np.any((m_f < 0.0) | (m_f > 1.0)):
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
