"""Inverse incomplete elliptic integral of the second kind."""
from __future__ import annotations
import numpy as np


def inverselliptic2(E_val, m, tol=1e-12):
    """Inverse of the incomplete elliptic integral of the second kind.

    Solves  E(phi | m) = E_val  for phi using period reduction, Boyd (2012)
    initialisation, and a Newton while-loop (issue #12: converges to *tol*).

    Period identity:  E(phi + pi | m) = E(phi | m) + 2*E(m)
    Symmetry:         E(pi - phi | m) = 2*E(m) - E(phi | m)

    Parameters
    ----------
    E_val : array_like   Target value(s) of E(phi | m).
    m     : array_like   Parameter, 0 <= m <= 1.
    tol   : float        Newton convergence tolerance (default 1e-12).

    Returns
    -------
    phi : ndarray   Amplitude in radians such that E(phi | m) ≈ E_val.
    """
    from scipy.special import ellipe, ellipeinc

    E_val = np.asarray(E_val, dtype=np.float64)
    m     = np.asarray(m,     dtype=np.float64)

    orig_shape = np.broadcast_shapes(E_val.shape, m.shape)
    E_flat = np.broadcast_to(E_val, orig_shape).ravel().copy()
    m_flat = np.broadcast_to(m,     orig_shape).ravel().copy()

    if np.any(m_flat < 0) or np.any(m_flat > 1):
        raise ValueError("m must be in [0, 1]")
    m_flat = np.where(m_flat < np.finfo(float).eps, 0.0, m_flat)

    # Complete integral E(m); each phi-period of π contributes 2*E1 to E.
    E1 = ellipe(m_flat)
    two_E1 = 2.0 * np.where(E1 > 0, E1, 1.0)

    # Step 1 — strip full periods:  phi = phi_base + k*pi
    k     = np.floor(E_flat / two_E1)
    z_red = E_flat - k * two_E1          # in [0, 2*E1)

    # Step 2 — fold second half-period using E(pi-phi|m) = 2E1 - E(phi|m)
    over   = z_red > E1
    z_red2 = np.where(over, two_E1 - z_red, z_red)   # in [0, E1]

    # Boyd (2012) empirical initialisation for phi in [0, pi/2]
    mu    = 1.0 - m_flat
    zeta  = 1.0 - z_red2 / np.where(E1 > 0, E1, 1.0)
    r     = np.sqrt(zeta**2 + mu**2)
    theta = np.arctan2(mu, z_red2 + 1e-300)
    phi   = np.pi / 2.0 + np.sqrt(r) * (theta - np.pi / 2.0)
    phi   = np.clip(phi, 0.0, np.pi / 2.0)

    # Newton while-loop until converged (issue #12: fixed 4 iters insufficient)
    for _ in range(200):
        E_cur = ellipeinc(phi, m_flat)
        res   = E_cur - z_red2
        if np.max(np.abs(res)) < tol:
            break
        denom = np.sqrt(np.maximum(1.0 - m_flat * np.sin(phi)**2, 0.0))
        phi   = phi - res / np.where(denom < 1e-15, 1e-15, denom)
        phi   = np.clip(phi, 0.0, np.pi / 2.0)

    # Step 3 — undo fold: phi_in_period = pi - phi (if over), else phi
    phi = np.where(over, np.pi - phi, phi)   # in [0, pi)

    # Step 4 — undo period strips: each strip adds pi to phi
    phi = phi + k * np.pi

    return phi.reshape(orig_shape) if orig_shape else phi.squeeze()
