"""Inverse incomplete elliptic integral of the second kind."""
from __future__ import annotations
import numpy as np
from array_api_compat import array_namespace


def inverselliptic2(E_val, m, tol=None):
    """Inverse of the incomplete elliptic integral of the second kind.

    Solves  E(phi | m) = E_val  for phi using the empirical initialisation
    of Boyd (2012) followed by Newton iterations.

    Parameters
    ----------
    E_val : array_like
        Target value(s) of E(phi | m).
    m : array_like
        Parameter, 0 <= m <= 1.  Scalar or same shape as *E_val*.
    tol : float, optional
        Not used; kept for API compatibility with the MATLAB signature.

    Returns
    -------
    phi : array
        Amplitude phi (radians) such that E(phi | m) ≈ E_val.
    """
    from scipy.special import ellipe, ellipeinc

    E_val = np.asarray(E_val, dtype=np.float64)
    m     = np.asarray(m,     dtype=np.float64)
    xp    = array_namespace(E_val, m)

    E_flat = np.asarray(E_val).ravel()
    m_flat = np.broadcast_to(np.asarray(m), np.broadcast_shapes(np.asarray(E_val).shape, np.asarray(m).shape)).ravel()
    m_flat = m_flat.copy()

    if np.any(m_flat < 0) or np.any(m_flat > 1):
        raise ValueError("m must be in [0, 1]")

    m_flat = np.where(m_flat < np.finfo(float).eps, 0.0, m_flat)

    # Complete second-kind integral for initialisation
    E1 = ellipe(m_flat)

    z      = E_flat
    mu     = 1.0 - m_flat
    zeta   = 1.0 - z / E1
    r      = np.sqrt(zeta**2 + mu**2)
    theta  = np.arctan2(mu, z + np.finfo(float).eps)

    # Boyd "empirical" initialisation
    phi = np.pi / 2.0 + np.sqrt(r) * (theta - np.pi / 2.0)

    # Four Newton iterations:  phi_{n+1} = phi_n - (E(phi_n) - z) / sqrt(1 - m sin^2(phi_n))
    for _ in range(4):
        E_cur = ellipeinc(phi, m_flat)
        denom = np.sqrt(np.maximum(1.0 - m_flat * np.sin(phi)**2, 0.0))
        phi   = phi - (E_cur - z) / np.where(denom < 1e-15, 1e-15, denom)

    orig_shape = np.broadcast_shapes(np.asarray(E_val).shape, np.asarray(m).shape)
    phi = phi.reshape(orig_shape) if orig_shape else phi.squeeze()
    return xp.asarray(phi)
