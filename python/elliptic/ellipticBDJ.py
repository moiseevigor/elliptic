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

import numpy as np
from array_api_compat import array_namespace
from .carlson import _rf_numpy, _rd_numpy, _rj_numpy


def ellipticBDJ(phi, m, n=None):
    """Associate incomplete elliptic integrals B, D and (optionally) J.

    Parameters
    ----------
    phi : array_like
        Amplitude in radians.
    m : array_like
        Parameter, 0 <= m < 1.
    n : array_like, optional
        Characteristic for J.  If omitted, J is not computed (returns None).

    Returns
    -------
    B, D : arrays
    J : array or None
    """
    compute_J = n is not None
    phi = np.asarray(phi, dtype=np.float64)
    m   = np.asarray(m,   dtype=np.float64)
    if compute_J:
        n   = np.asarray(n,   dtype=np.float64)
        xp = array_namespace(phi, m, n)
        phi, m, n = xp.broadcast_arrays(phi, m, n)
    else:
        xp = array_namespace(phi, m)
        phi, m = xp.broadcast_arrays(phi, m)

    orig_shape = phi.shape
    phi_np = np.asarray(phi).ravel()
    m_np   = np.asarray(m).ravel()
    n_np   = np.asarray(n).ravel() if compute_J else None

    B_np, D_np, J_np = _bdj_numpy(phi_np, m_np, n_np, compute_J)

    B = xp.asarray(B_np.reshape(orig_shape))
    D = xp.asarray(D_np.reshape(orig_shape))
    J = xp.asarray(J_np.reshape(orig_shape)) if compute_J else None
    return B, D, J


def _bdj_numpy(phi, m, n, compute_J):
    s  = np.sin(phi)
    c  = np.cos(phi)
    d2 = 1.0 - m * s ** 2
    s3o3 = s ** 3 / 3.0

    RF = _rf_numpy(c ** 2, d2, np.ones_like(phi))
    RD = _rd_numpy(c ** 2, d2, np.ones_like(phi))

    F_val = s * RF
    D_val = s3o3 * RD
    B_val = F_val - D_val

    zero = s == 0.0
    B_val[zero] = 0.0
    D_val[zero] = 0.0

    if compute_J:
        p = 1.0 - n * s ** 2
        RJ = _rj_numpy(c ** 2, d2, np.ones_like(phi), p)
        J_val = s3o3 * RJ
        J_val[zero] = 0.0
    else:
        J_val = None

    return B_val, D_val, J_val
