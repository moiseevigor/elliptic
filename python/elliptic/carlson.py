"""Carlson symmetric elliptic integrals RF, RD, RJ, RC.

All use Carlson's duplication algorithm with fixed iteration counts
(20 for RF, 30 for RD/RJ) so they are JAX-traceable.

References
----------
NIST DLMF §19.36 — https://dlmf.nist.gov/19.36
B.C. Carlson, Numer. Algorithms 10 (1995), 13-26.
"""
from __future__ import annotations

import numpy as np
from array_api_compat import array_namespace


# ---------------------------------------------------------------------------
# RC — degenerate case RF(x, y, y)  (closed-form)
# ---------------------------------------------------------------------------

def carlsonRC(x, y):
    """Carlson RC(x, y) = RF(x, y, y).

    Closed-form (DLMF 19.2.17-19.2.19).
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    xp = array_namespace(x, y)
    x, y = xp.broadcast_arrays(x, y)
    orig_shape = x.shape
    rc = _rc_numpy(np.asarray(x).ravel(), np.asarray(y).ravel())
    return xp.asarray(rc.reshape(orig_shape))


def _rc_numpy(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    N = x.size
    rc = np.zeros(N)
    pole   = y == 0.0
    zero_x = (x == 0.0) & ~pole
    eq     = (x == y) & ~pole & ~zero_x
    gt     = (y > x) & ~pole & ~zero_x & ~eq
    lt     = (y < x) & ~pole & ~zero_x & ~eq

    rc[pole]   = np.inf
    rc[zero_x] = np.pi / (2.0 * np.sqrt(y[zero_x]))
    rc[eq]     = 1.0 / np.sqrt(x[eq])

    if np.any(gt):
        d = np.sqrt((y[gt] - x[gt]) / x[gt])
        rc[gt] = np.arctan(d) / np.sqrt(y[gt] - x[gt])
    if np.any(lt):
        d = np.sqrt((x[lt] - y[lt]) / x[lt])
        rc[lt] = np.arctanh(d) / np.sqrt(x[lt] - y[lt])
    return rc


# ---------------------------------------------------------------------------
# RF — symmetric first kind
# ---------------------------------------------------------------------------

def carlsonRF(x, y, z):
    """Carlson RF(x, y, z) — symmetric elliptic integral of the first kind.

    Connection: F(phi|m) = sin(phi) * RF(cos^2 phi, 1 - m sin^2 phi, 1).
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    xp = array_namespace(x, y, z)
    x, y, z = xp.broadcast_arrays(x, y, z)
    orig_shape = x.shape
    rf = _rf_numpy(np.asarray(x).ravel(), np.asarray(y).ravel(), np.asarray(z).ravel())
    return xp.asarray(rf.reshape(orig_shape))


def _rf_numpy(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    for _ in range(20):
        lam = np.sqrt(x * y) + np.sqrt(y * z) + np.sqrt(z * x)
        x = (x + lam) / 4.0
        y = (y + lam) / 4.0
        z = (z + lam) / 4.0

    A  = (x + y + z) / 3.0
    X  = (A - x) / A
    Y  = (A - y) / A
    Z  = -X - Y
    E2 = X * Y - Z ** 2
    E3 = X * Y * Z
    return A ** (-0.5) * (1.0 - E2 / 10.0 + E3 / 14.0 + E2 ** 2 / 24.0 - 3.0 * E2 * E3 / 44.0)


# ---------------------------------------------------------------------------
# RD — symmetric second kind (degenerate RJ(x,y,z,z))
# ---------------------------------------------------------------------------

def carlsonRD(x, y, z):
    """Carlson RD(x, y, z) — symmetric elliptic integral of the second kind.

    Connection: D(phi|m) = (sin^3 phi / 3) * RD(cos^2 phi, 1 - m sin^2 phi, 1).
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    xp = array_namespace(x, y, z)
    x, y, z = xp.broadcast_arrays(x, y, z)
    orig_shape = x.shape
    rd = _rd_numpy(np.asarray(x).ravel(), np.asarray(y).ravel(), np.asarray(z).ravel())
    return xp.asarray(rd.reshape(orig_shape))


def _rd_numpy(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    S   = np.zeros_like(x)
    fac = np.ones_like(x)
    for _ in range(30):
        lam = np.sqrt(x * y) + np.sqrt(y * z) + np.sqrt(z * x)
        S  += fac / (np.sqrt(z) * (z + lam))
        fac /= 4.0
        x = (x + lam) / 4.0
        y = (y + lam) / 4.0
        z = (z + lam) / 4.0

    A  = (x + y + 3.0 * z) / 5.0
    X  = (A - x) / A
    Y  = (A - y) / A
    Z  = -(X + Y) / 3.0
    E2 = X * Y - 6.0 * Z ** 2
    E3 = (3.0 * X * Y - 8.0 * Z ** 2) * Z
    E4 = 3.0 * (X * Y - Z ** 2) * Z ** 2
    E5 = X * Y * Z ** 3
    poly = (1.0 - 3.0 * E2 / 14.0 + E3 / 6.0 + 9.0 * E2 ** 2 / 88.0
            - 3.0 * E4 / 22.0 - 9.0 * E2 * E3 / 52.0 + 3.0 * E5 / 26.0)
    return 3.0 * S + fac * A ** (-1.5) * poly


# ---------------------------------------------------------------------------
# RJ — symmetric third kind
# ---------------------------------------------------------------------------

def carlsonRJ(x, y, z, p):
    """Carlson RJ(x, y, z, p) — symmetric elliptic integral of the third kind.

    Connection: J(phi,n|m) = (sin^3 phi / 3) * RJ(cos^2 phi, 1-m sin^2 phi, 1, 1-n sin^2 phi).
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)
    xp = array_namespace(x, y, z, p)
    x, y, z, p = xp.broadcast_arrays(x, y, z, p)
    orig_shape = x.shape
    rj = _rj_numpy(np.asarray(x).ravel(), np.asarray(y).ravel(),
                   np.asarray(z).ravel(), np.asarray(p).ravel())
    return xp.asarray(rj.reshape(orig_shape))


def _rj_numpy(x: np.ndarray, y: np.ndarray, z: np.ndarray, p: np.ndarray) -> np.ndarray:
    S   = np.zeros_like(x)
    fac = np.ones_like(x)
    for _ in range(30):
        lam   = np.sqrt(x * y) + np.sqrt(y * z) + np.sqrt(z * x)
        alpha = (p * (np.sqrt(x) + np.sqrt(y) + np.sqrt(z)) + np.sqrt(x * y * z)) ** 2
        beta  = p * (p + lam) ** 2
        S    += fac * _rc_numpy(alpha, beta)
        fac  /= 4.0
        x = (x + lam) / 4.0
        y = (y + lam) / 4.0
        z = (z + lam) / 4.0
        p = (p + lam) / 4.0

    A  = (x + y + z + 2.0 * p) / 5.0
    X  = (A - x) / A
    Y  = (A - y) / A
    Z  = (A - z) / A
    P  = -(X + Y + Z) / 2.0
    E2 = X * Y + X * Z + Y * Z - 3.0 * P ** 2
    E3 = X * Y * Z + 2.0 * E2 * P + 3.0 * P ** 3
    E4 = (2.0 * X * Y * Z + E2 * P + 3.0 * P ** 3) * P
    E5 = X * Y * Z * P ** 2
    poly = (1.0 - 3.0 * E2 / 14.0 + E3 / 6.0 + 9.0 * E2 ** 2 / 88.0
            - 3.0 * E4 / 22.0 - 9.0 * E2 * E3 / 52.0 + 3.0 * E5 / 26.0)
    return 3.0 * S + fac * A ** (-1.5) * poly
