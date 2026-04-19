"""Carlson symmetric elliptic integrals RF, RD, RJ, RC.

All use Carlson's duplication algorithm with fixed iteration counts
(20 for RF, 30 for RD/RJ) so they are JAX-traceable and run natively on
any array backend (NumPy, PyTorch CUDA, JAX).

References
----------
NIST DLMF §19.36 — https://dlmf.nist.gov/19.36
B.C. Carlson, Numer. Algorithms 10 (1995), 13–26.
"""
from __future__ import annotations

import math
import numpy as np

from ._xputils import get_xp


# ---------------------------------------------------------------------------
# RC — degenerate RF(x, y, y)  (closed-form)
# ---------------------------------------------------------------------------

def carlsonRC(x, y):
    """Carlson RC(x, y) = RF(x, y, y).  Closed-form (DLMF 19.2.17)."""
    xp = get_xp(x, y)
    x = xp.asarray(x, dtype=xp.float64)
    y = xp.asarray(y, dtype=xp.float64)
    x, y = xp.broadcast_arrays(x, y)
    return _rc_xp(xp, x, y)


def _rc_xp(xp, x, y):
    EPS = 1e-300
    diff = y - x

    # safe arguments for each branch (avoid div-by-zero when not selected)
    x_safe   = xp.where(x > EPS, x, xp.full_like(x, 1.0))
    yd_safe  = xp.where(diff > EPS, diff, xp.full_like(diff, 1.0))
    yd_safe2 = xp.where(-diff > EPS, -diff, xp.full_like(diff, 1.0))
    y_safe   = xp.where(y > EPS, y, xp.full_like(y, 1.0))

    rc_gt  = xp.arctan(xp.sqrt(xp.clip(diff / x_safe, 0.0, None))) / xp.sqrt(yd_safe)
    rc_lt  = xp.arctanh(xp.sqrt(xp.clip(-diff / x_safe, 0.0, None))) / xp.sqrt(yd_safe2)
    rc_eq  = 1.0 / xp.sqrt(x_safe)
    rc_x0  = (math.pi * 0.5) / xp.sqrt(y_safe)

    TOL = 1e-14
    out = xp.where(diff > TOL, rc_gt, xp.where(diff < -TOL, rc_lt, rc_eq))
    out = xp.where(x < EPS, rc_x0, out)
    out = xp.where(y < EPS, xp.full_like(out, math.inf), out)
    return out


# keep old numpy version for any legacy callers
def _rc_numpy(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return _rc_xp(np, x, y)


# ---------------------------------------------------------------------------
# RF — symmetric first kind
# ---------------------------------------------------------------------------

def carlsonRF(x, y, z):
    """Carlson RF(x, y, z) — symmetric elliptic integral of the first kind.

    F(phi|m) = sin(phi) * RF(cos²phi, 1 - m sin²phi, 1).
    """
    xp = get_xp(x, y, z)
    x = xp.asarray(x, dtype=xp.float64)
    y = xp.asarray(y, dtype=xp.float64)
    z = xp.asarray(z, dtype=xp.float64)
    x, y, z = xp.broadcast_arrays(x, y, z)
    return _rf_xp(xp, x, y, z)


def _rf_xp(xp, x, y, z):
    for _ in range(20):
        lam = xp.sqrt(x * y) + xp.sqrt(y * z) + xp.sqrt(z * x)
        x = (x + lam) * 0.25
        y = (y + lam) * 0.25
        z = (z + lam) * 0.25
    A  = (x + y + z) / 3.0
    X  = (A - x) / A
    Y  = (A - y) / A
    Z  = -X - Y
    E2 = X * Y - Z * Z
    E3 = X * Y * Z
    return A ** (-0.5) * (1.0 - E2/10.0 + E3/14.0 + E2**2/24.0 - 3.0*E2*E3/44.0)


def _rf_numpy(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    return _rf_xp(np, x, y, z)


# ---------------------------------------------------------------------------
# RD — symmetric second kind
# ---------------------------------------------------------------------------

def carlsonRD(x, y, z):
    """Carlson RD(x, y, z) — symmetric elliptic integral of the second kind.

    D(phi|m) = (sin³phi/3) * RD(cos²phi, 1 - m sin²phi, 1).
    """
    xp = get_xp(x, y, z)
    x = xp.asarray(x, dtype=xp.float64)
    y = xp.asarray(y, dtype=xp.float64)
    z = xp.asarray(z, dtype=xp.float64)
    x, y, z = xp.broadcast_arrays(x, y, z)
    return _rd_xp(xp, x, y, z)


def _rd_xp(xp, x, y, z):
    S   = xp.zeros_like(x)
    fac = xp.ones_like(x)
    for _ in range(30):
        lam = xp.sqrt(x * y) + xp.sqrt(y * z) + xp.sqrt(z * x)
        S   = S + fac / (xp.sqrt(z) * (z + lam))
        fac = fac * 0.25
        x = (x + lam) * 0.25
        y = (y + lam) * 0.25
        z = (z + lam) * 0.25
    A  = (x + y + 3.0 * z) / 5.0
    X  = (A - x) / A
    Y  = (A - y) / A
    Z  = -(X + Y) / 3.0
    E2 = X * Y - 6.0 * Z**2
    E3 = (3.0 * X * Y - 8.0 * Z**2) * Z
    E4 = 3.0 * (X * Y - Z**2) * Z**2
    E5 = X * Y * Z**3
    poly = (1.0 - 3.0*E2/14.0 + E3/6.0 + 9.0*E2**2/88.0
            - 3.0*E4/22.0 - 9.0*E2*E3/52.0 + 3.0*E5/26.0)
    return 3.0 * S + fac * A**(-1.5) * poly


def _rd_numpy(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    return _rd_xp(np, x, y, z)


# ---------------------------------------------------------------------------
# RJ — symmetric third kind
# ---------------------------------------------------------------------------

def carlsonRJ(x, y, z, p):
    """Carlson RJ(x, y, z, p) — symmetric elliptic integral of the third kind.

    J(phi,n|m) = (sin³phi/3) * RJ(cos²phi, 1-m sin²phi, 1, 1-n sin²phi).
    """
    xp = get_xp(x, y, z, p)
    x = xp.asarray(x, dtype=xp.float64)
    y = xp.asarray(y, dtype=xp.float64)
    z = xp.asarray(z, dtype=xp.float64)
    p = xp.asarray(p, dtype=xp.float64)
    x, y, z, p = xp.broadcast_arrays(x, y, z, p)
    return _rj_xp(xp, x, y, z, p)


def _rj_xp(xp, x, y, z, p):
    S   = xp.zeros_like(x)
    fac = xp.ones_like(x)
    for _ in range(30):
        lam   = xp.sqrt(x * y) + xp.sqrt(y * z) + xp.sqrt(z * x)
        alpha = (p * (xp.sqrt(x) + xp.sqrt(y) + xp.sqrt(z)) + xp.sqrt(x * y * z)) ** 2
        beta  = p * (p + lam) ** 2
        S     = S + fac * _rc_xp(xp, alpha, beta)
        fac   = fac * 0.25
        x = (x + lam) * 0.25
        y = (y + lam) * 0.25
        z = (z + lam) * 0.25
        p = (p + lam) * 0.25
    A  = (x + y + z + 2.0 * p) / 5.0
    X  = (A - x) / A
    Y  = (A - y) / A
    Z  = (A - z) / A
    P  = -(X + Y + Z) / 2.0
    E2 = X*Y + X*Z + Y*Z - 3.0*P**2
    E3 = X*Y*Z + 2.0*E2*P + 3.0*P**3
    E4 = (2.0*X*Y*Z + E2*P + 3.0*P**3) * P
    E5 = X*Y*Z * P**2
    poly = (1.0 - 3.0*E2/14.0 + E3/6.0 + 9.0*E2**2/88.0
            - 3.0*E4/22.0 - 9.0*E2*E3/52.0 + 3.0*E5/26.0)
    return 3.0 * S + fac * A**(-1.5) * poly


def _rj_numpy(x: np.ndarray, y: np.ndarray, z: np.ndarray, p: np.ndarray) -> np.ndarray:
    return _rj_xp(np, x, y, z, p)
