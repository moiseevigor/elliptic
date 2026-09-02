"""Carlson symmetric elliptic integrals RF, RD, RJ, RC.

All use Carlson's duplication algorithm with fixed iteration counts
(20 for RF, 30 for RD, 60 for RJ) so they are JAX-traceable and run natively on
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
    diff = y - x
    # Branch selection must be RELATIVE: R_C is homogeneous of degree -1/2
    # (DLMF 19.20.3), and an absolute |y-x| < 1e-14 window sent every
    # small-scale input down the degenerate x==y branch -- RC(1e-20, 2e-20)
    # returned 1/sqrt(x), a 27% error.
    scale = xp.maximum(xp.abs(x), xp.abs(y))
    tol = 1e-14 * scale

    # safe arguments for each branch (avoid div-by-zero when not selected)
    x_safe   = xp.where(x > 0, x, xp.full_like(x, 1.0))
    yd_safe  = xp.where(diff > 0, diff, xp.full_like(diff, 1.0))
    yd_safe2 = xp.where(-diff > 0, -diff, xp.full_like(diff, 1.0))
    y_safe   = xp.where(y > 0, y, xp.full_like(y, 1.0))

    rc_gt  = xp.arctan(xp.sqrt(xp.clip(diff / x_safe, 0.0, None))) / xp.sqrt(yd_safe)
    lt_active = (diff < -tol) & (y > 0) & (x > 0)
    lt_ratio = xp.where(lt_active, -diff / x_safe, xp.full_like(diff, 0.5))
    # x > y: log((sqrt(x)+sqrt(x-y))/sqrt(y))/sqrt(x-y) -- algebraically
    # arctanh(sqrt(1-y/x))/sqrt(x-y), but without the 1 - sqrt(1-eps)
    # cancellation that lost 8 digits at RC(3, 1e-10).
    lt_y   = xp.where(lt_active, y, xp.ones_like(y))
    lt_x   = xp.where(lt_active, x, xp.full_like(x, 2.0))
    # ... and as log1p, since log((sx+sxy)/sy) = log1p(((x-y)/(sx+sy) + sxy)/sy):
    # the plain log lost 9 digits again when x - y was tiny (RC(1+1e-13, 1)).
    sx, sy = xp.sqrt(lt_x), xp.sqrt(lt_y)
    sxy    = xp.sqrt(xp.clip(lt_x - lt_y, 0.0, None))
    rc_lt  = xp.log1p(((lt_x - lt_y) / (sx + sy) + sxy) / sy) / xp.sqrt(yd_safe2)
    rc_eq  = 1.0 / xp.sqrt(x_safe)
    rc_x0  = (math.pi * 0.5) / xp.sqrt(y_safe)

    out = xp.where(diff > tol, rc_gt, xp.where(diff < -tol, rc_lt, rc_eq))
    out = xp.where(x == 0, rc_x0, out)
    out = xp.where(y == 0, xp.full_like(out, math.inf), out)
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
    out = _rf_xp(xp, x, y, z)
    # Two zero arguments: the integral diverges (DLMF 19.16.1); the fixed
    # duplication count otherwise returns a finite number.
    # pure boolean algebra: torch tensors have no .astype, and this must
    # stay backend-native (caught on the L4 hardware run)
    two0 = ((x == 0) & (y == 0)) | ((x == 0) & (z == 0)) | ((y == 0) & (z == 0))
    return xp.where(two0, xp.full_like(out, math.inf), out)


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
    out = _rd_xp(xp, x, y, z)
    return xp.where((x == 0) & (y == 0), xp.full_like(out, math.inf), out)  # DLMF 19.16.5


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

    Requires p > 0. For p < 0 the integral is a Cauchy principal value
    (DLMF 19.20.14); not implemented — raises ValueError.
    """
    xp = get_xp(x, y, z, p)
    x = xp.asarray(x, dtype=xp.float64)
    y = xp.asarray(y, dtype=xp.float64)
    z = xp.asarray(z, dtype=xp.float64)
    p = xp.asarray(p, dtype=xp.float64)
    x, y, z, p = xp.broadcast_arrays(x, y, z, p)

    if xp is np and np.any(p <= 0.0):
        raise ValueError(
            "carlsonRJ: p must be > 0. For p < 0 the integral is a Cauchy "
            "principal value (DLMF 19.20.14); use the transformation to "
            "a q > 0 argument before calling."
        )
    out = _rj_xp(xp, x, y, z, p)
    # pure boolean algebra: torch tensors have no .astype, and this must
    # stay backend-native (caught on the L4 hardware run)
    two0 = ((x == 0) & (y == 0)) | ((x == 0) & (z == 0)) | ((y == 0) & (z == 0))
    return xp.where(two0, xp.full_like(out, math.inf), out)  # DLMF 19.16.2


def _rj_xp(xp, x, y, z, p):
    S   = xp.zeros_like(x)
    fac = xp.ones_like(x)
    # 60 duplications: each halves the argument-ratio exponent (base 4), so
    # the series is valid for max/min argument ratios up to ~4^54 = 3e32.
    # 30 covered only ~1e16 -- RJ(1e-20, 2e-20, 3e-20, 0.5) was 11% off.
    for _ in range(60):
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
