"""Inverse incomplete elliptic integral of the second kind."""
from __future__ import annotations
import math
import numpy as np
from ._xputils import get_xp, is_numpy
from .elliptic12 import _elliptic12_xp


def inverselliptic2(E_val, m, tol=1e-12):
    """Inverse of the incomplete elliptic integral of the second kind.

    Solves E(phi | m) = E_val using period reduction and fixed-step Newton
    refinement against the library's own Carlson implementation.

    Period identity:  E(phi + pi | m) = E(phi | m) + 2*E(m)
    Symmetry:         E(pi - phi | m) = 2*E(m) - E(phi | m)

    Parameters
    ----------
    E_val : array_like   Target value(s) of E(phi | m).
    m     : array_like   Parameter, 0 <= m <= 1.
    tol   : float        Newton convergence tolerance (default 1e-12).

    Returns
    -------
    phi : array   Amplitude in radians such that E(phi | m) ≈ E_val.
    """
    xp = get_xp(E_val, m)
    E_val = xp.asarray(E_val, dtype=xp.float64)
    m = xp.asarray(m, dtype=xp.float64)
    E_val, m = xp.broadcast_arrays(E_val, m)

    if is_numpy(xp) and np.any((m < 0.0) | (m > 1.0)):
        raise ValueError("m must be in [0, 1]")

    # Complete integral E(m); each phi-period of π contributes 2*E1 to E.
    half_pi = xp.full_like(m, math.pi * 0.5)
    _, E1, _ = _elliptic12_xp(xp, half_pi, m)
    two_E1 = 2.0 * E1

    # Oddness first: E(-phi) = -E(phi).  Folding a tiny negative z through
    # 2*E1 - (z + 2*E1) lost all its digits (rel 1e-7 at z = -1e-9*E1).
    sgn   = xp.where(E_val < 0.0, -xp.ones_like(E_val), xp.ones_like(E_val))
    E_val = xp.abs(E_val)

    # Step 1 — strip full periods:  phi = phi_base + k*pi
    k = xp.floor(E_val / two_E1)
    z_red = E_val - k * two_E1          # in [0, 2*E1)

    # Step 2 — fold second half-period using E(pi-phi|m) = 2E1 - E(phi|m)
    over = z_red > E1
    z_red2 = xp.where(over, two_E1 - z_red, z_red)   # in [0, E1]

    # Monotone linear seed in [0, pi/2].  Fixed iteration count keeps the
    # routine JIT-safe; converged elements simply receive zero-sized updates.
    phi = xp.clip((z_red2 / E1) * (math.pi * 0.5), 0.0, math.pi * 0.5)

    for _ in range(24):
        _, E_cur, _ = _elliptic12_xp(xp, phi, m)
        res = E_cur - z_red2
        denom = xp.sqrt(xp.clip(1.0 - m * xp.sin(phi) ** 2, 0.0, None))
        safe_denom = xp.where(denom > tol, denom, xp.ones_like(denom))
        # unconditional step: gating on |res| > tol froze the solver at an
        # absolute 1e-12, i.e. only ~3 relative digits for |z| ~ 1e-9
        step = res / safe_denom
        phi = xp.clip(phi - step, 0.0, math.pi * 0.5)

    # Step 3 — undo fold: phi_in_period = pi - phi (if over), else phi
    phi = xp.where(over, math.pi - phi, phi)   # in [0, pi)

    # Step 4 — undo period strips and the sign
    return sgn * (phi + k * math.pi)
