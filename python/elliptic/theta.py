"""Jacobi theta functions, their derivatives, and jacobiThetaEta.

Conventions follow Abramowitz & Stegun §16–17:
    θ₁(v, q) = 2 Σ_{n=0}^∞ (-1)^n q^{(n+1/2)^2} sin((2n+1)v)
    θ₂(v, q) = 2 Σ_{n=0}^∞ q^{(n+1/2)^2} cos((2n+1)v)
    θ₃(v, q) = 1 + 2 Σ_{n=1}^∞ q^{n^2} cos(2nv)
    θ₄(v, q) = 1 + 2 Σ_{n=1}^∞ (-1)^n q^{n^2} cos(2nv)

Jacobi theta Θ(u|m) = θ₄(πu/(2K), q)  [A&S 16.27]
Jacobi eta   H(u|m) = θ₁(πu/(2K), q)  [A&S 16.27]

The public `theta(j, v, m)` function takes v in radians (angle argument)
and maps to the standard θⱼ(v, q(m)).
"""
from __future__ import annotations
import math
from ._xputils import get_xp
from .carlson import _rf_xp

_N_TERMS = 30   # sufficient for |q| ≤ 0.8 (m ≤ ~0.9997)


def _q_from_m_xp(xp, m):
    """Backend-native nome q(m), with invalid/singular inputs masked."""
    valid = (m >= 0.0) & (m < 1.0)
    m_safe = xp.where(valid, m, xp.full_like(m, 0.5))
    zero = xp.zeros_like(m_safe)
    one = xp.ones_like(m_safe)
    K = _rf_xp(xp, zero, 1.0 - m_safe, one)
    Kp = _rf_xp(xp, zero, m_safe, one)
    q = xp.exp(-math.pi * Kp / K)
    return xp.where(valid, q, xp.full_like(q, math.nan))


# -----------------------------------------------------------------------
# Low-level series (flat 1-D numpy, v in radians, q scalar or array)
# -----------------------------------------------------------------------

def _th1(xp, v, q):
    """θ₁(v, q)."""
    s = xp.zeros_like(v)
    for n in range(_N_TERMS):
        s = s + (-1)**n * q**((n + 0.5)**2) * xp.sin((2*n + 1) * v)
    return 2.0 * s


def _th2(xp, v, q):
    """θ₂(v, q)."""
    s = xp.zeros_like(v)
    for n in range(_N_TERMS):
        s = s + q**((n + 0.5)**2) * xp.cos((2*n + 1) * v)
    return 2.0 * s


def _th3(xp, v, q):
    """θ₃(v, q)."""
    s = xp.ones_like(v)
    for n in range(1, _N_TERMS + 1):
        s = s + 2.0 * q**(n**2) * xp.cos(2*n * v)
    return s


def _th4(xp, v, q):
    """θ₄(v, q)."""
    s = xp.ones_like(v)
    for n in range(1, _N_TERMS + 1):
        s = s + 2.0 * (-1)**n * q**(n**2) * xp.cos(2*n * v)
    return s


# Derivatives dθⱼ/dv
def _dth1(xp, v, q):
    s = xp.zeros_like(v)
    for n in range(_N_TERMS):
        s = s + (-1)**n * (2*n+1) * q**((n + 0.5)**2) * xp.cos((2*n + 1) * v)
    return 2.0 * s


def _dth2(xp, v, q):
    s = xp.zeros_like(v)
    for n in range(_N_TERMS):
        s = s - (2*n+1) * q**((n + 0.5)**2) * xp.sin((2*n + 1) * v)
    return 2.0 * s


def _dth3(xp, v, q):
    s = xp.zeros_like(v)
    for n in range(1, _N_TERMS + 1):
        s = s - 4.0 * n * q**(n**2) * xp.sin(2*n * v)
    return s


def _dth4(xp, v, q):
    s = xp.zeros_like(v)
    for n in range(1, _N_TERMS + 1):
        s = s - 4.0 * n * (-1)**n * q**(n**2) * xp.sin(2*n * v)
    return s


_TH_FNS  = [None, _th1,  _th2,  _th3,  _th4]
_DTH_FNS = [None, _dth1, _dth2, _dth3, _dth4]


# -----------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------

def jacobiThetaEta(u, m):
    """Jacobi theta Θ(u|m) and eta H(u|m).

    Θ(u|m) = θ₄(πu/(2K(m)), q(m))   [A&S 16.27]
    H(u|m) = θ₁(πu/(2K(m)), q(m))   [A&S 16.27]

    Parameters
    ----------
    u : array_like
        Argument in Jacobi scale (same unit as K(m)).
    m : array_like
        Parameter, 0 <= m <= 1.

    Returns
    -------
    Th, H : arrays
        Jacobi theta and eta values.
    """
    xp = get_xp(u, m)
    u = xp.asarray(u, dtype=xp.float64)
    m = xp.asarray(m, dtype=xp.float64)
    u, m = xp.broadcast_arrays(u, m)
    valid = (m >= 0.0) & (m < 1.0)
    m_safe = xp.where(valid, m, xp.full_like(m, 0.5))
    zero = xp.zeros_like(m_safe)
    K = _rf_xp(xp, zero, 1.0 - m_safe, xp.ones_like(m_safe))
    q = _q_from_m_xp(xp, m)
    v = math.pi * u / (2.0 * K)
    Th = _th4(xp, v, q)
    H = _th1(xp, v, q)
    nan = xp.full_like(Th, math.nan)
    return xp.where(valid, Th, nan), xp.where(valid, H, nan)


def theta(j, v, m):
    """Jacobi theta function of type j = 1, 2, 3, or 4.

    Takes the angle argument v in radians and parameter m in [0, 1].

    θ₁(v, q(m)) — odd, zeros at 0 and nπ
    θ₂(v, q(m)) — zeros at (n+1/2)π
    θ₃(v, q(m)) — zeros at (n+1/2)π + K'τ
    θ₄(v, q(m)) — zeros at nπ + K'τ   [= Jacobi Θ with argument rescaling]

    Parameters
    ----------
    j : int  (1–4)
    v : array_like  angle in radians
    m : array_like  parameter 0 <= m <= 1

    Returns
    -------
    Th : array
    """
    if j not in (1, 2, 3, 4):
        raise ValueError("j must be 1, 2, 3, or 4")

    xp = get_xp(v, m)
    v = xp.asarray(v, dtype=xp.float64)
    m = xp.asarray(m, dtype=xp.float64)
    v, m = xp.broadcast_arrays(v, m)
    q = _q_from_m_xp(xp, m)
    return _TH_FNS[j](xp, v, q)


def theta_prime(j, v, m):
    """Jacobi theta function and its derivative with respect to v.

    Parameters
    ----------
    j : int  (1–4)
    v : array_like  angle in radians
    m : array_like  parameter 0 <= m <= 1

    Returns
    -------
    th  : array  theta_j(v, m)
    thp : array  d/dv theta_j(v, m)
    """
    if j not in (1, 2, 3, 4):
        raise ValueError("j must be 1, 2, 3, or 4")

    xp = get_xp(v, m)
    v = xp.asarray(v, dtype=xp.float64)
    m = xp.asarray(m, dtype=xp.float64)
    v, m = xp.broadcast_arrays(v, m)
    q = _q_from_m_xp(xp, m)
    return _TH_FNS[j](xp, v, q), _DTH_FNS[j](xp, v, q)
