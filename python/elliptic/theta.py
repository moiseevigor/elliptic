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
import numpy as np
from array_api_compat import array_namespace

_N_TERMS = 30   # sufficient for |q| ≤ 0.8 (m ≤ ~0.9997)


def _q_from_m(m_np: np.ndarray) -> np.ndarray:
    """Nome q(m) — uses scipy K; kept internal to avoid circular import."""
    from scipy.special import ellipk
    K  = ellipk(m_np)
    Kp = ellipk(1.0 - m_np)
    return np.exp(-np.pi * Kp / K)


# -----------------------------------------------------------------------
# Low-level series (flat 1-D numpy, v in radians, q scalar or array)
# -----------------------------------------------------------------------

def _th1(v: np.ndarray, q: np.ndarray) -> np.ndarray:
    """θ₁(v, q)."""
    s = np.zeros_like(v)
    for n in range(_N_TERMS):
        s += (-1)**n * q**((n + 0.5)**2) * np.sin((2*n + 1) * v)
    return 2.0 * s


def _th2(v: np.ndarray, q: np.ndarray) -> np.ndarray:
    """θ₂(v, q)."""
    s = np.zeros_like(v)
    for n in range(_N_TERMS):
        s += q**((n + 0.5)**2) * np.cos((2*n + 1) * v)
    return 2.0 * s


def _th3(v: np.ndarray, q: np.ndarray) -> np.ndarray:
    """θ₃(v, q)."""
    s = np.ones_like(v)
    for n in range(1, _N_TERMS + 1):
        s += 2.0 * q**(n**2) * np.cos(2*n * v)
    return s


def _th4(v: np.ndarray, q: np.ndarray) -> np.ndarray:
    """θ₄(v, q)."""
    s = np.ones_like(v)
    for n in range(1, _N_TERMS + 1):
        s += 2.0 * (-1)**n * q**(n**2) * np.cos(2*n * v)
    return s


# Derivatives dθⱼ/dv
def _dth1(v: np.ndarray, q: np.ndarray) -> np.ndarray:
    s = np.zeros_like(v)
    for n in range(_N_TERMS):
        s += (-1)**n * (2*n+1) * q**((n + 0.5)**2) * np.cos((2*n + 1) * v)
    return 2.0 * s


def _dth2(v: np.ndarray, q: np.ndarray) -> np.ndarray:
    s = np.zeros_like(v)
    for n in range(_N_TERMS):
        s += -(2*n+1) * q**((n + 0.5)**2) * np.sin((2*n + 1) * v)
    return 2.0 * s


def _dth3(v: np.ndarray, q: np.ndarray) -> np.ndarray:
    s = np.zeros_like(v)
    for n in range(1, _N_TERMS + 1):
        s += -2.0 * 2*n * q**(n**2) * np.sin(2*n * v)
    return s


def _dth4(v: np.ndarray, q: np.ndarray) -> np.ndarray:
    s = np.zeros_like(v)
    for n in range(1, _N_TERMS + 1):
        s += 2.0 * (-1)**n * (-2*n) * q**(n**2) * np.sin(2*n * v)
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
    from scipy.special import ellipk

    u = np.asarray(u, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64)
    xp = array_namespace(u, m)
    u, m = xp.broadcast_arrays(u, m)
    orig_shape = np.asarray(u).shape

    u_np = np.asarray(u).ravel()
    m_np = np.asarray(m).ravel()

    Th = np.ones_like(u_np)
    H  = np.zeros_like(u_np)

    # m = 1: undefined
    mask1 = m_np >= 1.0 - 1e-14
    Th[mask1] = np.nan
    H[mask1]  = np.nan

    # m = 0: Th = 1, H = 0
    mask0 = m_np < 1e-14
    # already set by np.ones / np.zeros

    maskN = ~mask0 & ~mask1
    if np.any(maskN):
        u_g = u_np[maskN]
        m_g = m_np[maskN]
        K_g = ellipk(m_g)
        q_g = _q_from_m(m_g)
        v_g = np.pi * u_g / (2.0 * K_g)       # normalised angle argument
        Th[maskN] = _th4(v_g, q_g)
        H[maskN]  = _th1(v_g, q_g)

    Th = Th.reshape(orig_shape)
    H  = H.reshape(orig_shape)
    return xp.asarray(Th), xp.asarray(H)


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

    v = np.asarray(v, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64)
    xp = array_namespace(v, m)
    v, m = xp.broadcast_arrays(v, m)
    orig_shape = np.asarray(v).shape

    v_np = np.asarray(v).ravel()
    m_np = np.asarray(m).ravel()

    Th = np.zeros_like(v_np)
    maskN = (m_np >= 1e-14) & (m_np < 1.0 - 1e-14)

    if np.any(maskN):
        q_g  = _q_from_m(m_np[maskN])
        Th[maskN] = _TH_FNS[j](v_np[maskN], q_g)

    # special cases
    mask0 = m_np < 1e-14
    if np.any(mask0):
        if j == 1:   Th[mask0] = np.sin(v_np[mask0])  # θ₁(v,0) = sin v (leading term)
        elif j == 2: Th[mask0] = 0.0                   # q→0 all terms vanish except trivially 0
        elif j == 3: Th[mask0] = 1.0
        elif j == 4: Th[mask0] = 1.0

    mask1 = m_np >= 1.0 - 1e-14
    Th[mask1] = np.nan

    return xp.asarray(Th.reshape(orig_shape))


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

    v = np.asarray(v, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64)
    xp = array_namespace(v, m)
    v, m = xp.broadcast_arrays(v, m)
    orig_shape = np.asarray(v).shape

    v_np = np.asarray(v).ravel()
    m_np = np.asarray(m).ravel()

    th_np  = np.zeros_like(v_np)
    thp_np = np.zeros_like(v_np)

    maskN = (m_np >= 1e-14) & (m_np < 1.0 - 1e-14)
    if np.any(maskN):
        q_g = _q_from_m(m_np[maskN])
        th_np[maskN]  = _TH_FNS[j](v_np[maskN], q_g)
        thp_np[maskN] = _DTH_FNS[j](v_np[maskN], q_g)

    th  = th_np.reshape(orig_shape)
    thp = thp_np.reshape(orig_shape)
    return xp.asarray(th), xp.asarray(thp)
