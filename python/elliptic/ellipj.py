"""Jacobi elliptic functions sn(u|m), cn(u|m), dn(u|m) and amplitude am(u|m).

Algorithm: Arithmetic-Geometric Mean + descending Landen back-substitution
(Abramowitz & Stegun §16.4).  Fixed _AGM_ITERS = 25 iterations makes the
loop JAX-traceable; convergence always occurs within 25 steps for float64.
"""
from __future__ import annotations

import numpy as np
from array_api_compat import array_namespace

_AGM_ITERS = 25


def ellipj(u, m):
    """Jacobi elliptic functions.

    Parameters
    ----------
    u : array_like
        Argument.
    m : array_like
        Parameter, 0 <= m <= 1.

    Returns
    -------
    sn, cn, dn, am : arrays (broadcast shape of u and m)
    """
    u = np.asarray(u, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64)
    xp = array_namespace(u, m)
    u, m = xp.broadcast_arrays(u, m)
    orig_shape = u.shape

    sn, cn, dn, am = _ellipj_numpy(np.asarray(u).ravel(), np.asarray(m).ravel())
    return (xp.asarray(sn.reshape(orig_shape)),
            xp.asarray(cn.reshape(orig_shape)),
            xp.asarray(dn.reshape(orig_shape)),
            xp.asarray(am.reshape(orig_shape)))


def _ellipj_numpy(u: np.ndarray, m: np.ndarray):
    N = u.size
    sn = np.zeros(N)
    cn = np.zeros(N)
    dn = np.zeros(N)
    am = np.zeros(N)

    eps = np.finfo(np.float64).eps

    mask0 = m == 0.0
    mask1 = m == 1.0
    maskN = ~mask0 & ~mask1 & (m >= 0.0) & (m <= 1.0)

    # m == 0
    am[mask0] = u[mask0]
    sn[mask0] = np.sin(u[mask0])
    cn[mask0] = np.cos(u[mask0])
    dn[mask0] = 1.0

    # m == 1
    am[mask1] = np.arcsin(np.tanh(u[mask1]))
    sn[mask1] = np.tanh(u[mask1])
    cn[mask1] = 1.0 / np.cosh(u[mask1])
    dn[mask1] = 1.0 / np.cosh(u[mask1])

    if not np.any(maskN):
        return sn, cn, dn, am

    u_g = u[maskN]
    m_g = m[maskN]
    n_g = u_g.size

    # --- AGM: build sequences a, b, c ---
    # a[0]=1, b[0]=sqrt(1-m), c[0]=sqrt(m)
    a_arr = np.ones((n_g, _AGM_ITERS + 1))
    b_arr = np.empty((n_g, _AGM_ITERS + 1))
    c_arr = np.empty((n_g, _AGM_ITERS + 1))
    b_arr[:, 0] = np.sqrt(1.0 - m_g)
    c_arr[:, 0] = np.sqrt(m_g)

    for i in range(1, _AGM_ITERS + 1):
        a_arr[:, i] = 0.5 * (a_arr[:, i - 1] + b_arr[:, i - 1])
        b_arr[:, i] = np.sqrt(a_arr[:, i - 1] * b_arr[:, i - 1])
        c_arr[:, i] = 0.5 * (a_arr[:, i - 1] - b_arr[:, i - 1])

    # per-element convergence index: first i where |c[i]| <= eps
    n_conv = np.full(n_g, _AGM_ITERS, dtype=int)
    for i in range(1, _AGM_ITERS + 1):
        just = (np.abs(c_arr[:, i]) <= eps) & (np.abs(c_arr[:, i - 1]) > eps)
        n_conv[just & (n_conv == _AGM_ITERS)] = i - 1

    # a_final[k] = a at convergence step for element k
    a_final = a_arr[np.arange(n_g), n_conv]

    # Starting amplitude for descending Landen:
    # phi_N = 2^N * a_N * u  (the "linearised" amplitude at convergence)
    phin = (2.0 ** n_conv) * a_final * u_g

    # --- Descending Landen back-substitution ---
    # phi_{i-1} = 0.5 * (asin(c_i/a_i * sin(phi_i)) + phi_i)  for i=N,...,1
    for i in range(_AGM_ITERS, 0, -1):
        active = n_conv >= i
        if not np.any(active):
            continue
        ratio = np.where(active,
                         c_arr[:, i] / a_arr[:, i],
                         np.zeros(n_g))
        arg = np.clip(ratio * np.sin(phin), -1.0, 1.0)
        phin_new = 0.5 * (np.arcsin(arg) + phin)
        phin = np.where(active, phin_new, phin)

    am[maskN] = phin
    sn[maskN] = np.sin(phin)
    cn[maskN] = np.cos(phin)
    dn[maskN] = np.sqrt(np.clip(1.0 - m_g * np.sin(phin) ** 2, 0.0, None))

    return sn, cn, dn, am
