"""Nome q(m) and its inverse m(q)."""
from __future__ import annotations
import math
import numpy as np
from ._xputils import get_xp
from .theta import _q_from_m_xp


def nomeq(m):
    """Nome q = exp(-π K'(m) / K(m)).

    Parameters
    ----------
    m : array_like
        Parameter, 0 <= m <= 1.

    Returns
    -------
    q : array
        Nome in [0, 1).
    """
    xp = get_xp(m)
    m = xp.asarray(m, dtype=xp.float64)
    q = _q_from_m_xp(xp, m)
    return xp.where(m == 1.0, xp.ones_like(q), q)


def inversenomeq(q):
    """Inverse nome: parameter m from nome q.

    Closed form m = (theta2(0,q)/theta3(0,q))^4 (DLMF 20.9.1), exact at every
    scale. In double precision the representable range is q ∈ [0, q_max] with
    q_max = nomeq(nextafter(1, 0)) ≈ 0.7789; beyond this m(q) exceeds 1 - 2⁻⁵³
    and cannot be represented.

    Parameters
    ----------
    q : array_like
        Nome values in [0, q_max) with q_max ≈ 0.7789534...

    Returns
    -------
    m : array
        Parameter m = m(q) in [0, 1).
    """
    xp = get_xp(q)
    q = xp.asarray(q, dtype=xp.float64)

    m_hi_scalar = np.nextafter(1.0, 0.0)
    q_max = float(_q_from_m_xp(np, np.asarray(m_hi_scalar)))

    if xp is np:
        if np.any((q < 0.0) | (q >= 1.0)):
            raise ValueError("q must be in [0, 1)")
        if np.any(q > q_max):
            raise ValueError(
                f"inversenomeq: q must be <= {q_max:.15f} in double precision "
                "(the essential singularity of m(q) at q=1 cannot be resolved in f64)"
            )

    # Closed form, DLMF 20.9.1:  m = (theta2(0,q) / theta3(0,q))^4.
    # Exact at every scale -- the previous 64-step bisection in m had an
    # absolute resolution floor of 2^-64, so m(1e-30) came back 2.7e-20
    # instead of 1.6e-29 (nine orders of magnitude off).
    #   theta2(0,q) = 2 q^(1/4) sum q^(n(n+1)),  theta3(0,q) = 1 + 2 sum q^(n^2)
    # The q^(1/4) factor is kept outside the ratio so tiny q cannot underflow.
    valid = (q >= 0.0) & (q <= q_max)
    q_safe = xp.where(valid, q, xp.zeros_like(q))
    s2 = xp.ones_like(q_safe)                   # sum q^(n(n+1)), n >= 0
    s3 = xp.ones_like(q_safe)                   # theta3 = 1 + 2 sum q^(n^2)
    for n in range(1, 31):
        s2 = s2 + q_safe ** (n * (n + 1))
        s3 = s3 + 2.0 * q_safe ** (n * n)
    result = 16.0 * q_safe * (s2 / s3) ** 4
    return xp.where(valid, result, xp.full_like(result, math.nan))
