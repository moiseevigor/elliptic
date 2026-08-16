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

    Uses a fixed-iteration bisection against the library's own Carlson-based
    ``nomeq`` implementation. In double precision the representable range is
    roughly q ∈ [0, 0.779]; beyond this, m(q) exceeds 1 - 2⁻⁵³ and cannot be
    represented.

    Parameters
    ----------
    q : array_like
        Nome values in [0, q_max) with q_max ≈ 0.7789534...

    Returns
    -------
    m : array
        Parameter m = m(q) in [0, 1).
    """
    import warnings

    xp = get_xp(q)
    q = xp.asarray(q, dtype=xp.float64)

    m_hi_scalar = np.nextafter(1.0, 0.0)
    q_max = float(_q_from_m_xp(np, np.asarray(m_hi_scalar)))

    if xp is np:
        if np.any((q < 0.0) | (q >= 1.0)):
            raise ValueError("q must be in [0, 1)")
        if np.any(q >= q_max):
            raise ValueError(
                f"inversenomeq: q must be < {q_max:.15f} in double precision "
                "(the essential singularity of m(q) at q=1 cannot be resolved in f64)"
            )
        if np.any(q > 0.76):
            warnings.warn(
                "inversenomeq: accuracy degrades for q > 0.76 (near m=1 singularity)",
                RuntimeWarning,
                stacklevel=2,
            )

    valid = (q >= 0.0) & (q < q_max)
    q_safe = xp.where(valid, q, xp.zeros_like(q))
    lo = xp.zeros_like(q_safe)
    hi = xp.full_like(q_safe, m_hi_scalar)
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        q_mid = _q_from_m_xp(xp, mid)
        lower = q_mid < q_safe
        lo = xp.where(lower, mid, lo)
        hi = xp.where(lower, hi, mid)
    result = 0.5 * (lo + hi)
    result = xp.where(q == 0.0, xp.zeros_like(result), result)
    return xp.where(valid, result, xp.full_like(result, math.nan))
