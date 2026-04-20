"""Nome q(m) and its inverse m(q)."""
from __future__ import annotations
import numpy as np
from array_api_compat import array_namespace


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
    from scipy.special import ellipk
    m = np.asarray(m, dtype=np.float64)
    xp = array_namespace(m)
    m_np = np.asarray(m).ravel()
    K  = ellipk(m_np)
    Kp = ellipk(1.0 - m_np)
    q  = np.exp(-np.pi * Kp / K)
    q  = q.reshape(np.asarray(m).shape)
    return xp.asarray(q)


def inversenomeq(q):
    """Inverse nome: parameter m from nome q.

    Uses ``scipy.optimize.brentq`` to invert ``nomeq``. In double precision the
    representable range is roughly q ∈ [0, 0.779]; beyond this, m(q) exceeds
    1 - 2⁻⁵³ and cannot be represented.

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
    from scipy.optimize import brentq
    from scipy.special import ellipk

    q = np.asarray(q, dtype=np.float64)
    xp = array_namespace(q)
    q_flat = np.asarray(q).ravel()

    if np.any(q_flat < 0) or np.any(q_flat >= 1):
        raise ValueError("q must be in [0, 1)")

    m_hi  = np.nextafter(1.0, 0.0)        # largest f64 strictly < 1
    q_max = float(np.exp(-np.pi * ellipk(1.0 - m_hi) / ellipk(m_hi)))
    if np.any(q_flat >= q_max):
        raise ValueError(
            f"inversenomeq: q must be < {q_max:.15f} in double precision "
            "(the essential singularity of m(q) at q=1 cannot be resolved in f64)"
        )
    if np.any(q_flat > 0.76):
        warnings.warn("inversenomeq: accuracy degrades for q > 0.76 (near m=1 singularity)",
                      RuntimeWarning, stacklevel=2)

    def _nomeq_scalar(m_val):
        K  = float(ellipk(m_val))
        Kp = float(ellipk(1.0 - m_val))
        return float(np.exp(-np.pi * Kp / K))

    m_out = np.empty_like(q_flat)
    for i, qi in enumerate(q_flat):
        if qi == 0.0:
            m_out[i] = 0.0
        else:
            m_out[i] = brentq(lambda m: _nomeq_scalar(m) - qi, 0.0, m_hi, xtol=1e-14)

    m_out = m_out.reshape(np.asarray(q).shape)
    return xp.asarray(m_out)
