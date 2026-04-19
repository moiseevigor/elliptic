"""AGM kernel shared by elliptic12 and ellipj.

All loops run a *fixed* 25 iterations so the code is JAX-traceable
(no dynamic stopping condition). For float64 the AGM converges in at
most ~25 halvings, so the extra no-op iterations cost nothing.
"""
from __future__ import annotations

_AGM_ITERS = 25   # safe upper bound for float64


def agm_coeffs(m, xp):
    """Compute AGM sequences a, b, c of shape (N, _AGM_ITERS+1).

    Returns
    -------
    a, b, c : (N, iters+1) arrays
    n       : (N,) int array — converged iteration index per element
    """
    N = m.shape[0]
    iters = _AGM_ITERS

    # Allocate: we'll build columns one by one (works for numpy/torch/jax)
    a = [xp.ones(N, dtype=xp.float64)]
    b = [xp.sqrt(1.0 - m)]
    c = [xp.sqrt(m)]

    for _ in range(iters):
        a_new = 0.5 * (a[-1] + b[-1])
        b_new = xp.sqrt(a[-1] * b[-1])
        c_new = 0.5 * (a[-1] - b[-1])
        a.append(a_new)
        b.append(b_new)
        c.append(c_new)

    return a, b, c


def agm_n(c, xp):
    """Return per-element convergence index n from the c sequence."""
    import numpy as _np
    # n[k] = first i where |c[i]| ≈ 0 (i.e., c[i+1] < eps)
    eps = float(xp.finfo(xp.float64).eps) if hasattr(xp, 'finfo') else 2.2e-16
    iters = len(c) - 1
    # Build as numpy for indexing simplicity; convert back if needed
    try:
        c_np = _np.stack([_np.asarray(ci) for ci in c], axis=0)  # (iters+1, N)
    except Exception:
        c_np = _np.array([_np.array(ci) for ci in c])
    N = c_np.shape[1]
    n = _np.zeros(N, dtype=_np.intp)
    for i in range(1, iters + 1):
        just_converged = (abs(c_np[i]) <= eps) & (abs(c_np[i - 1]) > eps)
        n[just_converged & (n == 0)] = i - 1
    # anything still 0 gets maximum
    n[n == 0] = iters - 1
    return n
