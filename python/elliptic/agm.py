"""Arithmetic-geometric mean."""
from __future__ import annotations
import numpy as np
from array_api_compat import array_namespace


def agm(a0, b0):
    """Arithmetic-geometric mean of *a0* and *b0*.

    Iterates  a_{n+1} = (a_n + b_n)/2,  b_{n+1} = sqrt(a_n * b_n)
    until convergence (25 fixed steps — safe for float64).

    Parameters
    ----------
    a0, b0 : array_like
        Non-negative inputs of the same shape (or scalar).

    Returns
    -------
    agm_val : array
        AGM value with shape broadcast(a0, b0).
    """
    a0 = np.asarray(a0, dtype=np.float64)
    b0 = np.asarray(b0, dtype=np.float64)
    xp = array_namespace(a0, b0)
    a0, b0 = xp.broadcast_arrays(a0, b0)
    a = np.asarray(a0, dtype=np.float64).copy()
    b = np.asarray(b0, dtype=np.float64).copy()
    for _ in range(25):
        a, b = 0.5 * (a + b), np.sqrt(a * b)
    return xp.asarray(a)
