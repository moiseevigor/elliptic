"""Arithmetic-geometric mean — native on any array backend."""
from __future__ import annotations
from ._xputils import get_xp


def agm(a0, b0):
    """Arithmetic-geometric mean of *a0* and *b0*.

    Iterates  a_{n+1} = (a_n + b_n)/2,  b_{n+1} = sqrt(a_n * b_n)
    for 25 fixed steps (safe for float64).  Runs natively on NumPy,
    PyTorch CUDA, and JAX.
    """
    xp = get_xp(a0, b0)
    a = xp.asarray(a0, dtype=xp.float64)
    b = xp.asarray(b0, dtype=xp.float64)
    a, b = xp.broadcast_arrays(a, b)
    for _ in range(25):
        a, b = 0.5 * (a + b), xp.sqrt(a * b)
    return a
