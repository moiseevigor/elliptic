"""Shared utilities: backend detection, type coercion, parallel chunking."""
from __future__ import annotations

import os
import math
from concurrent.futures import ProcessPoolExecutor
from typing import Any

from array_api_compat import array_namespace, is_numpy_array
import numpy as _np


def get_xp(*args: Any):
    """Return the array namespace for the inputs (numpy, torch, jax.numpy, …).

    Scalars and Python lists are treated as numpy by default.
    """
    # Convert any non-array args to numpy so array_namespace can detect them
    converted = []
    for a in args:
        if isinstance(a, _np.ndarray):
            converted.append(a)
        else:
            try:
                # torch.Tensor, jax.Array, etc.
                _ = a.shape
                converted.append(a)
            except AttributeError:
                converted.append(_np.asarray(a))
    if not converted:
        return _np
    try:
        return array_namespace(*converted)
    except TypeError:
        return _np


def to_float64(x: Any, xp):
    """Cast *x* to float64 in namespace *xp*, preserving device."""
    return xp.asarray(x, dtype=xp.float64)


def broadcast_arrays(*args, xp):
    """Broadcast all arrays to a common shape using *xp*."""
    import numpy as np  # used only for shape computation
    shapes = [xp.asarray(a).shape for a in args]
    try:
        out_shape = np.broadcast_shapes(*shapes)
    except AttributeError:  # numpy < 1.20
        out_shape = np.broadcast(*[np.empty(s) for s in shapes]).shape
    return tuple(xp.broadcast_to(xp.asarray(a, dtype=xp.float64), out_shape) for a in args)


# ---------------------------------------------------------------------------
# Parallel chunking for NumPy path
# ---------------------------------------------------------------------------
_PARALLEL_THRESHOLD = 50_000


def _n_workers() -> int:
    return os.cpu_count() or 1


def parallel_apply(fn, *arrays, threshold: int = _PARALLEL_THRESHOLD):
    """Run *fn* across chunks if the input is large enough to benefit.

    Falls back to serial when the array is small or only 1 CPU is available.
    *fn* must accept flat numpy arrays and return a tuple of flat arrays.
    """
    import numpy as np
    N = arrays[0].size
    n = _n_workers()
    if n <= 1 or N < threshold:
        return fn(*arrays)
    orig_shape = arrays[0].shape
    flat = [a.ravel() for a in arrays]
    chunk = math.ceil(N / n)
    slices = [slice(i * chunk, min((i + 1) * chunk, N)) for i in range(n)]
    chunks = [[f[s] for f in flat] for s in slices]
    with ProcessPoolExecutor(max_workers=n) as pool:
        results = list(pool.map(lambda c: fn(*c), chunks))
    # results is a list of tuples (or single arrays); re-assemble
    if isinstance(results[0], tuple):
        return tuple(np.concatenate([r[i] for r in results]).reshape(orig_shape)
                     for i in range(len(results[0])))
    return np.concatenate(results).reshape(orig_shape)
