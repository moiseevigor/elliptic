"""Array-namespace detection helper shared across all modules."""
from __future__ import annotations

from array_api_compat import array_namespace, is_array_api_obj
import numpy as np


def is_numpy(xp):
    """True for the eager numpy namespace, whichever module object represents
    it: ``numpy`` itself, or ``array_api_compat.numpy`` (what
    array_namespace() returns for numpy arrays and numpy scalars).  Every
    ``xp is np`` test in the package used to miss the latter, so eager domain
    checks were silently skipped for ndarray inputs."""
    return xp is np or getattr(xp, "__name__", "") in ("numpy", "array_api_compat.numpy")


def check_range(xp, x, lo, hi, what):
    """Domain check that is honest on every backend.

    numpy (eager): raise ValueError for any value outside [lo, hi] -- NaN is
    let through and propagates.  Traced/device backends (JAX, torch): return
    a validity mask so callers can emit NaN instead of a silent placeholder
    value (ellipj(0.3, 1.5) used to return sn(0.3 | 0.5)).
    """
    valid = ~((x < lo) | (x > hi))
    if is_numpy(xp):
        bad = ~valid & ~np.isnan(x)
        if np.any(bad):
            raise ValueError(f"{what} must be in [{lo}, {hi}]")
    return valid


def get_xp(*args):
    """Return the array namespace for *args*, defaulting to numpy for plain scalars."""
    api_objs = [a for a in args if is_array_api_obj(a)]
    if api_objs:
        return array_namespace(*api_objs)
    return np
