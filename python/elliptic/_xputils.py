"""Array-namespace detection helper shared across all modules."""
from __future__ import annotations

from array_api_compat import array_namespace, is_array_api_obj
import numpy as np


# Three-term split of pi: _PI_A and _PI_B carry 25 significant bits each, so
# k*_PI_A and k*_PI_B are exact in double for |k| < 2**28 (|u| < 8e8); the
# remaining k*_PI_C rounds at eps*|k|*1.6e-8, far below eps*|r|.  Using
# float(pi) as the leading term ("Cody-Waite" with a 53-bit head) does not
# help: k*float(pi) already rounds by eps*|u| (2.3e-10 at u = 1e6), which
# Jacobi Zeta and E inherit.  Residual pi - A - B - C = 1.3e-24.
_PI_A = 3.1415926218032837        # 0x1.921fb5p+1
_PI_B = 1.5893254712295857e-08    # 0x1.110b46p-26
_PI_C = 1.5893254834760535e-08


def sub_kpi(u, k):
    """u - k*pi for integer-valued k, accurate to eps*|result| (any backend)."""
    return ((u - k * _PI_A) - k * _PI_B) - k * _PI_C


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
