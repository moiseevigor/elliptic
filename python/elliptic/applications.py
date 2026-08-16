"""Application-level helpers built on top of the core elliptic functions."""
from __future__ import annotations
import numpy as np

from ._xputils import get_xp
from .elliptic12 import elliptic12


def arclength_ellipse(a, b, theta0=0.0, theta1=None):
    """Arc length of an ellipse from angle theta0 to theta1.

    The ellipse is parameterised as  x = a cos t,  y = b sin t.
    Angle t is measured from the positive a-axis (semi-major or semi-minor).

    Parameters
    ----------
    a : float
        First semi-axis.
    b : float
        Second semi-axis.
    theta0 : float, optional
        Start angle in radians. Default 0.
    theta1 : float, optional
        End angle in radians. Default 2π (full perimeter).

    Returns
    -------
    arc : float
        Arc length.

    Examples
    --------
    Full perimeter of ellipse with a=5, b=10 (matches Mathematica):

    >>> arclength_ellipse(5, 10)  # doctest: +ELLIPSIS
    48.4422...

    Notes
    -----
    For a circle (a == b) the formula reduces to a*(theta1 - theta0).
    When b > a the standard formula is used:  b * E(theta1|1-(a/b)^2).
    When a > b the complement formula applies: a * E(π/2-theta|1-(b/a)^2).
    """
    if theta1 is None:
        theta1 = 2.0 * np.pi

    xp = get_xp(a, b, theta0, theta1)
    a = xp.asarray(a, dtype=xp.float64)
    b = xp.asarray(b, dtype=xp.float64)
    theta0 = xp.asarray(theta0, dtype=xp.float64)
    theta1 = xp.asarray(theta1, dtype=xp.float64)
    a, b, theta0, theta1 = xp.broadcast_arrays(a, b, theta0, theta1)

    # Give ordinary NumPy callers an explicit domain error.  Traced backends
    # cannot branch on array values, so invalid elements are marked NaN below.
    if xp is np and (np.any(a <= 0.0) or np.any(b <= 0.0)):
        raise ValueError("ellipse semi-axes must be strictly positive")

    valid = (a > 0.0) & (b > 0.0)
    a_safe = xp.where(valid, a, xp.ones_like(a))
    b_safe = xp.where(valid, b, xp.ones_like(b))

    # Evaluate both orientations elementwise.  The old scalar ``float`` casts
    # rejected arrays and JAX tracers even though the rest of the public API is
    # backend-native.
    m_b = 1.0 - (a_safe / b_safe) ** 2
    m_a = 1.0 - (b_safe / a_safe) ** 2

    _, E1_b, _ = elliptic12(theta1, xp.where(b > a, m_b, xp.zeros_like(m_b)))
    _, E0_b, _ = elliptic12(theta0, xp.where(b > a, m_b, xp.zeros_like(m_b)))

    comp1 = np.pi / 2.0 - theta1
    comp0 = np.pi / 2.0 - theta0
    _, E1_a, _ = elliptic12(comp1, xp.where(a > b, m_a, xp.zeros_like(m_a)))
    _, E0_a, _ = elliptic12(comp0, xp.where(a > b, m_a, xp.zeros_like(m_a)))

    arc_b = b_safe * (E1_b - E0_b)
    arc_a = a_safe * (E0_a - E1_a)
    arc_circle = a_safe * xp.abs(theta1 - theta0)
    arc = xp.where(b > a, arc_b, xp.where(a > b, arc_a, arc_circle))
    return xp.where(valid, arc, xp.full_like(arc, np.nan))
