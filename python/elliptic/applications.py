"""Application-level helpers built on top of the core elliptic functions."""
from __future__ import annotations
import numpy as np
from array_api_compat import array_namespace
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

    a = float(a); b = float(b)
    theta0 = float(theta0); theta1 = float(theta1)

    if a == b:
        return a * abs(theta1 - theta0)

    if b > a:
        m = 1.0 - (a / b) ** 2
        _, E1, _ = elliptic12(np.asarray(theta1), np.asarray(m))
        _, E0, _ = elliptic12(np.asarray(theta0), np.asarray(m))
        return float(b * (float(E1) - float(E0)))
    else:  # a > b
        m = 1.0 - (b / a) ** 2
        _, E1, _ = elliptic12(np.asarray(np.pi / 2.0 - theta1), np.asarray(m))
        _, E0, _ = elliptic12(np.asarray(np.pi / 2.0 - theta0), np.asarray(m))
        return float(a * (float(E0) - float(E1)))
