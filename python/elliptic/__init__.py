"""elliptic — Elliptic integrals and functions for NumPy, PyTorch, and JAX.

Standalone implementation (no scipy runtime dependency).
Pass a torch.Tensor or jax.Array for automatic GPU/TPU dispatch.

Public API
----------
Incomplete integrals:
    elliptic12(u, m)             -> F, E, Z
    elliptic3(u, m, n)           -> Pi
    elliptic12i(u, m)            -> Fi, Ei, Zi  (complex u)

Jacobi elliptic functions:
    ellipj(u, m)                 -> sn, cn, dn, am
    ellipji(u, m)                -> sn, cn, dn  (complex u)

Associate integrals (incomplete):
    ellipticBDJ(phi, m[, n])     -> B, D, J

Complete associate integrals:
    ellipticBD(m)                -> B, D, S

Jacobi-argument forms:
    jacobiEDJ(u, m[, n])         -> Eu, Du, Ju

Carlson symmetric forms:
    carlsonRF(x, y, z)
    carlsonRD(x, y, z)
    carlsonRJ(x, y, z, p)
    carlsonRC(x, y)

Bulirsch generalised complete integrals:
    cel(kc, p, a, b)
    cel1(kc)
    cel2(kc, a, b)
    cel3(kc, p)

Jacobi theta functions:
    jacobiThetaEta(u, m)         -> Th, H
    theta(j, v, m)               -> Th_j
    theta_prime(j, v, m)         -> Th_j, dTh_j/dv

Weierstrass functions:
    weierstrassP(z, e1, e2, e3)
    weierstrassZeta(z, e1, e2, e3)
    weierstrassSigma(z, e1, e2, e3)
    weierstrassPPrime(z, e1, e2, e3)
    weierstrassInvariants(e1, e2, e3)  -> g2, g3, Delta

Nome and inverse:
    nomeq(m)                     -> q
    inversenomeq(q)              -> m

Inverse integral:
    inverselliptic2(E, m)        -> phi

AGM:
    agm(a, b)                    -> agm(a, b)

Applications:
    arclength_ellipse(a, b[, theta0, theta1])
"""

from .elliptic12         import elliptic12
from .elliptic3          import elliptic3
from .ellipj             import ellipj
from .ellipticBDJ        import ellipticBDJ
from .ellipticBD         import ellipticBD
from .jacobi_edj         import jacobiEDJ
from .carlson            import carlsonRF, carlsonRD, carlsonRJ, carlsonRC
from .bulirsch           import cel, cel1, cel2, cel3
from .weierstrass        import (weierstrassP, weierstrassZeta, weierstrassSigma,
                                  weierstrassPPrime, weierstrassInvariants)
from .theta              import jacobiThetaEta, theta, theta_prime
from .complex_elliptic   import elliptic12i, ellipji
from .nome               import nomeq, inversenomeq
from .inverse            import inverselliptic2
from .agm                import agm
from .applications       import arclength_ellipse

__all__ = [
    # integrals
    "elliptic12", "elliptic3", "elliptic12i",
    # Jacobi
    "ellipj", "ellipji",
    # associate
    "ellipticBDJ", "ellipticBD", "jacobiEDJ",
    # Carlson
    "carlsonRF", "carlsonRD", "carlsonRJ", "carlsonRC",
    # Bulirsch
    "cel", "cel1", "cel2", "cel3",
    # theta
    "jacobiThetaEta", "theta", "theta_prime",
    # Weierstrass
    "weierstrassP", "weierstrassZeta", "weierstrassSigma",
    "weierstrassPPrime", "weierstrassInvariants",
    # nome / inverse
    "nomeq", "inversenomeq", "inverselliptic2",
    # misc
    "agm", "arclength_ellipse",
]
