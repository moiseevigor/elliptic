"""elliptic — Elliptic integrals and functions for NumPy, PyTorch, and JAX.

Standalone implementation (no scipy runtime dependency).
Pass a torch.Tensor or jax.Array for automatic GPU/TPU dispatch.

Public API
----------
Incomplete integrals:
    elliptic12(u, m)           -> F, E, Z
    elliptic3(u, m, n)         -> Pi

Jacobi elliptic functions:
    ellipj(u, m)               -> sn, cn, dn, am

Associate integrals (incomplete):
    ellipticBDJ(phi, m[, n])   -> B, D, J

Complete associate integrals:
    ellipticBD(m)              -> B, D, S

Jacobi-argument forms:
    jacobiEDJ(u, m[, n])       -> Eu, Du, Ju

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

Weierstrass functions:
    weierstrassP(z, e1, e2, e3)
    weierstrassZeta(z, e1, e2, e3)
    weierstrassSigma(z, e1, e2, e3)
"""

from .elliptic12    import elliptic12
from .elliptic3     import elliptic3
from .ellipj        import ellipj
from .ellipticBDJ   import ellipticBDJ
from .ellipticBD    import ellipticBD
from .jacobi_edj    import jacobiEDJ
from .carlson       import carlsonRF, carlsonRD, carlsonRJ, carlsonRC
from .bulirsch      import cel, cel1, cel2, cel3
from .weierstrass   import weierstrassP, weierstrassZeta, weierstrassSigma

__all__ = [
    "elliptic12", "elliptic3", "ellipj",
    "ellipticBDJ", "ellipticBD", "jacobiEDJ",
    "carlsonRF", "carlsonRD", "carlsonRJ", "carlsonRC",
    "cel", "cel1", "cel2", "cel3",
    "weierstrassP", "weierstrassZeta", "weierstrassSigma",
]
