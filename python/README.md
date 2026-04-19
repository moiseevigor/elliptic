
# elliptic (Python)

Standalone Python implementation of elliptic integrals and functions.
Works identically on NumPy arrays, PyTorch tensors (CPU and CUDA), and JAX arrays — no scipy runtime dependency.

[![tests](https://github.com/moiseevigor/elliptic/actions/workflows/python.yml/badge.svg)](https://github.com/moiseevigor/elliptic/actions/workflows/python.yml)

## Install

```bash
pip install elliptic-functions                   # NumPy only
pip install "elliptic-functions[torch]"          # + PyTorch backend
pip install "elliptic-functions[jax]"            # + JAX backend
pip install "elliptic-functions[dev]"            # + test deps (scipy, mpmath, pytest)
```

## Quick start

```python
import numpy as np
import elliptic

phi = np.linspace(0, np.pi/2, 200)
m   = 0.7

# Incomplete integrals F(φ|m), E(φ|m), Jacobi zeta Z(φ|m)
F, E, Z = elliptic.elliptic12(phi, m)

# Third kind Π(φ,m,n)  — missing from scipy (issue #4452)
Pi = elliptic.elliptic3(phi, m, 0.3)

# Jacobi elliptic functions sn, cn, dn, am
sn, cn, dn, am = elliptic.ellipj(np.linspace(0, 3, 100), m)

# Associate incomplete integrals B, D, J  (DLMF §19.2)
B, D, J = elliptic.ellipticBDJ(phi, m, 0.3)
# identities: F = B+D,  E = B+(1-m)D,  Pi = B+D+n*J

# Complete associate integrals B(m), D(m), S(m)
B_c, D_c, S_c = elliptic.ellipticBD(np.array([0.3, 0.5, 0.7]))

# Carlson symmetric forms  (DLMF §19.16)
RF = elliptic.carlsonRF(0.0, 1.0 - m, 1.0)    # = K(m)
RD = elliptic.carlsonRD(0.0, 1.0 - m, 1.0)    # = 3D(m)

# Bulirsch generalised complete integrals
K  = elliptic.cel1(np.sqrt(1 - m))             # = K(m)
E2 = elliptic.cel2(np.sqrt(1 - m), 1.0, 1-m)  # = E(m)

# Jacobi theta functions (4 types)
Th3 = elliptic.theta(3, phi, m)
Th,  H  = elliptic.jacobiThetaEta(np.linspace(0, 1, 50), m)
Th3, dTh3 = elliptic.theta_prime(3, phi, m)

# Weierstrass functions
e1, e2, e3 = 2.0, 0.5, -2.5
P   = elliptic.weierstrassP(0.4, e1, e2, e3)
Zf  = elliptic.weierstrassZeta(0.4, e1, e2, e3)
S   = elliptic.weierstrassSigma(0.4, e1, e2, e3)
dP  = elliptic.weierstrassPPrime(0.4, e1, e2, e3)
g2, g3, Delta = elliptic.weierstrassInvariants(e1, e2, e3)

# Nome and inverse
q = elliptic.nomeq(m)
m_back = elliptic.inversenomeq(q)

# Inverse E(phi|m): given E value, solve for phi
phi_inv = elliptic.inverselliptic2(E, m)

# Arithmetic-geometric mean
agm_val = elliptic.agm(1.0, np.sqrt(1 - m))   # = π/(2K(m))

# Complex arguments
import numpy as np
u_c = np.array([0.4 + 0.3j, 0.8 - 0.2j])
Fi, Ei, Zi = elliptic.elliptic12i(u_c, m)
sni, cni, dni = elliptic.ellipji(u_c, m)

# Application: ellipse arc length
arc = elliptic.arclength_ellipse(5, 10)       # full perimeter ≈ 48.44
arc = elliptic.arclength_ellipse(5, 10, 0, np.pi/4)  # quarter arc
```

## GPU (PyTorch)

Pass a `torch.Tensor` on any device:

```python
import torch, elliptic

phi = torch.linspace(0.01, 1.5, 1_000_000, dtype=torch.float64, device="cuda")
F, E, Z = elliptic.elliptic12(phi, torch.full_like(phi, 0.7))
```

## JAX (jit / vmap)

```python
import jax, jax.numpy as jnp, elliptic

phi = jnp.linspace(0.01, 1.5, 10_000)
F, E, Z = jax.jit(elliptic.elliptic12)(phi, 0.7)
```

## Comparison with scipy.special

| Function | scipy | elliptic |
|---|:---:|:---:|
| K(m), E(m) | ✓ | ✓ |
| F(φ,m), E(φ,m) | ✓ | ✓ |
| Jacobi sn/cn/dn (real) | ✓ | ✓ |
| Jacobi sn/cn/dn (complex) | ✗ | ✓ |
| F(u,m), E(u,m) (complex) | ✗ | ✓ |
| **Π(φ,m,n) — 3rd kind** | ✗ (open since 2015) | ✓ |
| Carlson RF/RD/RJ/RC | ✓ since v1.8 | ✓ |
| **B, D, J associate** | ✗ | ✓ |
| **Bulirsch cel** | ✗ | ✓ |
| **Jacobi theta θ₁…θ₄** | ✗ | ✓ |
| **Weierstrass P/ζ/σ/℘′** | ✗ | ✓ |
| **Weierstrass invariants g₂,g₃** | ✗ | ✓ |
| **Nome q(m) and m(q)** | ✗ | ✓ |
| **Inverse E(φ,m)** | ✗ | ✓ |
| **AGM** | ✗ | ✓ |
| **Ellipse arc length** | ✗ | ✓ |
| **PyTorch GPU** | ✗ | ✓ |
| **JAX jit/vmap** | ✗ | ✓ |

All algorithms run identically across backends — the same Python code dispatches to NumPy, PyTorch CUDA, or JAX TPU.

## References

- Abramowitz & Stegun, *Handbook of Mathematical Functions*, §16–18
- NIST DLMF §19, §22, §23 — https://dlmf.nist.gov
- Fukushima (2015), "Elliptic functions and elliptic integrals for celestial mechanics"
- Bulirsch (1965), "Numerical computation of elliptic integrals"
- Carlson (1995), "Numerical Computation of Real or Complex Elliptic Integrals"
- Boyd (2012), "Numerical inversion of the incomplete elliptic integral of the second kind"
