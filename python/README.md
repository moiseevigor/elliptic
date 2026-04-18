# elliptic (Python)

Standalone Python implementation of elliptic integrals and functions.
Works identically on NumPy arrays, PyTorch tensors (CPU and CUDA), and JAX arrays — no scipy runtime dependency.

## Install

```bash
pip install elliptic                   # NumPy only
pip install "elliptic[torch]"          # + PyTorch backend
pip install "elliptic[jax]"            # + JAX backend
pip install "elliptic[dev]"            # + test deps (scipy, mpmath, pytest)
```

## Quick start

```python
import numpy as np
import elliptic

# F(phi, m), E(phi, m), Jacobi Zeta Z(phi, m)
phi = np.linspace(0, np.pi/2, 200)
F, E, Z = elliptic.elliptic12(phi, 0.7)

# Jacobi elliptic functions
sn, cn, dn, am = elliptic.ellipj(np.linspace(0, 3, 100), 0.7)

# Third kind Pi(phi, m, n)  — missing from scipy (issue #4452)
Pi = elliptic.elliptic3(phi, 0.7, 0.3)

# Associate integrals B, D, J  (DLMF §19.2)
B, D, J = elliptic.ellipticBDJ(phi, 0.7, 0.3)
# F = B+D,  E = B+(1-m)*D,  Pi = B+D+n*J

# Complete associate integrals
B_c, D_c, S_c = elliptic.ellipticBD(np.array([0.3, 0.5, 0.7]))

# Carlson symmetric forms  (DLMF §19.16)
RF = elliptic.carlsonRF(0.0, 0.3, 1.0)   # = K(0.7)
RD = elliptic.carlsonRD(c2, d2, 1.0)

# Bulirsch generalised complete integrals
C  = elliptic.cel(kc, p, a, b)
K  = elliptic.cel1(kc)                    # = K(1 - kc^2)
E2 = elliptic.cel2(kc, 1.0, 1 - kc**2)  # = E(1 - kc^2)
Pi = elliptic.cel3(kc, 1 - n)            # = Pi(n | 1 - kc^2)

# Weierstrass functions
P  = elliptic.weierstrassP(z, e1, e2, e3)
Zf = elliptic.weierstrassZeta(z, e1, e2, e3)
S  = elliptic.weierstrassSigma(z, e1, e2, e3)
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
| Jacobi sn/cn/dn | ✓ (real only) | ✓ |
| **Π(φ,m,n) — 3rd kind** | ✗ (open since 2015) | ✓ |
| Carlson RF/RD/RJ/RC | ✓ since v1.8 | ✓ |
| **B, D, J associate** | ✗ | ✓ |
| **Bulirsch cel** | ✗ | ✓ |
| **Weierstrass P/ζ/σ** | ✗ | ✓ |
| **PyTorch GPU** | ✗ | ✓ |
| **JAX jit/vmap** | ✗ | ✓ |

All algorithms run identically across backends — the same Python code dispatches to NumPy, PyTorch CUDA, or JAX TPU.

## References

- Abramowitz & Stegun, *Handbook of Mathematical Functions*, §17–18
- NIST DLMF §19, §22, §23 — https://dlmf.nist.gov
- Fukushima (2015), "Elliptic functions and elliptic integrals for celestial mechanics"
- Bulirsch (1965), "Numerical computation of elliptic integrals"
- Carlson (1995), "Numerical Computation of Real or Complex Elliptic Integrals"
