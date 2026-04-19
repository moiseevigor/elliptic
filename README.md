
# elliptic — Elliptic integrals and functions

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/moiseevigor/elliptic/tree/master.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/moiseevigor/elliptic/tree/master) [![DOI](https://zenodo.org/badge/5762/moiseevigor/elliptic.svg)](https://zenodo.org/badge/latestdoi/5762/moiseevigor/elliptic) [![tests](https://github.com/moiseevigor/elliptic/actions/workflows/python.yml/badge.svg)](https://github.com/moiseevigor/elliptic/actions/workflows/python.yml)

Full-precision implementations of elliptic integrals and functions for **MATLAB/Octave** and **Python**.
No external toolboxes or heavy dependencies — works standalone.

---

## Choose your language

| | MATLAB / Octave | Python |
|---|---|---|
| **Source** | [`matlab/`](matlab/) | [`python/`](python/) |
| **README** | [`matlab/README.md`](matlab/README.md) | [`python/README.md`](python/README.md) |
| **Install** | `git clone` + `setup` | `pip install elliptic-functions` |
| **Backends** | CPU · parfor · gpuArray | NumPy · PyTorch · JAX |
| **GPU** | ✓ CUDA via gpuArray | ✓ CUDA via PyTorch / JAX |

---

## MATLAB / Octave

```matlab
git clone https://github.com/moiseevigor/elliptic.git
cd elliptic/matlab
setup    % adds matlab/src to the MATLAB/Octave path
```

```matlab
phi = linspace(0, pi/2, 200);  m = 0.7;

[F, E, Z] = elliptic12(phi, m);           % F(φ|m), E(φ|m), Jacobi zeta
Pi        = elliptic3(phi, m, 0.3);       % Π(φ,n|m)
[sn,cn,dn,am] = ellipj(linspace(0,3,100), m);
[B, D, J] = ellipticBDJ(phi, m, 0.3);    % associate integrals
arc       = arclength_ellipse(5, 10);     % ≈ 48.44
```

Full API, GPU and parallel examples → **[`matlab/README.md`](matlab/README.md)**

---

## Python

```bash
pip install elliptic-functions               # NumPy only
pip install "elliptic-functions[torch]"      # + PyTorch (CPU and CUDA)
pip install "elliptic-functions[jax]"        # + JAX (jit / vmap / TPU)
```

```python
import numpy as np, elliptic

phi = np.linspace(0, np.pi/2, 200);  m = 0.7

F, E, Z = elliptic.elliptic12(phi, m)          # F(φ|m), E(φ|m), Jacobi zeta
Pi      = elliptic.elliptic3(phi, m, 0.3)      # Π(φ,n|m)
sn, cn, dn, am = elliptic.ellipj(np.linspace(0, 3, 100), m)
B, D, J = elliptic.ellipticBDJ(phi, m, 0.3)   # associate integrals
arc     = elliptic.arclength_ellipse(5, 10)    # ≈ 48.44
```

The same code runs unchanged on NumPy arrays, PyTorch tensors, or JAX arrays:

```python
import torch, elliptic

phi_gpu = torch.linspace(0, torch.pi/2, 200, device="cuda")
F, E, Z = elliptic.elliptic12(phi_gpu, 0.7)   # computed on GPU
```

Full API, GPU benchmarks, and JAX examples → **[`python/README.md`](python/README.md)**

---

## Function coverage

| Function | MATLAB / Octave | Python |
|---|:---:|:---:|
| F(φ\|m), E(φ\|m), Z(φ\|m) | ✓ | ✓ |
| Π(φ,n\|m) — third kind | ✓ | ✓ |
| Jacobi sn / cn / dn / am (real) | ✓ | ✓ |
| Jacobi sn / cn / dn (complex) | ✓ | ✓ |
| F, E, Z (complex argument) | ✓ | ✓ |
| Jacobi theta ϑ₁…ϑ₄ | ✓ | ✓ |
| Associate B, D, J integrals | ✓ | ✓ |
| Carlson RF / RD / RJ / RC | ✓ | ✓ |
| Bulirsch cel / cel1 / cel2 / cel3 | ✓ | ✓ |
| Weierstrass ℘ / ζ / σ / ℘′ | ✓ | ✓ |
| Weierstrass invariants g₂, g₃ | ✓ | ✓ |
| Nome q(m) and m(q) | ✓ | ✓ |
| Inverse E(φ,m) | ✓ | ✓ |
| AGM | ✓ | ✓ |
| Ellipse arc length | ✓ | ✓ |
| Multi-core CPU | ✓ parfor | ✓ ProcessPoolExecutor |
| GPU acceleration | ✓ gpuArray | ✓ PyTorch CUDA / JAX |
| JAX jit / vmap | — | ✓ |

---

## Use cases

Elliptic functions arise across physics, astronomy, geometry, and signal processing wherever periodic or nonlinear dynamics appear.
The interactive examples at **https://moiseeviorg.github.io/elliptic/examples** show the functions in action.

### Celestial mechanics

| Problem | Functions |
|---|---|
| Kepler orbit arc length | `elliptic12` |
| Mercury perihelion precession | `elliptic12`, `ellipticBDJ` |
| Gravitational potential of oblate spheroid | `carlsonRF`, `carlsonRD` |

```matlab
% MATLAB: Kepler orbit arc length from periapsis to true anomaly phi
a = 1.496e11;  e = 0.0167;       % Earth: semi-major axis, eccentricity
phi = linspace(0, 2*pi, 1000);
[~, E] = elliptic12(phi, e^2);
arc_length = a * E;               % metres from periapsis
```

```python
# Python: same calculation
import elliptic, numpy as np
a, e = 1.496e11, 0.0167
phi  = np.linspace(0, 2*np.pi, 1000)
_, E, _ = elliptic.elliptic12(phi, e**2)
arc_length = a * np.asarray(E)
```

### Physical pendulum

```matlab
% MATLAB: exact period for large-amplitude pendulum
g = 9.81;  L = 1.0;  theta0 = 2.5;   % near separatrix (143°)
k = sin(theta0 / 2);
K = cel1(sqrt(1 - k^2));              % = K(k²)
T = 4 * K * sqrt(L / g);

% exact trajectory
t = linspace(0, T/2, 500);
[sn, ~, ~] = ellipj(t * sqrt(g/L), k^2);
theta = 2 * asin(k * sn);
```

```python
# Python: same calculation
import elliptic, numpy as np
g, L, theta0 = 9.81, 1.0, 2.5
k  = np.sin(theta0 / 2)
K  = float(elliptic.cel1(np.sqrt(1 - k**2)))
T  = 4 * K * np.sqrt(L / g)

t  = np.linspace(0, T / 2, 500)
sn, _, _, _ = elliptic.ellipj(t * np.sqrt(g / L), k**2)
theta = 2 * np.arcsin(k * np.asarray(sn))
```

### Rigid body rotation (Euler–Poinsot)

```python
# Python: torque-free symmetric top
import elliptic, numpy as np
I1, I2, I3 = 1.0, 1.5, 2.0   # principal moments of inertia
omega0 = 1.0
m  = I2*(I3 - I2) / (I1*(I3 - I1)) * omega0**2
t  = np.linspace(0, 10, 1000)
sn, cn, dn, _ = elliptic.ellipj(t, m)
omega1, omega2, omega3 = omega0*np.asarray(cn), omega0*np.asarray(sn), omega0*np.asarray(dn)
```

### Cnoidal waves (KdV)

```python
# Python: KdV cnoidal wave solution
import elliptic, numpy as np
m = 0.99; kappa = 1.0
x = np.linspace(-10, 10, 500)
_, cn, _, _ = elliptic.ellipj(kappa * x, m)
u = 2 * m * kappa**2 * np.asarray(cn)**2
```

---

## Citation

```bibtex
@misc{elliptic,
  author       = {Moiseev I.},
  title        = {Elliptic functions for Matlab and Octave},
  year         = {2008},
  publisher    = {GitHub},
  howpublished = {\url{https://github.com/moiseevigor/elliptic}},
  doi          = {10.5281/zenodo.48264},
}
```

## References

- Abramowitz & Stegun, *Handbook of Mathematical Functions*, §16–18
- NIST DLMF §19, §22, §23 — https://dlmf.nist.gov
- Fukushima (2015), "Elliptic functions and elliptic integrals for celestial mechanics"
- Carlson (1995), "Numerical Computation of Real or Complex Elliptic Integrals"
- Bulirsch (1965), "Numerical computation of elliptic integrals and elliptic functions"
