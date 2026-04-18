
# elliptic — Elliptic integrals and functions

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/moiseevigor/elliptic/tree/master.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/moiseevigor/elliptic/tree/master) [![DOI](https://zenodo.org/badge/5762/moiseevigor/elliptic.svg)](https://zenodo.org/badge/latestdoi/5762/moiseevigor/elliptic) [![tests](https://github.com/moiseevigor/elliptic/actions/workflows/python.yml/badge.svg)](https://github.com/moiseevigor/elliptic/actions/workflows/python.yml)

Full-precision implementations of elliptic integrals and functions for **MATLAB/Octave** and **Python** (NumPy · PyTorch · JAX).
No external toolboxes required at runtime — works standalone.

## Function coverage

| | MATLAB / Octave | Python |
|---|:---:|:---:|
| F(φ,m), E(φ,m), Z(φ,m) | ✓ | ✓ |
| Π(φ,m,n) — 3rd kind | ✓ | ✓ |
| Jacobi sn/cn/dn/am (real) | ✓ | ✓ |
| Jacobi sn/cn/dn (complex) | ✓ | ✓ |
| F, E, Z (complex argument) | ✓ | ✓ |
| Jacobi theta θ₁…θ₄ | ✓ | ✓ |
| Associate B, D, J integrals | ✓ | ✓ |
| Carlson RF / RD / RJ / RC | ✓ | ✓ |
| Bulirsch cel / cel1 / cel2 / cel3 | ✓ | ✓ |
| Weierstrass P / ζ / σ / ℘′ | ✓ | ✓ |
| Weierstrass invariants g₂, g₃ | ✓ | ✓ |
| Nome q(m) and m(q) | ✓ | ✓ |
| Inverse E(φ,m) | ✓ | ✓ |
| AGM | ✓ | ✓ |
| Ellipse arc length | ✓ | ✓ |
| Multi-core CPU | ✓ parfor | ✓ ProcessPoolExecutor |
| GPU | ✓ gpuArray | ✓ PyTorch CUDA / JAX |
| JAX jit / vmap | — | ✓ |

## MATLAB / Octave

```matlab
git clone https://github.com/moiseevigor/elliptic.git
cd elliptic/matlab
setup    % adds matlab/src to path
```

See [`matlab/README.md`](matlab/README.md) for full API reference, examples, and GPU/parallel usage.

Source in `matlab/src/`, tests in `matlab/tests/`, benchmarks in `matlab/bench/`.

## Python

```bash
pip install elliptic               # NumPy only
pip install "elliptic[torch]"      # + PyTorch (GPU via CUDA)
pip install "elliptic[jax]"        # + JAX (jit / vmap / TPU)
```

```python
import numpy as np, elliptic

phi = np.linspace(0, np.pi/2, 200)
F, E, Z = elliptic.elliptic12(phi, 0.7)
Pi      = elliptic.elliptic3(phi, 0.7, 0.3)
sn, cn, dn, am = elliptic.ellipj(np.linspace(0, 3, 100), 0.7)
B, D, J = elliptic.ellipticBDJ(phi, 0.7, 0.3)
Th3     = elliptic.theta(3, phi, 0.7)
arc     = elliptic.arclength_ellipse(5, 10)       # ≈ 48.44
```

See [`python/README.md`](python/README.md) for the full API, JAX/PyTorch examples, and comparison with scipy.

---

## Use cases and examples

Elliptic functions appear across physics, astronomy, geometry, and signal processing wherever periodic or nonlinear dynamics arise. The interactive examples at **https://moiseevigor.github.io/elliptic/examples** show the functions in action; the table below gives the mathematical mapping.

### Celestial mechanics

| Problem | Functions | Key formula |
|---|---|---|
| [Kepler orbit arc length](https://moiseevigor.github.io/elliptic/examples/arclength-celestial-mechanics/) | `elliptic12` | $s = a\,E(\phi\|e^2)$, eccentricity $e$ |
| Mercury perihelion precession | `elliptic12`, `ellipticBDJ` | Relativistic correction to elliptic arc |
| Gravitational potential of oblate spheroid | `carlsonRF`, `carlsonRD` | Closed-form via Carlson symmetric forms |
| Two-body problem period | `ellipk` → `elliptic12` | $T = 2\pi a^{3/2}/\sqrt{GM}$ generalised for eccentricity |

```python
# Kepler orbit: arc length from periapsis to true anomaly phi
import elliptic, numpy as np
a, e = 1.496e11, 0.0167          # Earth: semi-major axis, eccentricity
phi  = np.linspace(0, 2*np.pi, 1000)
_, E, _ = elliptic.elliptic12(phi, e**2)
arc_length = a * np.asarray(E)   # metres from periapsis
```

### Physical pendulum

| Problem | Functions | Key formula |
|---|---|---|
| Period for large amplitude | `ellipk` | $T = 4K(k^2)\sqrt{L/g}$, $k = \sin(\theta_0/2)$ |
| Phase portrait (exact trajectory) | `ellipj` | $\theta(t) = 2\arcsin(k\cdot\text{sn}(t\sqrt{g/L},k^2))$ |
| Escape time near separatrix | `ellipk` | $T\to\infty$ as $\theta_0\to\pi$ (logarithmic divergence) |

```python
import elliptic, numpy as np

g, L, theta0 = 9.81, 1.0, 2.5   # near separatrix (2.5 rad ≈ 143°)
k  = np.sin(theta0 / 2)
K  = float(elliptic.cel1(np.sqrt(1 - k**2)))   # = K(k^2)
T  = 4 * K * np.sqrt(L / g)                    # exact period

# exact trajectory
t  = np.linspace(0, T / 2, 500)
sn, _, _, _ = elliptic.ellipj(t * np.sqrt(g / L), k**2)
theta = 2 * np.arcsin(k * np.asarray(sn))
```

### Rigid body rotation (Euler / Euler-Poinsot)

| Problem | Functions | Key formula |
|---|---|---|
| Torque-free rotation of asymmetric top | `ellipj` | $\omega_1(t) = \omega_0\,\text{cn}(t\|m)$, $\omega_2\,\text{sn}$, $\omega_3\,\text{dn}$ |
| Polhode period | `ellipk` | Period $T = 4K(m)/\Omega$ |
| Attitude angle vs time | `elliptic12`, `jacobiEDJ` | Euler angles via Euler integrals |

```python
# Torque-free symmetric top — dn gives the precession rate
import elliptic, numpy as np
I1, I2, I3 = 1.0, 1.5, 2.0   # principal moments of inertia
omega0 = 1.0
m  = I2*(I3 - I2) / (I1*(I3 - I1)) * omega0**2  # modulus²
t  = np.linspace(0, 10, 1000)
sn, cn, dn, _ = elliptic.ellipj(t, m)
omega1 = omega0 * np.asarray(cn)
omega2 = omega0 * np.asarray(sn)
omega3 = omega0 * np.asarray(dn)
```

### Geometry — ellipse, lemniscate, geodesics

| Problem | Functions | Key formula |
|---|---|---|
| Ellipse perimeter / arc length | `arclength_ellipse` → `elliptic12` | $L = b\,E(2\pi\|e^2)$ |
| Lemniscate arc length (historic) | `elliptic12` | $s = \int_0^r (1-t^4)^{-1/2}dt$ — first known elliptic integral |
| Geodesic distance on ellipsoid | `ellipticBDJ` | Bessel / Helmert series via B, D |
| Cassini oval area | `elliptic12`, `ellipticBD` | Cassini ovals $r^2 = a^2\cos 2\phi$ |

```python
# Ellipse with semi-axes a=3, b=5: full perimeter
print(elliptic.arclength_ellipse(3, 5))   # ≈ 25.527

# Lemniscate of Bernoulli: half arc-length (the "lemniscate constant")
import scipy.special
omega = scipy.special.gamma(1/4)**2 / (4 * np.sqrt(np.pi))   # ≈ 2.622
```

### Signal processing and nonlinear waves

| Problem | Functions | Key formula |
|---|---|---|
| Jacobi elliptic filter (Cauer) | `ellipj`, `ellipk`, `nomeq` | Poles placed at $\text{sn}(k\cdot K/n\|m)$ |
| Cnoidal waves (KdV solitons) | `ellipj` | $u(x,t) = A\,\text{cn}^2(\kappa(x-vt)\|m)$ |
| Nonlinear lattice vibrations | `ellipj`, `ellipticBDJ` | Toda lattice exact solution |

```python
# Cnoidal wave: KdV 1-soliton in the cnoidal limit
import elliptic, numpy as np
m = 0.99                                  # near soliton limit
kappa = 1.0
x = np.linspace(-10, 10, 500)
sn, cn, dn, _ = elliptic.ellipj(kappa * x, m)
u_cnoidal = 2 * m * kappa**2 * np.asarray(cn)**2  # cnoidal wave amplitude
```

### Lattice sums, number theory, and modular forms

| Problem | Functions | Key formula |
|---|---|---|
| Theta function identities | `theta` | $\theta_3(0,q)^4 = \theta_2(0,q)^4 + \theta_4(0,q)^4$ (Jacobi) |
| Counting representations as sums of squares | `theta` | $r_2(n) = 4(\theta_3(0,q)^2)$ coefficient |
| AGM algorithm for π | `agm` | $\pi = 4\,\text{agm}(1,\sqrt{2}/2)^2 / \sum 2^n c_n^2$ (Gauss) |
| Weierstrass elliptic curve cryptography | `weierstrassP`, `weierstrassInvariants` | Group law on $y^2 = 4x^3 - g_2 x - g_3$ |

```python
# Jacobi identity: θ₃⁴ = θ₂⁴ + θ₄⁴  (holds for all q)
import elliptic, numpy as np
q = 0.3
th2 = float(elliptic.theta(2, 0.0, float(elliptic.inversenomeq(np.array([q]))[0])))
th3 = float(elliptic.theta(3, 0.0, float(elliptic.inversenomeq(np.array([q]))[0])))
th4 = float(elliptic.theta(4, 0.0, float(elliptic.inversenomeq(np.array([q]))[0])))
print(abs(th3**4 - th2**4 - th4**4))   # < 1e-13
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
- Boyd (2012), "Numerical inversion of the incomplete elliptic integral of the second kind"
