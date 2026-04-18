
# Elliptic functions for MATLAB and Octave

![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-blue?logo=mathworks) ![Octave](https://img.shields.io/badge/Octave-6.0%2B-blue?logo=octave) ![GPU](https://img.shields.io/badge/GPU-OpenCL%20%7C%20CUDA-brightgreen) ![Parallel](https://img.shields.io/badge/CPU-multi--core-brightgreen)

Pure MATLAB/Octave implementation of elliptic integrals and functions — no external toolboxes required. Scales from a single core to multi-core CPU or GPU with a one-line config call.

## Installation

```matlab
git clone https://github.com/moiseevigor/elliptic.git
cd elliptic/matlab
setup    % adds matlab/src to the MATLAB/Octave path
```

## Functions

| Function | Description |
|---|---|
| `elliptic12(u, m)` | F(φ\|m), E(φ\|m), Jacobi Zeta Z(φ\|m) |
| `elliptic12i(u, m)` | Complex-argument F, E |
| `elliptic3(u, m, n)` | Π(φ,n\|m) — elliptic integral of the third kind |
| `ellipj(u, m)` | Jacobi sn, cn, dn, am |
| `ellipji(u, m)` | Jacobi sn, cn, dn for complex argument |
| `jacobiThetaEta(u, m)` | Jacobi theta and eta functions |
| `theta(n, u, m)` | Theta functions ϑ₁…ϑ₄ |
| `theta_prime(n, u, m)` | Derivatives of theta functions |
| `ellipticBDJ(phi, m, n)` | Associate incomplete B, D, J integrals |
| `ellipticBD(m)` | Complete associate B, D, S integrals |
| `carlsonRF(x,y,z)` | Carlson symmetric RF |
| `carlsonRD(x,y,z)` | Carlson symmetric RD |
| `carlsonRJ(x,y,z,p)` | Carlson symmetric RJ |
| `carlsonRC(x,y)` | Carlson degenerate RC |
| `jacobiEDJ(u, m, n)` | Jacobi-argument E_u, D_u, J_u |
| `cel(kc,p,a,b)` | Bulirsch generalised complete integral |
| `cel1(kc)`, `cel2(kc,a,b)`, `cel3(kc,p)` | Bulirsch wrappers |
| `weierstrassP(z,e1,e2,e3)` | Weierstrass ℘ function |
| `weierstrassZeta(z,e1,e2,e3)` | Weierstrass ζ function |
| `weierstrassSigma(z,e1,e2,e3)` | Weierstrass σ function |
| `agm(a, b)` | Arithmetic-geometric mean |

## Multi-core and GPU

```matlab
% Enable parallel computation (requires Parallel Computing Toolbox)
elliptic_config('parallel', true);

% Enable GPU acceleration (requires Parallel Computing Toolbox + CUDA)
elliptic_config('gpu', true);

% Disable
elliptic_config('parallel', false);
```

Speedups measured on 8-core machine: **2.3–2.5×** (parfor), up to **13×** (GPU, `elliptic3`).

## Key identities

```
F(φ|m) = B + D
E(φ|m) = B + (1−m)·D
Π(φ,n|m) = B + D + n·J

K(m) = RF(0, 1-m, 1)
E(m) = K - m/3 · RD(0, 1-m, 1)

cel1(kc) = K(1−kc²)
```

## Running tests

```matlab
cd matlab
setup
cd tests
test testElliptic12
test testElliptic3
test testEllipj
% ... etc
```

Or via Octave in CI: see `../.circleci/config.yml`.

## Benchmarks

```matlab
cd matlab/bench
bench        % serial benchmarks
bench_gpu    % GPU benchmarks (requires GPU)
```
