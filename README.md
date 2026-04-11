
# Elliptic functions for Matlab and Octave

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/moiseevigor/elliptic/tree/master.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/moiseevigor/elliptic/tree/master) [![DOI](https://zenodo.org/badge/5762/moiseevigor/elliptic.svg)](https://zenodo.org/badge/latestdoi/5762/moiseevigor/elliptic) ![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-blue?logo=mathworks) ![Octave](https://img.shields.io/badge/Octave-6.0%2B-blue?logo=octave) ![GPU](https://img.shields.io/badge/GPU-OpenCL%20%7C%20CUDA-brightgreen) ![Parallel](https://img.shields.io/badge/CPU-multi--core-brightgreen)

The Matlab/Octave implementation of [Elliptic integrals of three types](http://en.wikipedia.org/wiki/Elliptic_integral), [Jacobi's elliptic functions](http://en.wikipedia.org/wiki/Jacobi%27s_elliptic_functions) and [Jacobi theta functions](http://en.wikipedia.org/wiki/Theta_function) of four types with their derivatives.

> **Goal:** Pure MATLAB/Octave scripts for elliptic integrals and functions — **no external libraries** (no Maple, no Mathematica toolboxes). Every function accepts tensors of any shape, nearly all have complex-argument variants, and the library scales from a single core to multi-core CPU or GPU with a one-line config call.

| Feature | Detail |
|---|---|
| **Functions** | `elliptic12`, `elliptic3`, `ellipj`, `jacobiThetaEta`, `theta`, `theta_prime`, inverses |
| **Input shapes** | Scalars, vectors, matrices — tensors of any shape |
| **Complex support** | `ellipj` → `ellipji` · `elliptic12` → `elliptic12i` |
| **Multi-core CPU** | `elliptic_config('parallel', true)` — **2.3–2.5×** speedup on 8 cores |
| **GPU** | `elliptic_config('gpu', true)` — up to **13×** (`elliptic3` via OpenCL / CUDA) |
| **Platforms** | MATLAB R2019b+ · GNU Octave 6.0+ |

# Installation

Clone the repository and run the setup script to add the library to your path:

```matlab
git clone https://github.com/moiseevigor/elliptic.git
cd elliptic
setup    % adds src/ to the MATLAB/Octave path
```

All source files are located in the `src/` directory. Tests are in the `tests/` directory.

# Citations and references

If you've used any of the routines in this package please cite and support the effort. Here is the example of the BibTeX entry

```
@misc{elliptic,
  author = {Moiseev I.},
  title = {Elliptic functions for Matlab and Octave},
  year = {2008},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/moiseevigor/elliptic}},
  commit = {98181c4c0d8992746bcc6bea75740bb11b74b51b},
  doi = {10.5281/zenodo.48264},
  url = {http://dx.doi.org/10.5281/zenodo.48264}
}
```

or simply

```
Moiseev I., Elliptic functions for Matlab and Octave, (2008), GitHub repository, DOI: http://dx.doi.org/10.5281/zenodo.48264
```

# Contents of the package

  - [Elliptic Functions](#elliptic-functions)
    - [ELLIPJ: Jacobi's elliptic functions](#ellipj-jacobis-elliptic-functions)
    - [ELLIPJI: Jacobi's elliptic functions of the complex argument](#ellipji-jacobis-elliptic-functions-of-the-complex-argument)
    - [JACOBITHETAETA: Jacobi's Theta and Eta Functions](#jacobithetaeta-jacobis-theta-and-eta-functions)
    - [THETA: Theta Functions of Four Types](#theta-theta-functions-of-four-types)
    - [THETA_PRIME: Theta Functions and Their Derivatives](#theta_prime-theta-functions-and-their-derivatives)
  - [Elliptic Integrals](#elliptic-integrals)
    - [ELLIPTIC12: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function](#elliptic12-incomplete-elliptic-integrals-of-the-first-second-kind-and-jacobis-zeta-function)
    - [ELLIPTIC12I: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function of the complex argument](#elliptic12i-incomplete-elliptic-integrals-of-the-first-second-kind-and-jacobis-zeta-function-of-the-complex-argument)
    - [ELLIPTIC3: Incomplete Elliptic Integral of the Third Kind](#elliptic3-incomplete-elliptic-integral-of-the-third-kind)
    - [ELLIPTIC123: Complete and Incomplete Elliptic Integrals of the First, Second, and Third Kind](#elliptic123-complete-and-incomplete-elliptic-integrals-of-the-first-second-and-third-kind)
    - [INVERSELLIPTIC2: INVERSE Incomplete Elliptic Integrals of the Second Kind](#inverselliptic2-inverse-incomplete-elliptic-integrals-of-the-second-kind)

  - [Weierstrass's elliptic functions (in development)](#weierstrasss-elliptic-functions-in-development)
  - [Elliptic Related Functions](#elliptic-related-functions)
    - [AGM: Arithmetic Geometric Mean](#agm-arithmetic-geometric-mean)
    - [NOMEQ: The Value of Nome q = q(m)](#nomeq-the-value-of-nome-q--qm)
    - [INVERSENOMEQ: The Value of Nome m = m(q)](#inversenomeq-the-value-of-nome-m--mq)
  - [Multi-Core Parallel Execution](#multi-core-parallel-execution)
  - [GPU Acceleration](#gpu-acceleration)
  - [Contributors](#contributors)
  - [References](#references)

# Elliptic Functions

**The Jacobi's elliptic functions** are a set of basic elliptic functions, and auxiliary theta functions, that have historical importance with also many features that show up important structure, and have direct relevance to some applications (e.g. the equation of the pendulum). They also have useful analogies to the functions of trigonometry, as indicated by the matching notation `SN` for `SIN`. They are not the simplest way to develop a general theory, as now seen: that can be said for the Weierstrass elliptic functions. They are not, however, outmoded. They were introduced by Carl Gustav Jakob Jacobi, around 1830.

**Theta functions** are special functions of several complex variables. They are important in several areas, including the theories of abelian varieties and moduli spaces, and of quadratic forms. They have also been applied to soliton theory. When generalized to a Grassmann algebra, they also appear in quantum field theory, specifically string theory and D-branes.

## ELLIPJ: Jacobi's elliptic functions

`ELLIPJ` evaluates the [Jacobi's elliptic functions](http://en.wikipedia.org/wiki/Jacobi%27s_elliptic_functions) and [Jacobi's amplitude](http://mathworld.wolfram.com/JacobiAmplitude.html).

`[Sn,Cn,Dn,Am] = ELLIPJ(U,M)` returns the values of the Jacobi elliptic functions `SN`, `CN`, `DN` and `AM` evaluated for corresponding elements of argument U and parameter M.  The arrays U and M must be of the same size (or either can be scalar).  As currently implemented, M is limited to `0 <= M <= 1`.

*General definition:*
```
u = Integral(1/sqrt(1-m*sin(theta)^2), 0, phi);
Sn(u) = sin(phi);
Cn(u) = cos(phi);
Dn(u) = sqrt(1-m*sin(phi)^2);
```

_Depends on_  `AGM`, `ELLIPKE`.<br>
_Used by_ `THETA`.<br>
_See also_ `ELLIPKE`.


## ELLIPJI: Jacobi's elliptic functions of the complex argument

`ELLIPJI` evaluates the Jacobi elliptic functions of complex phase `U`.

`[Sni,Cni,Dni] = ELLIPJ(U,M)` returns the values of the Jacobi elliptic functions `SNI`, `CNI` and `DNI` evaluated for corresponding  elements of argument `U` and parameter `M`. The arrays `U` and `M` must  be of the same size (or either can be scalar).  As currently implemented, `M` is real and limited to `0 <= M <= 1`.


### Example:

```
[phi1,phi2] = meshgrid(-pi:3/20:pi, -pi:3/20:pi);
phi = phi1 + phi2*i;
[Sni,Cni,Dni]= ellipji(phi, 0.99);
```

_Depends on_ `AGM`, `ELLIPJ`, `ELLIPKE`<br>
_See also_ `ELLIPTIC12`, `ELLIPTIC12I`

## JACOBITHETAETA: Jacobi's Theta and Eta Functions

`JACOBITHETAETA` evaluates Jacobi's theta and eta functions.

`[Th, H] = JACOBITHETAETA(U,M)` returns the values of the Jacobi's theta and eta elliptic functions `TH` and `H` evaluated for corresponding elements of argument `U` and parameter `M`.  The arrays `U` and `M` must be the same size (or either can be scalar).  As currently implemented, `M` is real and limited to `0 <= M <= 1`.

### Example:
```
[phi,alpha]= meshgrid(0:5:90, 0:2:90);
[Th, H] = jacobiThetaEta(pi/180*phi, sin(pi/180*alpha).^2);
```

*Depends on* `AGM`, `ELLIPJ`, `ELLIPKE`<br>
*See also* `ELLIPTIC12`, `ELLIPTIC12I`, `THETA`

## THETA: Theta Functions of Four Types

`THETA` evaluates theta functions of four types.

`Th = THETA(TYPE,V,M)` returns values of theta functions
evaluated for corresponding values of argument `V` and parameter `M`. `TYPE` is a type of the theta function, there are four numbered types. The arrays `V` and `M` must be the same size (or either can be scalar). As currently implemented, `M` is limited to `0 <= M <= 1`.

`Th = THETA(TYPE,V,M,TOL)` computes the theta and eta elliptic functions to the accuracy `TOL` instead of the default `TOL = EPS`.

The parameter `M` is related to the nome `Q` as `Q = exp(-pi*K(1-M)/K(M))`. Some definitions of the Jacobi's elliptic functions use the modulus `k` instead of the parameter `m`.  They are related by `m = k^2`.

### Example:
```
[phi,alpha] = meshgrid(0:5:90, 0:2:90);
Th1 = theta(1, pi/180*phi, sin(pi/180*alpha).^2);
Th2 = theta(2, pi/180*phi, sin(pi/180*alpha).^2);
Th3 = theta(3, pi/180*phi, sin(pi/180*alpha).^2);
Th4 = theta(4, pi/180*phi, sin(pi/180*alpha).^2);
```

_Depends on_ `AGM`, `ELLIPJ`, `ELLIPKE`, `JACOBITHETAETA`<br>
_See also_ `ELLIPTIC12`, `ELLIPTIC12I`

## THETA_PRIME: Theta Functions and Their Derivatives

`THETA_PRIME` evaluates theta functions and their derivatives with respect to the argument.

`[TH, THP] = THETA_PRIME(J, Z, M)` returns values of the Jacobi theta function `TH` and its derivative `THP` with respect to the argument `Z`. `J` is the type of theta function (1, 2, 3, or 4), `Z` is the argument, and `M` is the parameter (0 <= M <= 1).

`[TH, THP] = THETA_PRIME(J, Z, M, TOL)` computes the theta function and its derivative to the accuracy `TOL` instead of the default `TOL = EPS`.

The arrays `Z` and `M` must be the same size (or either can be scalar). As currently implemented, `M` is limited to `0 <= M <= 1`.

The derivatives are computed using the relation:
```
theta'_j(z,m) = theta_j(z,m) * (2K/pi) * (Z + delta_j)
```
where `K` is the complete elliptic integral, `Z` is the Jacobi zeta function, and `delta_j` depends on the theta function type.

### Example:
```
j = 1; z = 0.5; m = 0.8;
[th, thp] = theta_prime(j, z, m);
```

**Note:** To reproduce results in Mathematica, use the `inversenomeq` function to convert from Mathematica's nome parameter to the parameter `m`:
```
[th_4, thp_4] = theta_prime(4, z, inversenomeq(0.1));
```

_Depends on_ `THETA`, `ELLIPJ`, `ELLIPTIC12`, `ELLIPKE`<br>
_See also_ `THETA`, `JACOBITHETAETA`, `INVERSENOMEQ`

# Elliptic Integrals

Elliptic integrals originally arose in connection with the problem of giving the arc length of an ellipse. They were first studied by Giulio Fagnano and Leonhard Euler. Modern mathematics defines an elliptic integral as any function f which can be expressed in the form

```
f(x) = Integral(R(t,P(t), c, x)dt,
```

where `R` is a rational function of its two arguments, `P` is the square root of a polynomial of degree `3` or `4` with no repeated roots, and `c` is a constant.
In general, elliptic integrals cannot be expressed in terms of elementary functions. Exceptions to this general rule are when `P` has repeated roots, or when `R(x,y)` contains no odd powers of `y`. However, with the appropriate reduction formula, every elliptic integral can be brought into a form that involves integrals over rational functions and the three canonical forms (i.e. the elliptic integrals of the first, second and third kind).

## ELLIPTIC12: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function

`ELLIPTIC12` evaluates the value of the Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function.

`[F,E,Z] = ELLIPTIC12(U,M,TOL)` uses the method of the Arithmetic-Geometric Mean and Descending Landen Transformation described in [1](#references) Ch. 17.6, to determine the value of the Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function (see [1](#references), [2](#references)).

### General definition:

```
F(phi,m) = int(1/sqrt(1-m*sin(t)^2), t=0..phi);
E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
Z(phi,m) = E(u,m) - E(m)/K(m)*F(phi,m).
```

Tables generating code (see [1](#references), pp. 613-621):

```
[phi,alpha] = meshgrid(0:5:90, 0:2:90);                  % modulus and phase in degrees
[F,E,Z] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);  % values of integrals
```

*Depends on* `AGM`<br>
*See also* `ELLIPKE`, `ELLIPJ`, `ELLIPTIC12I`, `ELLIPTIC3`, `THETA`.

## ELLIPTIC12I: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function of the complex argument

`ELLIPTIC12i` evaluates the Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function for the complex value of phase `U`. Parameter `M` must be in the range `0 <= M <= 1`.

`[Fi,Ei,Zi] = ELLIPTIC12i(U,M,TOL)` where `U` is a complex phase in radians, `M` is the real parameter and `TOL` is the tolerance (optional). Default value for the tolerance is `eps = 2.220e-16`.

`ELLIPTIC12i` uses the function `ELLIPTIC12` to evaluate the values of corresponding integrals.

### Example:

```
[phi1,phi2] = meshgrid(-2*pi:3/20:2*pi, -2*pi:3/20:2*pi);
phi = phi1 + phi2*i;
[Fi,Ei,Zi] = elliptic12i(phi, 0.5);
```

_Depends on_ `ELLIPTIC12`, `AGM`<br>
_See also_ `ELLIPKE`, `ELLIPJ`, `ELLIPTIC3`, `THETA`.

## ELLIPTIC3: Incomplete Elliptic Integral of the Third Kind

`ELLIPTIC3` evaluates incomplete elliptic integral of the third kind `Pi = ELLIPTIC3(U,M,C)` where `U` is a phase in radians, `0 < M < 1` is the module and `0 < C < 1` is a parameter.

`ELLIPTIC3` uses Gauss-Legendre 10 points quadrature template described in [3] to determine the value of the Incomplete Elliptic Integral of the Third Kind (see [1, 2]).

### General definition:

```
Pi(u,m,c) = int(1/((1 - c*sin(t)^2)*sqrt(1 - m*sin(t)^2)), t=0..u)
```

Tables generating code ([1](#references), pp. 625-626):
```
[phi,alpha,c] = meshgrid(0:15:90, 0:15:90, 0:0.1:1);
Pi = elliptic3(pi/180*phi, sin(pi/180*alpha).^2, c);  % values of integrals
```

## ELLIPTIC123: Complete and Incomplete Elliptic Integrals of the First, Second, and Third Kind

`ELLIPTIC123` is a wrapper around the different elliptic integral functions, providing a unified interface and greater range of input parameters. (Unlike ELLIPKE, ELLIPTIC12 and ELLIPTIC3, which all require a phase between zero and pi/2 and a parameter between zero and one.)

`[F,E] = ELLIPTIC123(m)` -- complete Elliptic Integrals of the first and second kind.

`[F,E] = ELLIPTIC123(b,m)` -- incomplete Elliptic Integrals of the first and second kind.

`[F,E,PI] = ELLIPTIC123(m,n)` -- complete Elliptic Integrals of the first to third kind.

`[F,E,PI] = ELLIPTIC123(b,m,n)` -- incomplete Elliptic Integrals of the first to third kind.

The order of the input arguments has been chosen to be consistent with the pre-existing `elliptic12` and `elliptic3` functions.

This function is still under development and its results are not always well-defined or even able to be calculated (especially for the third elliptic integral with `n>1`). Please see the documentation for further details.

## INVERSELLIPTIC2: INVERSE Incomplete Elliptic Integrals of the Second Kind

`INVERSELLIPTIC2` evaluates the value of the INVERSE Incomplete Elliptic Integrals of the Second Kind.

`INVERSELLIPTIC2` uses the method described by Boyd J. P. to determine the value of the inverse Incomplete Elliptic Integrals of the Second Kind using the "Empirical" initialization to the Newton's iteration method [7](#references).

Elliptic integral of the second kind:

```
E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
```

"Empirical" initialization [7](#references):

```
T0(z,m) = pi/2 + sqrt(r)/(theta - pi/2)
```

where

```
z in (-E(pi/2,m), E(pi/2,m)) x (0, 1) - value of the entire parameter space
r = sqrt((1-m)^2 + zeta^2)
zeta = 1 - z/E(pi/2,m)
theta = atan((1 - m)/zeta)
```

Example:
```
% modulus and phase in degrees
[phi,alpha] = meshgrid(0:5:90, 0:2:90);
% values of integrals
[F,E] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);
% values of inverse
invE = inverselliptic2(E, sin(pi/180*alpha).^2);
% the difference between phase phi and invE should close to zero
phi - invE * 180/pi
```

# Weierstrass's elliptic functions (in development)

Weierstrass's elliptic functions are elliptic functions that take a particularly simple form (cf Jacobi's elliptic functions); they are named for Karl Weierstrass. This class of functions are also referred to as p-functions and generally written using the symbol &#8472; (a stylised letter p called Weierstrass p).

The Weierstrass elliptic function can be defined in three closely related ways, each of which possesses certain advantages. One is as a function of a complex variable z and a lattice &#923; in the complex plane. Another is in terms of z and two complex numbers &#969;1 and &#969;2 defining a pair of generators, or periods, for the lattice. The third is in terms of z and a modulus &#964; in the upper half-plane. This is related to the previous definition by `tau = omega2 / omega1`, which by the conventional choice on the pair of periods is in the upper half-plane. Using this approach, for fixed z the Weierstrass functions become modular functions of &#964;.


# Elliptic Related Functions

## AGM: Arithmetic Geometric Mean

`AGM` calculates the [Arithmetic Geometric Mean](http://en.wikipedia.org/wiki/Arithmetic-geometric_mean) of `A` and `B` (see [1](#references)).

`[A,B,C,N]= AGM(A0,B0,C0,TOL)` carry out the process of the arithmetic geometric mean, starting with a given positive numbers triple `(A0, B0, C0)` and returns in
`(A, B, C)` the generated sequence. `N` is a number of steps (returns in the value `uint32`).

The general scheme of the procedure:

```
A(i) = 1/2*( A(i-1)+B(i-1) );     A(0) = A0;
B(i) = sqrt( A(i-1)*B(i-1) );     B(0) = B0;
C(i) = 1/2*( A(i-1)-B(i-1) );     C(0) = C0;
```

Stop at the `N`-th step when `A(N) = B(N)`, i.e., when `C(N) = 0`.

_Used by_  `ELLIPJ` and `ELLIPTIC12`.<br>_See also_ `ELLIPKE`, `ELLIPTIC3`, `THETA`.

## NOMEQ: The Value of Nome `q = q(m)`

`NOMEQ` gives the value of Nome `q = q(m)`.

Nome `Q = nomeq(M,TOL)`, where `0<=M<=1` is the module and `TOL` is the tolerance (optional). Default value for the tolerance is `eps = 2.220e-16`.

*Used by*  `ELLIPJ`.<br>
*Depends on* `ELLIPKE`<br>
*See also* `ELLIPTIC12I`, `ELLIPTIC3`, `THETA`.

## INVERSENOMEQ: The Value of Nome `m = m(q)`

`INVERSENOMEQ` gives the value of Nome `m = m(q)`.

`M = inversenomeq(q)`, where `Q` is the Nome of q-series.

**WARNING**. The function `INVERSENOMEQ` does not return correct values of `M` for `Q < 0.00001` or `Q > 0.76`, because of computer precision limitation. The function `NomeQ(m)` has an essential singularity at `M = 1`, so it cannot be inverted at this point and is very hard to invert in its neighbourhood.

More precisely:

```
nomeq(1) = 1
nomeq(1-eps) = 0.77548641878026
```

### Example:

```
nomeq(inversenomeq([0.3 0.4 0.5 0.6 0.7 0.8]))
warning: The function INVERSENOMEQ does not return
correct values of M for Q < 0.00001 and Q > 0.76, because of computer precision limitation.
warning: called from
    inversenomeq at line 54 column 5

ans =

   0.3000   0.4000   0.5000   0.5948   0.7000   0.7447
```

*Used by*  `ELLIPJ`.<br>
*Depends on* `ELLIPKE`<br>
*See also* `ELLIPTIC12I`, `ELLIPTIC3`, `THETA`.

# Project structure

```
elliptic/
  src/                      % source files (added to path by setup.m)
    elliptic12.m
    elliptic12i.m
    elliptic123.m
    elliptic3.m
    ellipj.m
    ellipji.m
    jacobiThetaEta.m
    theta.m
    theta_prime.m
    agm.m
    nomeq.m
    inversenomeq.m
    inverselliptic2.m
    arclength_ellipse.m
    uniquetol_compat.m
    elliptic_config.m       % get/set library config (parallel, gpu, chunk_size)
    get_nworkers.m          % detect available parallel workers
    has_gpu.m               % detect GPU availability (MATLAB PCT or Octave ocl)
    par_worker.m            % generic parallel worker (Octave parcellfun)
  tests/                    % test files
    testElliptic12.m
    testElliptic3.m
    testEllipj.m
    testThetaPrime.m
    testAgm.m
    testJacobiThetaEta.m
    testGpu.m               % GPU correctness tests (skipped when no GPU)
  docs/                     % documentation
    GPU.md                  % GPU setup guide (MATLAB + Octave/OpenCL)
    benchmark_gpu.md        % three-way benchmark results with hardware utilisation
    benchmark_comparison.md % before/after vectorization comparison
  setup.m                   % run this to add src/ to path
  bench_gpu.m               % three-way benchmark: serial / parallel / GPU
  bench_gpu_results.csv     % raw benchmark data
```

# Running tests

From root folder run `setup` first, then `runtests`, or check the CI build status [![CircleCI](https://dl.circleci.com/status-badge/img/gh/moiseevigor/elliptic/tree/master.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/moiseevigor/elliptic/tree/master)

### MATLAB

```matlab
setup
runtests("tests")
```

### Octave

```matlab
setup
cd tests
test testElliptic12
test testElliptic3
test testEllipj
test testThetaPrime
test testAgm
test testJacobiThetaEta
test testGpu          % skipped automatically when no GPU present
```

### Expected output (Octave)

```
PASSES 5 out of 5 tests    % testElliptic12
PASSES 6 out of 6 tests    % testElliptic3
PASSES 5 out of 5 tests    % testEllipj
PASSES 4 out of 4 tests    % testAgm
PASSES 7 out of 7 tests    % testJacobiThetaEta
PASSES 19 out of 19 tests  % testThetaPrime
```

**Note:** Timing thresholds in the benchmark tests are set conservatively to accommodate CI machines. All correctness tests pass on both MATLAB and Octave.

# Multi-Core Parallel Execution

Parallelism is built into the main functions — no separate API needed. When enabled, inputs larger than the chunk size threshold are automatically split across CPU cores. The API stays the same; only two setup steps are required.

### MATLAB

Requires [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html).

```matlab
setup;
parpool(8);                          % 1. start a parallel pool (once per session)
elliptic_config('parallel', true);   % 2. enable parallel mode

% All calls now run in parallel automatically
[F,E,Z] = elliptic12(u, m);
[Sn,Cn,Dn,Am] = ellipj(u, m);
[Th,H] = jacobiThetaEta(u, m);
Pi = elliptic3(u, m, c);
```

### Octave

Requires the `parallel` package. Install once:

```bash
# Install build dependencies (Debian/Ubuntu)
sudo apt-get install -y liboctave-dev

# Install Octave packages
octave --eval "pkg install -forge struct"
octave --eval "pkg install -forge parallel"
```

Then in Octave:

```matlab
pkg load parallel;                   % 1. load the parallel package
setup;                               %    add src/ to path
elliptic_config('parallel', true);   % 2. enable parallel mode

% All calls now run in parallel automatically
[F,E,Z] = elliptic12(u, m);
[Sn,Cn,Dn,Am] = ellipj(u, m);
[Th,H] = jacobiThetaEta(u, m);
Pi = elliptic3(u, m, c);
```

### How it works

1. `elliptic_config('parallel', true)` enables parallel mode
2. Each function checks `get_nworkers()` — returns available cores (Octave) or pool size (MATLAB)
3. If workers > 1 and input size >= `chunk_size`, the input is split into chunks
4. Chunks are computed in parallel via `parfor` (MATLAB) or `parcellfun` (Octave)
5. Results are reassembled — bit-exact same output as serial

When parallel is disabled (default), overhead is zero — just one `if` branch skipped.

### Configuration

```matlab
elliptic_config('parallel', true);     % enable parallel (default: false)
elliptic_config('chunk_size', 20000);  % min elements per chunk (default: 10000)
elliptic_config('parallel', false);    % disable and return to serial
```

### Benchmark results (8-core i7-7700K @ 4.2 GHz, Octave 6.4)

Measured with `bench_gpu.m`, 3 repetitions, minimum wall-clock time. Parallel uses 8 workers with chunk size `ceil(N / nWorkers)`.

| Function | N | Serial (s) | Parallel (s) | Speedup |
|---|---|---|---|---|
| `elliptic12` | 1 M | 1.296 | 0.529 | **2.4×** |
| `ellipj` | 1 M | 1.089 | 0.465 | **2.3×** |
| `elliptic3` | 4 M | 2.668 | 1.167 | **2.3×** |
| `jacobiThetaEta` | 1 M | 1.460 | 0.583 | **2.5×** |

The parallel path incurs ~35 ms dispatch overhead (`parcellfun` spawn cost in Octave). It pays off for N > 50 000–100 000 elements; below that use serial (the default). Full results across seven input sizes: [`docs/benchmark_gpu.md`](docs/benchmark_gpu.md).

Run the benchmark yourself:
```matlab
pkg load parallel ocl;               % Octave only
setup;
bench_gpu;
```

# GPU Acceleration

All four core functions support GPU-accelerated evaluation via an opt-in configuration flag. The same code paths work in both **MATLAB** (Parallel Computing Toolbox / CUDA) and **Octave** (ocl Forge package / OpenCL).

### Enable GPU mode

```matlab
% Octave: load the package first
pkg load ocl          % Octave only

setup;
elliptic_config('gpu', true);   % enable GPU for all functions

% All calls now run on the GPU automatically
[F, E, Z] = elliptic12(u, m);
Pi        = elliptic3(u, m, c);
[sn, cn]  = ellipj(u, m);
[Th, H]   = jacobiThetaEta(u, m);

elliptic_config('gpu', false);  % disable and return to CPU
```

GPU and parallel modes are mutually exclusive — GPU takes precedence when both are enabled.

### Benchmark results (GTX 1080 Ti, Octave 6.4 / OpenCL 3.0)

GPU time includes host↔device memory transfer.

| Function | Serial (s) | Parallel 8-core (s) | GPU (s) | GPU speedup | GPU Mpts/s |
|---|---|---|---|---|---|
| `elliptic12` | 1.296 | 0.529 | 0.280 | **4.6×** | 3.6 |
| `ellipj` | 1.089 | 0.465 | 0.285 | **3.8×** | 3.5 |
| `elliptic3` | 2.668 | 1.167 | 0.199 | **13.4×** | 20.1 |
| `jacobiThetaEta` | 1.460 | 0.583 | 0.494 | **3.1×** | 2.0 |

`elliptic3` benefits most from the GPU (13×) because it uses pure Gauss-Legendre quadrature with no sequential AGM convergence loop — all 10 iterations execute fully on the GPU. AGM-based functions require a `gather()` call each iteration to check convergence on the CPU, limiting GPU occupancy to 27–41%.

### When to use GPU

| Scenario | Recommendation |
|---|---|
| N < 50 000 | Serial CPU — transfer overhead dominates |
| 50 000 ≤ N < 200 000 | GPU starts to break even (2–4×) |
| N ≥ 500 000 | GPU recommended (3–13×) |
| `elliptic3`, any large N | GPU strongly recommended (10–13×) |

Full installation guide, troubleshooting, and memory budget: [`docs/GPU.md`](docs/GPU.md)  
Full benchmark tables with hardware utilisation: [`docs/benchmark_gpu.md`](docs/benchmark_gpu.md)

# Compatibility

The library is tested on both **MATLAB** and **GNU Octave** via CircleCI:
- **MATLAB**: Full support via the MathWorks MATLAB orb
- **Octave**: Full support via `uniquetol_compat.m` compatibility wrapper (Octave lacks the built-in `uniquetol` function)

# Contributors

Contributors

- [@moiseevigor](http://github.com/moiseevigor/) - maintainer (Igor Moiseev)
- [@drbitboy](https://github.com/drbitboy) (Brian Carcich)
- [@wspr](https://github.com/wspr) (Will Robertson)


# References

  1. [NIST Digital Library of Mathematical Functions](http://dlmf.nist.gov/)
  1. M. Abramowitz and I.A. Stegun, "[Handbook of Mathematical Functions](https://personal.math.ubc.ca/~cbm/aands/)" Dover Publications", 1965, Ch. 17.1 - 17.6.
  1. D. F. Lawden, "[Elliptic Functions and Applications](https://www.amazon.com/Elliptic-Functions-Applications-Mathematical-Sciences/dp/0387969659)" Springer-Verlag, vol. 80, 1989
  1. S. Zhang, J. Jin, "[Computation of Special Functions](https://www.amazon.com/Computation-Special-Functions-Shanjie-Zhang/dp/0471119636)" (Wiley, 1996).
  1. B. Carlson, "[Three Improvements in Reduction and Computation of Elliptic Integrals](https://pmc.ncbi.nlm.nih.gov/articles/PMC4861378/)", J. Res. Natl. Inst. Stand. Technol. 107 (2002) 413-418.
  1. N. H. Abel, "[Studies on Elliptic Functions](https://old.maa.org/sites/default/files/images/upload_library/1/abeltranslation.pdf)", english translation from french by Marcus Emmanuel Barnes. Original "Recherches sur les fonctions elliptiques", Journal fr die reine und angewandte Mathematik, Vol. 2, 1827. pp. 101-181.
  1. B. C. Berndt, H. H. Chan, S.-S. Huang, "[Incomplete Elliptic Integrals in Ramanujan's Lost Notebook](https://faculty.math.illinois.edu/~berndt/publications.html)" in q-series from a Contemporary Perspective, M. E. H. Ismail and D. Stanton, eds., Amer. Math. Soc., 2000, pp. 79-126.
  1. J. P. Boyd, "[Numerical, Perturbative and Chebyshev Inversion of the Incomplete Elliptic Integral of the Second Kind](https://doi.org/10.1016/j.amc.2011.12.021)", Applied Mathematics and Computation, vol. 218, no. 13, 2012, pp. 7005-7013.
