# Elliptic functions for Matlab and Octave

[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/moiseevigor/elliptic?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge) [![DOI](https://zenodo.org/badge/5762/moiseevigor/elliptic.svg)](https://zenodo.org/badge/latestdoi/5762/moiseevigor/elliptic)

The [Matlab](http://www.mathworks.com/) script implementations of [Elliptic integrals of three types](http://en.wikipedia.org/wiki/Elliptic_integral), [Jacobi's elliptic functions](http://en.wikipedia.org/wiki/Jacobi%27s_elliptic_functions) and [Jacobi theta functions](http://en.wikipedia.org/wiki/Theta_function) of four types.

The main *GOAL* of the project is to provide the natural Matlab scripts *WITHOUT* external library calls like Maple and others. All scripts are developed to accept tensors as arguments and almost all of them have their complex versions. Performance and complete control on the execution are the main features.

# Interactive Examples

Explore the mathematics behind the functions with live, browser-based visualisations — all computations run client-side using the Carlson duplication algorithm ported to JavaScript.

| Example | Topics |
|---------|--------|
| [**Arc Length & Celestial Mechanics**](examples/arclength-celestial-mechanics/) | Ellipse arc length · Keplerian orbits · Kepler's equation · Associate integrals B, D, J · Carlson convergence |
| [**Physical Pendulum & Exact Trajectories**](examples/physical-pendulum/) | Exact period via K(k²) · Jacobi sn/cn/dn · Phase portrait · Separatrix · Period divergence · Compound pendula |
| [**GPU Asteroid Swarm**](examples/gpu-asteroid-swarm/) | Real SBDB catalog · Batched Keplerian propagation · Incomplete E(φ\|m) arc lengths · PyTorch CUDA / JAX / gpuArray · Horizons validation |

[**Browse all interactive examples →**](examples/)

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
  - [Elliptic Integrals](#elliptic-integrals)
    - [ELLIPTIC12: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function](#elliptic12-incomplete-elliptic-integrals-of-the-first-second-kind-and-jacobis-zeta-function)
    - [ELLIPTIC12I: Incomplete Elliptic Integrals of the First, Second Kind and Jacobi's Zeta Function of the complex argument](#elliptic12i-incomplete-elliptic-integrals-of-the-first-second-kind-and-jacobis-zeta-function-of-the-complex-argument)
    - [ELLIPTIC3: Incomplete Elliptic Integral of the Third Kind](#elliptic3-incomplete-elliptic-integral-of-the-third-kind)
    - [ELLIPTIC123: Complete and Incomplete Elliptic Integrals of the First, Second, and Third Kind](#elliptic123-complete-and-incomplete-elliptic-integrals-of-the-first-second-and-third-kind)
    - [INVERSELLIPTIC2: INVERSE Incomplete Elliptic Integrals of the Second Kind](#inverselliptic2-inverse-incomplete-elliptic-integrals-of-the-second-kind)
  - [Weierstrass's Elliptic Functions](#weierstrasss-elliptic-functions)
    - [WEIERSTRASSINVARIANTS: Lattice Invariants](#weierstrassinvariants-lattice-invariants)
    - [WEIERSTRASSP: Weierstrass P-function](#weierstrassp-weierstrass-p-function)
    - [WEIERSTRASSPPRIME: Derivative of the P-function](#weierstrasspprime-derivative-of-the-p-function)
    - [WEIERSTRASSZETA: Weierstrass Zeta Function](#weierstrasszeta-weierstrass-zeta-function)
    - [WEIERSTRASSSIGMA: Weierstrass Sigma Function](#weierstrasssigma-weierstrass-sigma-function)
  - [Associate Elliptic Integrals (B, D, J)](#associate-elliptic-integrals-b-d-j)
    - [ELLIPTICBDJ: Incomplete Associate Integrals](#ellipticbdj-incomplete-associate-integrals)
    - [ELLIPTICBD: Complete Associate Integrals](#ellipticbd-complete-associate-integrals)
    - [JACOBIEDJ: Jacobi-Argument Associate Integrals](#jacobiedj-jacobi-argument-associate-integrals)
  - [Carlson Symmetric Elliptic Integrals](#carlson-symmetric-elliptic-integrals)
    - [CARLSONRF: Symmetric First Kind](#carlsonrf-symmetric-first-kind)
    - [CARLSONRD: Symmetric Second Kind](#carlsonrd-symmetric-second-kind)
    - [CARLSONRJ: Symmetric Third Kind](#carlsonrj-symmetric-third-kind)
    - [CARLSONRC: Degenerate Form](#carlsonrc-degenerate-form)
  - [Bulirsch Complete Elliptic Integral](#bulirsch-complete-elliptic-integral)
    - [CEL: Bulirsch Generalised Complete Integral](#cel-bulirsch-generalised-complete-integral)
    - [CEL1, CEL2, CEL3: Wrappers](#cel1-cel2-cel3-wrappers)
  - [Elliptic Related Functions](#elliptic-related-functions)
    - [AGM: Arithmetic Geometric Mean](#agm-arithmetic-geometric-mean)
    - [NOMEQ: The Value of Nome q = q(m)](#nomeq-the-value-of-nome-q--qm)
    - [INVERSENOMEQ: The Value of Nome m = m(q)](#inversenomeq-the-value-of-nome-m--mq)
  - [Interactive Examples](#interactive-examples)
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
u = Integral(1/sqrt(1-m^2*sin(theta)^2), 0, phi);
Sn(u) = sin(phi);
Cn(u) = cos(phi);
Dn(u) = sqrt(1-m^2*sin(phi)^2);
```

*Depends on* `AGM`, `ELLIPKE`.<br>
*Used by* `THETA`.<br>
*See also* `ELLIPKE`.

## ELLIPJI: Jacobi's elliptic functions of the complex argument

`ELLIPJI` evaluates the Jacobi elliptic functions of complex phase `U`.

`[Sni,Cni,Dni] = ELLIPJ(U,M)` returns the values of the Jacobi elliptic functions `SNI`, `CNI` and `DNI` evaluated for corresponding elements of argument `U` and parameter `M`. The arrays `U` and `M` must be of the same size (or either can be scalar).  As currently implemented, `M` is real and limited to `0 <= M <= 1`.

### Example:

```
[phi1,phi2] = meshgrid(-pi:3/20:pi, -pi:3/20:pi);
phi = phi1 + phi2*i;
[Sni,Cni,Dni]= ellipji(phi, 0.99);
```

*Depends on* `AGM`, `ELLIPJ`, `ELLIPKE`<br>
*See also* `ELLIPTIC12`, `ELLIPTIC12I`

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

`Th = THETA(TYPE,V,M)` returns values of theta functions evaluated for corresponding values of argument `V` and parameter `M`. `TYPE` is a type of the theta function, there are four numbered types. The arrays `V` and `M` must be the same size (or either can be scalar). As currently implemented, `M` is limited to `0 <= M <= 1`.

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

*Depends on* `AGM`, `ELLIPJ`, `ELLIPKE`, `JACOBITHETAETA`<br>
*See also* `ELLIPTIC12`, `ELLIPTIC12I`

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

*Depends on* `ELLIPTIC12`, `AGM`<br>
*See also* `ELLIPKE`, `ELLIPJ`, `ELLIPTIC3`, `THETA`.

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

`[F,E] = ELLIPTIC123(m)` — complete Elliptic Integrals of the first and second kind.

`[F,E] = ELLIPTIC123(b,m)` — incomplete Elliptic Integrals of the first and second kind.

`[F,E,PI] = ELLIPTIC123(m,n)` — complete Elliptic Integrals of the first to third kind.

`[F,E,PI] = ELLIPTIC123(b,m,n)` — incomplete Elliptic Integrals of the first to third kind.

The order of the input arguments has been chosen to be consistent with the pre-existing `elliptic12` and `elliptic3` functions.

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
z in (-E(pi/2,m), E(pi/2,m)) x (0,1) - value of the entire parameter space
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

# Weierstrass's Elliptic Functions

Weierstrass's elliptic functions are elliptic functions that take a particularly simple form; they are named for Karl Weierstrass. The four functions in this library parametrize elliptic curves in terms of the half-period roots `(e1, e2, e3)` — the real roots of `4t^3 - g2*t - g3 = 0`, ordered `e1 > e2 > e3`. Parallel and GPU dispatch follow the same `elliptic_config` mechanism as the rest of the library.

## WEIERSTRASSINVARIANTS: Lattice Invariants

`WEIERSTRASSINVARIANTS` converts the half-period roots to the Weierstrass lattice invariants.

`[g2, g3, Delta] = WEIERSTRASSINVARIANTS(e1, e2, e3)` returns the invariants `g2`, `g3`, and the discriminant `Delta = g2^3 - 27*g3^2` (A&S 18.1.3).

```
g2    = -4*(e1*e2 + e1*e3 + e2*e3)
g3    =  4* e1 * e2 * e3
Delta = g2^3 - 27*g3^2
```

### Example:

```matlab
[g2, g3, D] = weierstrassInvariants(1, 0, -1);
% g2 = 4,  g3 = 0,  D = 64
```

*See also* `WEIERSTRASSP`, `WEIERSTRASSPPRIME`.

## WEIERSTRASSP: Weierstrass P-function

`WEIERSTRASSP` evaluates the Weierstrass ℘-function (pe-function).

`P = WEIERSTRASSP(z, e1, e2, e3)` returns the value of `℘(z; e1, e2, e3)` using the connection formula to Jacobi elliptic functions (A&S 18.9.1):

```
℘(z) = e3 + (e1 - e3) / sn^2(z*sqrt(e1-e3), m),   m = (e2-e3)/(e1-e3)
```

The function has a double pole at `z = 0` and is even and doubly periodic with real half-period `omega1 = K(m)/sqrt(e1-e3)`.

### Example:

```matlab
e1 = 1; e2 = 0; e3 = -1;
z  = linspace(0.05, 1.2, 200);
P  = weierstrassP(z, e1, e2, e3);
plot(z, P);
```

*Depends on* `ELLIPJ`.<br>
*See also* `WEIERSTRASSPPRIME`, `WEIERSTRASSZETA`, `WEIERSTRASSSIGMA`.

## WEIERSTRASSPPRIME: Derivative of the P-function

`WEIERSTRASSPPRIME` evaluates the derivative `℘'(z; e1, e2, e3)`.

`dP = WEIERSTRASSPPRIME(z, e1, e2, e3)` uses the chain-rule formula (A&S 18.9.8):

```
℘'(z) = -2*(e1-e3)^(3/2) * cn(u,m)*dn(u,m) / sn^3(u,m),   u = z*sqrt(e1-e3)
```

The function is odd and satisfies the ODE identity `(℘')^2 = 4℘^3 - g2*℘ - g3`.

### Example:

```matlab
e1 = 1; e2 = 0; e3 = -1;
z  = linspace(0.05, 1.25, 200);
dP = weierstrassPPrime(z, e1, e2, e3);
```

*Depends on* `ELLIPJ`.<br>
*See also* `WEIERSTRASSP`.

## WEIERSTRASSZETA: Weierstrass Zeta Function

`WEIERSTRASSZETA` evaluates the Weierstrass zeta function `ζ(z; e1, e2, e3)`.

`Z = WEIERSTRASSZETA(z, e1, e2, e3)` computes `ζ` via numerical integration using a regularised Gauss-Legendre quadrature scheme. The quasi-period `η1 = ζ(ω1)` is computed once; then `ζ(z)` follows from the quasi-periodicity relation.

The zeta function is **not** doubly periodic but satisfies:

```
ζ'(z)        = -℘(z)
ζ(z + 2ω1)  = ζ(z) + 2η1
ζ(-z)        = -ζ(z)   (odd)
```

Accurate for `|z| < ~4*omega1`; returns `Inf` at lattice poles.

### Example:

```matlab
e1 = 1; e2 = 0; e3 = -1;
z  = linspace(0.05, 1.1, 200);
Z  = weierstrassZeta(z, e1, e2, e3);
```

*Depends on* `WEIERSTRASSP`, `ELLIPKE`.<br>
*See also* `WEIERSTRASSSIGMA`.

## WEIERSTRASSSIGMA: Weierstrass Sigma Function

`WEIERSTRASSSIGMA` evaluates the Weierstrass sigma function `σ(z; e1, e2, e3)`.

`S = WEIERSTRASSSIGMA(z, e1, e2, e3)` computes `σ` from its logarithmic derivative. The integrand `ζ(t) - 1/t` is `O(t^3)` near `t = 0`; a Laurent series handles the near-origin part and Gauss-Legendre quadrature handles the rest.

The sigma function is entire, odd, and satisfies `σ'(z)/σ(z) = ζ(z)` and `σ(0) = 0`, `σ'(0) = 1`.

### Example:

```matlab
e1 = 1; e2 = 0; e3 = -1;
z  = linspace(0.05, 1.0, 200);
S  = weierstrassSigma(z, e1, e2, e3);
```

*Depends on* `WEIERSTRASSZETA`, `ELLIPKE`.<br>
*See also* `WEIERSTRASSZETA`, `WEIERSTRASSP`.

# Associate Elliptic Integrals (B, D, J)

The associate integrals `B(φ|m)`, `D(φ|m)`, `J(φ,n|m)` are more fundamental than `F`, `E`, `Π`: the standard integrals decompose as `F = B+D`, `E = B+(1−m)D`, `Π = B+D+n·J`. Computing via `B/D/J` avoids precision loss when `F−E` or `Π−F` are small (e.g. near `m = 1`).

Reference: Fukushima, T. (2015). *Elliptic functions and elliptic integrals for celestial mechanics and dynamical astronomy*. MNRAS.

## ELLIPTICBDJ: Incomplete Associate Integrals

`ELLIPTICBDJ` evaluates the incomplete associate elliptic integrals simultaneously.

`[B, D, J] = ELLIPTICBDJ(PHI, M, N)` returns:

```
B(φ|m)   = ∫₀^φ cos²θ / sqrt(1-m*sin²θ) dθ
D(φ|m)   = ∫₀^φ sin²θ / sqrt(1-m*sin²θ) dθ
J(φ,n|m) = ∫₀^φ sin²θ / ((1-n*sin²θ)*sqrt(1-m*sin²θ)) dθ
```

Connection to standard integrals: `F = B+D`, `E = B+(1-m)*D`, `Π = B+D+n*J`.

```matlab
phi = 0.8;  m = 0.5;  n = 0.3;
[B, D, J] = ellipticBDJ(phi, m, n);
[F, E]    = elliptic12(phi, m);
assert(abs(B+D-F) < 1e-12)          % F = B+D
assert(abs(B+(1-m)*D-E) < 1e-12)    % E = B+(1-m)*D
Pi = elliptic3(phi, m, n);
assert(abs(B+D+n*J-Pi) < 1e-11)     % Pi = B+D+n*J
```

`[B, D] = ELLIPTICBDJ(PHI, M)` computes only B and D.

*Depends on* `CARLSONRF`, `CARLSONRD`, `CARLSONRJ`.

## ELLIPTICBD: Complete Associate Integrals

`ELLIPTICBD` evaluates the complete associate elliptic integrals `B(m)`, `D(m)`, `S(m)`.

`[B, D, S] = ELLIPTICBD(M)` returns:

```
D(m) = (K(m) - E(m)) / m
B(m) = K(m) - D(m)
S(m) = (D(m) - B(m)) / m
```

Identities: `K = B+D`, `E = B+(1-m)*D`.

```matlab
m = 0.7;
[B, D, S] = ellipticBD(m);
[K, E] = ellipke(m);
disp(abs(B+D - K))     % ≈ 0
disp(abs(B+(1-m)*D - E))  % ≈ 0
```

*Depends on* `ELLIPKE`.

## JACOBIEDJ: Jacobi-Argument Associate Integrals

`JACOBIEDJ` evaluates `E_u`, `D_u`, `J_u` in terms of the Jacobi argument `u = F(φ|m)`.

`[Eu, Du, Ju] = JACOBIEDJ(U, M, N)` returns:

```
E_u(u|m)   = u - m * D_u(u|m)
D_u(u|m)   = integral_0^u sn^2(v|m) dv
J_u(u,n|m) = J(am(u|m), n, m)
```

```matlab
m = 0.5;  n = 0.3;
[K, ~] = ellipke(m);  u = 0.6*K;
[Eu, Du, Ju] = jacobiEDJ(u, m, n);
assert(abs(Eu - (u - m*Du)) < 1e-14)   % E_u = u - m*D_u
```

*Depends on* `ELLIPJ`, `ELLIPTICBDJ`.

# Carlson Symmetric Elliptic Integrals

The Carlson symmetric forms are the DLMF §19 reference forms. They connect directly to the Legendre forms and are the basis for all `B/D/J` computations. All four functions accept scalar or array inputs (broadcasting applies).

*Algorithm:* Carlson duplication iteration (DLMF §19.36). Convergence in ~5 iterations for double precision.

## CARLSONRF: Symmetric First Kind

`RF = CARLSONRF(X, Y, Z)` evaluates the symmetric elliptic integral of the first kind:

```
R_F(x,y,z) = (1/2) * integral_0^inf dt / sqrt((t+x)*(t+y)*(t+z))
```

Connection to Legendre: `F(φ|m) = sin(φ) · R_F(cos²φ, 1−m·sin²φ, 1)`.

```matlab
m = 0.7;
[K, ~] = ellipke(m);
RF = carlsonRF(0, 1-m, 1);
assert(abs(RF - K) < 1e-13)   % RF(0, kc^2, 1) = K(m)
```

## CARLSONRD: Symmetric Second Kind

`RD = CARLSONRD(X, Y, Z)` evaluates the symmetric elliptic integral of the second kind (special case of `R_J` with `p = z`):

```
R_D(x,y,z) = (3/2) * integral_0^inf dt / sqrt((t+x)*(t+y)*(t+z)^3)
```

Connection to Legendre: `D(φ|m) = (sin³φ/3) · R_D(cos²φ, 1−m·sin²φ, 1)`.

## CARLSONRJ: Symmetric Third Kind

`RJ = CARLSONRJ(X, Y, Z, P)` evaluates:

```
R_J(x,y,z,p) = (3/2) * integral_0^inf dt / ((t+p)*sqrt((t+x)*(t+y)*(t+z)))
```

Connection: `Π(n,φ|m) = sin(φ)·R_F + (n·sin³φ/3)·R_J(cos²φ, 1−m·sin²φ, 1, 1−n·sin²φ)`.

Special case: `R_J(x,y,z,z) = R_D(x,y,z)`.

## CARLSONRC: Degenerate Form

`RC = CARLSONRC(X, Y)` evaluates the degenerate symmetric form `R_C(x,y) = R_F(x,y,y)` via closed-form arctan/arctanh (DLMF 19.2.17–18):

```
R_C(0, 1/4) = π
R_C(9/4, 2) = ln(2)
R_C(x, x)   = 1/sqrt(x)
```

# Bulirsch Complete Elliptic Integral

## CEL: Bulirsch Generalised Complete Integral

`C = CEL(KC, P, A, B)` evaluates Bulirsch's generalised complete elliptic integral:

```
cel(kc,p,a,b) = integral_0^{pi/2}
    (a*cos^2(phi) + b*sin^2(phi)) / ((cos^2(phi) + p*sin^2(phi)) * sqrt(cos^2(phi) + kc^2*sin^2(phi))) dphi
```

where `kc = sqrt(1-m)` is the complementary modulus.

Special cases:

| Call                    | Result  |
|-------------------------|---------|
| `cel(kc, 1, 1, 1)`      | K(m)    |
| `cel(kc, 1, 1, kc^2)`   | E(m)    |
| `cel(kc, 1, 1, 0)`      | B(m)    |
| `cel(kc, 1, 0, 1)`      | D(m)    |
| `cel(kc, 1-n, 1, 1)`    | Π(n\|m) |

```matlab
m = 0.5;  kc = sqrt(1-m);
[K, E] = ellipke(m);
assert(abs(cel(kc,1,1,1) - K) < 1e-13)       % K(m)
assert(abs(cel(kc,1,1,kc^2) - E) < 1e-13)    % E(m)
```

*Depends on* `ELLIPTICBD`, `CARLSONRJ`.

## CEL1, CEL2, CEL3: Wrappers

Thin wrappers over `cel`:

- `cel1(kc)` = `cel(kc, 1, 1, 1)` = K(m)
- `cel2(kc, a, b)` = `cel(kc, 1, a, b)`
- `cel3(kc, p)` = `cel(kc, p, 1, 1)` = Π(1−p|m)

```matlab
m = linspace(0.1, 0.9, 5);  kc = sqrt(1-m);
[K, ~] = ellipke(m);
assert(max(abs(cel1(kc) - K)) < 1e-13)   % cel1 = K
```

# Elliptic Related Functions

## AGM: Arithmetic Geometric Mean

`AGM` calculates the [Arithmetic Geometric Mean](http://en.wikipedia.org/wiki/Arithmetic-geometric_mean) of `A` and `B` (see [1](#references)).

`[A,B,C,N]= AGM(A0,B0,C0,TOL)` carry out the process of the arithmetic geometric mean, starting with a given positive numbers triple `(A0, B0, C0)` and returns in
`(A, B, C)` the generated sequence. `N` is a number of steps (returns in the value `uint32`).

The general scheme of the procedure:

```
A(i) = 1/2*( A(i-1)+B(i-1) );     A(0) = A0;
B(i) = sqrt( A(i-1)*B(i-1) );     B(0) = B0;
C(i) = 1/2*( A(i-1)+B(i-1) );     C(0) = C0;
```

Stop at the `N`-th step when `A(N) = B(N)`, i.e., when `C(N) = 0`.

*Used by* `ELLIPJ` and `ELLIPTIC12`.<br>*See also* `ELLIPKE`, `ELLIPTIC3`, `THETA`.

## NOMEQ: The Value of Nome `q = q(m)`

`NOMEQ` gives the value of Nome `q = q(m)`.

Nome `Q = nomeq(M,TOL)`, where `0<=M<=1` is the module and `TOL` is the tolerance (optional). Default value for the tolerance is `eps = 2.220e-16`.

*Used by* `ELLIPJ`.<br>
*Depends on* `ELLIPKE`<br>
*See also* `ELLIPTIC12I`, `ELLIPTIC3`, `THETA`.

## INVERSENOMEQ: The Value of Nome `m = m(q)`

`INVERSENOMEQ` gives the value of Nome `m = m(q)`.

`M = inversenomeq(q)`, where `Q` is the Nome of q-series.

**WARNING**. The function `INVERSENOMEQ` does not return correct values of `M` for `Q > 0.6`, because of computer precision limitation. The function `NomeQ(m)` has an essential singularity at `M = 1`, so it cannot be inverted at this point and actually it is very hard to find and inverse in the neighborhood also.

More precisely:

```
nomeq(1) = 1
nomeq(1-eps) = 0.77548641878026
```

### Example:

```
nomeq(inversenomeq([0.001 0.3 0.4 0.5 0.6 0.7 0.8]))
```

*Used by* `ELLIPJ`.<br>
*Depends on* `ELLIPKE`<br>
*See also* `ELLIPTIC12I`, `ELLIPTIC3`, `THETA`.

# Contributors

- [@moiseevigor](http://github.com/moiseevigor/) - maintainer (Igor Moiseev)
- [@drbitboy](https://github.com/drbitboy) (Brian Carcich)
- [@wspr](https://github.com/wspr) (Will Robertson)

# References

  1. [NIST Digital Library of Mathematical Functions](http://dlmf.nist.gov/)
  1. M. Abramowitz and I.A. Stegun, "[Handbook of Mathematical Functions](http://www.nrbook.com/abramowitz_and_stegun/)" Dover Publications", 1965, Ch. 17.1 - 17.6.
  1. D. F. Lawden, "[Elliptic Functions and Applications](http://www.amazon.com/Elliptic-Functions-Applications-Mathematical-Sciences/dp/0387969659)" Springer-Verlag, vol. 80, 1989
  1. S. Zhang, J. Jin, "[Computation of Special Functions](http://jin.ece.uiuc.edu/specfunc.html)" (Wiley, 1996).
  1. B. Carlson, "[Three Improvements in Reduction and Computation of Elliptic Integrals](https://nvlpubs.nist.gov/nistpubs/jres/107/5/j75car.pdf)", J. Res. Natl. Inst. Stand. Technol. 107 (2002) 413-418.
  1. N. H. Abel, "[Studies on Elliptic Functions](http://mathdl.maa.org/mathDL/46/?pa=content&sa=viewDocument&nodeId=1557)", english translation from french by Marcus Emmanuel Barnes.
  1. J. P. Boyd, "Numerical, Perturbative and Chebyshev Inversion of the Incomplete Elliptic Integral of the Second Kind", Applied Mathematics and Computation (January 2012).
  1. T. Fukushima, "[Elliptic functions and elliptic integrals for celestial mechanics and dynamical astronomy](https://doi.org/10.1093/mnras/stv1629)", MNRAS 452 (2015) 442---460.
