
# Elliptic functions for Matlab and Octave 

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/moiseevigor/elliptic/tree/master.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/moiseevigor/elliptic/tree/master) [![DOI](https://zenodo.org/badge/5762/moiseevigor/elliptic.svg)](https://zenodo.org/badge/latestdoi/5762/moiseevigor/elliptic)


The Matlab/Octave implementation of [Elliptic integrals of three types](http://en.wikipedia.org/wiki/Elliptic_integral), [Jacobi's elliptic functions](http://en.wikipedia.org/wiki/Jacobi%27s_elliptic_functions) and [Jacobi theta functions](http://en.wikipedia.org/wiki/Theta_function) of four types.

The main *GOAL* of the project is to provide the natural Matlab scripts *WITHOUT* external library calls like Maple and others. All scripts are developed to accept tensors as arguments and almost all of them have their complex versions. Performance and complete control on the execution are the main features.

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
    - <font color="red"></font> [INVERSELLIPTIC2: INVERSE Incomplete Elliptic Integrals of the Second Kind](#inverselliptic2-inverse-incomplete-elliptic-integrals-of-the-second-kind)

  - [Weierstrass's elliptic functions (in development)](#weierstrasss-elliptic-functions-in-development)
  - [Elliptic Related Functions](#elliptic-related-functions)
    - [AGM: Artihmetic Geometric Mean](#agm-artihmetic-geometric-mean)
    - [NOMEQ: The Value of Nome q = q(m)](#nomeq-the-value-of-nome-q--qm)
    - [INVERSENOMEQ: The Value of Nome m = m(q)](#inversenomeq-the-value-of-nome-m--mq)
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

`[F,E] = ELLIPTIC123(m)` — complete Elliptic Integrals of the first and second kind.  

`[F,E] = ELLIPTIC123(b,m)` — incomplete Elliptic Integrals of the first and second kind.

`[F,E,PI] = ELLIPTIC123(m,n)` — complete Elliptic Integrals of the first to third kind.  

`[F,E,PI] = ELLIPTIC123(b,m,n)` — incomplete Elliptic Integrals of the first to third kind.

The order of the input arguments has been chosen to be consistent with the pre-existing `elliptic12` and `elliptic3` functions.

This function is still under development and its results are not always well-defined or even able to be calculated (especially for the third elliptic integral with `n>1`). Please see the documentation for further details.

## INVERSELLIPTIC2: INVERSE Incomplete Elliptic Integrals of the Second Kind

`INVERSELLIPTIC2` evaluates the value of the INVERSE Incomplete Elliptic Integrals of the Second Kind.

`INVERSELLIPTIC2` uses the method described by Boyd J. P. to determine the value of the inverse Incomplete Elliptic Integrals of the Second Kind using the “Empirical” initialization to the Newton’s iteration method [7](#references). 

Elliptic integral of the second kind:

```
E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
```
 
“Empirical” initialization [7](#references):

```
T0(z,m) = pi/2 + sqrt(r)/(theta − pi/2)
```

where 

```
z \in [E(pi/2,m)](−E(pi/2,m),)x[1](0,) - value of the entire parameter space
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

<font color="red">!IN DEVELOPMENT, help needed!</font>.

Weierstrass's elliptic functions are elliptic functions that take a particularly simple form (cf Jacobi's elliptic functions); they are named for Karl Weierstrass. This class of functions are also referred to as p-functions and generally written using the symbol `℘` (a stylised letter p called Weierstrass p).

The Weierstrass elliptic function can be defined in three closely related ways, each of which possesses certain advantages. One is as a function of a complex variable z and a lattice Λ in the complex plane. Another is in terms of z and two complex numbers ω1 and ω2 defining a pair of generators, or periods, for the lattice. The third is in terms z and of a modulus τ in the upper half-plane. This is related to the previous definition by `τ = ω2 / ω1`, which by the conventional choice on the pair of periods is in the upper half-plane. Using this approach, for fixed z the Weierstrass functions become modular functions of τ.


# Elliptic Related Functions

## AGM: Artihmetic Geometric Mean

`AGM` calculates the [Artihmetic Geometric Mean](http://en.wikipedia.org/wiki/Arithmetic-geometric_mean) of `A` and `B` (see [1](#references)). 

`[A,B,C,N]= AGM(A0,B0,C0,TOL)` carry out the process of the arithmetic geometric mean, starting with a given positive numbers triple `(A0, B0, C0)` and returns in 
`(A, B, C)` the generated sequence. `N` is a number of steps (returns in the value`uint32`).

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

**WARNING**. The function `INVERSENOMEQ` does not return correct values of `M` for `Q > 0.6`, because of computer precision limitation. The function `NomeQ(m)` has an essential singularity at `M = 1`, so it cannot be inverted at this point and actually it is very hard to find and inverse in the neigborhood also.

More preciesly:

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

# Running tests

From root folder run `runtests` or check the CI build status [![CircleCI](https://dl.circleci.com/status-badge/img/gh/moiseevigor/elliptic/tree/master.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/moiseevigor/elliptic/tree/master)

```
runtests("tests")
Processing files in /home/igor/Work/elliptic/tests:

  testEllipj.m ................................................ PASS      5/5   
  testElliptic12.m ............................................ PASS      5/5   
  testElliptic3.m ............................................. PASS      5/5
```

# Contributors

Contributors

- [@moiseevigor](http://github.com/moiseevigor/) - maintainer (Igor Moiseev)
- [@drbitboy](https://github.com/drbitboy) (Brian Carcich)
- [@wspr](https://github.com/wspr) (Will Robertson)


# References

  1. [NIST Digital Library of Mathematical Functions](http://dlmf.nist.gov/)
  1. M. Abramowitz and I.A. Stegun, "[Handbook of Mathematical Functions](http://www.nrbook.com/abramowitz_and_stegun/)" Dover Publications", 1965, Ch. 17.1 - 17.6.
  1. D. F. Lawden, "[Elliptic Functions and Applications](http://www.amazon.com/Elliptic-Functions-Applications-Mathematical-Sciences/dp/0387969659)" Springer-Verlag, vol. 80, 1989
  1. S. Zhang, J. Jin, "[Computation of Special Functions](http://jin.ece.uiuc.edu/specfunc.html)" (Wiley, 1996).
  1. B. Carlson, "[Three Improvements in Reduction and Computation of Elliptic Integrals](http://nvl.nist.gov/pub/nistpubs/jres/107/5/j75car.pdf)", J. Res. Natl. Inst. Stand. Technol. 107 (2002) 413-418.
  1. N. H. Abel, "[Studies on Elliptic Functions](http://mathdl.maa.org/mathDL/46/?pa=content&sa=viewDocument&nodeId=1557)", english translation from french by Marcus Emmanuel Barnes. Original "Recherches sur les fonctions elliptiques", Journal fr die reine und angewandte Mathematik, Vol. 2, 1827. pp. 101-181.
  1. B. C. Berndt, H. H. Chan, S.-S. Huang, "[Incomplete Elliptic Integrals in Ramanujan's Lost Notebook](http://www.math.uiuc.edu/~berndt/articles/ellipticintegral.pdf)" in q-series from a Contemporary Perspective, M. E. H. Ismail and D. Stanton, eds., Amer. Math. Soc., 2000, pp. 79-126.
  1. J. P. Boyd, "Numerical, Perturbative and Chebyshev Inversion of the Incomplete Elliptic Integral of the Second Kind", Applied Mathematics and Computation (January 2012)
