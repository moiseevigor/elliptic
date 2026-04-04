# Elliptic Integrals

If `R(x, y)` is a rational function of `x` and `y`, where `y^2` is equal to a cubic or quadratic polynomial in `x`, the integral

```
f(x) = ∫ R(x, y) dx
```

is called an **elliptic integral**.

Elliptic integrals originally arose in connection with the problem of giving the arc length of an ellipse. They were first studied by Giulio Fagnano and Leonhard Euler.

In general, elliptic integrals cannot be expressed in terms of elementary functions. Exceptions to this general rule are when `P` has repeated roots, or when `R(x, y)` contains no odd powers of `y`. However, with the appropriate reduction formula, every elliptic integral can be brought into a form that involves integrals over rational functions and the three canonical forms (i.e. the elliptic integrals of the first, second and third kind).

## Applications

Elliptic integrals appear in many areas of mathematics and physics:

- Arc-length of plane curves (ellipse, hyperbola, Bernoulli's lemniscate)
- Surface area of ellipsoid in 3-dimensional space
- Electric and magnetic field associated with ellipsoid
- Periodicity of an anharmonic oscillator
- Mutual inductance of coaxial circles
- Age of the universe in the Friedman model
- [Euler's three-body problem](https://en.wikipedia.org/wiki/Euler%27s_three-body_problem)
- [And many others](https://arxiv.org/search/?query=elliptic+integral&searchtype=all)

---

# Elliptic Integral of the First Kind

## Definition

The incomplete elliptic integral of the first kind is defined as:

```
F(φ | m) = ∫₀^φ dθ / √(1 - m sin²θ)
         = ∫₀^x dt / √((1 - t²)(1 - m t²))
```

where `m = sin²α` is the parameter and `α` is the modular angle. The substitution `t = sin θ` transforms between the two forms. When `φ = π/2`, we obtain the **complete** elliptic integral of the first kind `K(m)`.

## Library Usage

```matlab
setup;  % add src/ to path

% Incomplete elliptic integral of the first kind
[F, E, Z] = elliptic12(phi, m);

% Generate Abramowitz & Stegun tables (pp. 613-621)
[phi, alpha] = meshgrid(0:5:90, 0:2:90);
[F, E, Z] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);

% Complete and incomplete integrals with extended interface
[F, E] = elliptic123(m);           % complete
[F, E] = elliptic123(m, phi);      % incomplete
```

---

# Elliptic Integral of the Second Kind

## Definition

The incomplete elliptic integral of the second kind is defined as:

```
E(φ | m) = ∫₀^φ √(1 - m sin²θ) dθ
```

When `φ = π/2`, we obtain the **complete** elliptic integral of the second kind `E(m)`.

The Jacobi zeta function is related to both integrals of the first and second kind:

```
Z(φ | m) = E(φ | m) - [E(m) / K(m)] F(φ | m)
```

## Library Usage

The function `elliptic12` returns all three quantities simultaneously:

```matlab
[F, E, Z] = elliptic12(phi, m);
```

For the inverse problem (given `E`, find `φ`):

```matlab
phi = inverselliptic2(E, m);
```

---

# Elliptic Integral of the Third Kind

## Definition

The incomplete elliptic integral of the third kind is defined as:

```
Π(φ, m, c) = ∫₀^φ dθ / ((1 - c sin²θ) √(1 - m sin²θ))
```

where `c` is an additional parameter. This integral reduces to the first kind when `c = 0`.

## Library Usage

```matlab
Pi = elliptic3(phi, m, c);

% Generate Abramowitz & Stegun tables (pp. 625-626)
[phi, alpha, c] = meshgrid(0:15:90, 0:15:90, 0:0.1:1);
Pi = elliptic3(pi/180*phi, sin(pi/180*alpha).^2, c);
```

---

# Complex Arguments

The library also supports complex-valued arguments for the elliptic integrals of the first and second kind:

```matlab
[Fi, Ei, Zi] = elliptic12i(phi_complex, m);
```

This computes the integrals along a path in the complex plane, which is useful for applications involving analytic continuation.

---

# Numerical Methods

## The Arithmetic-Geometric Mean

The arithmetic-geometric mean (AGM) of two positive real numbers `x` and `y` is computed by iterating:

```
a_{n+1} = (a_n + g_n) / 2     (arithmetic mean)
g_{n+1} = √(a_n * g_n)        (geometric mean)
```

starting from `a_0 = x`, `g_0 = y`. Both sequences converge quadratically to the same limit `M(x, y)`.

The complete elliptic integral of the first kind is related to the AGM by:

```
K(m) = π / (2 * M(1, √(1 - m)))
```

The library provides the AGM function directly:

```matlab
[a, b, c, n] = agm(a0, b0, c0);
```

## Landen's Transformation

[Landen's transformation](https://en.wikipedia.org/wiki/Landen%27s_transformation), independently rediscovered by Gauss, is a mapping of the parameters of an elliptic integral that leaves the value of the integral unchanged. It provides the basis for the descending transformation used in `elliptic12`.

## Carlson's Method

The conventional methods for computing elliptic integrals are Gauss and Landen transformations, which converge quadratically and work well for elliptic integrals of the first and second kinds. Unfortunately they suffer from loss of significant digits for the third kind. Carlson's algorithm provides a unified method for all three kinds with satisfactory precision. The third kind integral in this library uses a Gauss-Legendre 10-point quadrature instead.

---

# References

1. M. Abramowitz and I.A. Stegun, "[Handbook of Mathematical Functions](https://personal.math.ubc.ca/~cbm/aands/)" Dover Publications, 1965, Ch. 17.1 - 17.6.
2. D. F. Lawden, "[Elliptic Functions and Applications](https://www.amazon.com/Elliptic-Functions-Applications-Mathematical-Sciences/dp/0387969659)" Springer-Verlag, vol. 80, 1989.
3. S. Zhang, J. Jin, "[Computation of Special Functions](https://www.amazon.com/Computation-Special-Functions-Shanjie-Zhang/dp/0471119636)" (Wiley, 1996).
4. B. Carlson, "[Three Improvements in Reduction and Computation of Elliptic Integrals](https://pmc.ncbi.nlm.nih.gov/articles/PMC4861378/)", J. Res. Natl. Inst. Stand. Technol. 107 (2002) 413-418.
5. [NIST Digital Library of Mathematical Functions, Ch. 19](https://dlmf.nist.gov/19)

# Wiki References

1. [Landen's transformation](https://en.wikipedia.org/wiki/Landen%27s_transformation) — used in elliptic integral numerical evaluations.
2. [Elliptic integral](https://en.wikipedia.org/wiki/Elliptic_integral) — overview of theoretical aspects.
3. [Arithmetic-geometric mean](https://en.wikipedia.org/wiki/Arithmetic%E2%80%93geometric_mean) — foundation for fast computation.

# External Links

- [ALGLIB](https://www.alglib.net/) — Cross-platform numerical analysis and data processing library.
- [Boost C++ Elliptic Integrals](https://www.boost.org/doc/libs/release/libs/math/doc/html/math_toolkit/ellint.html) — Elliptic integral overview in the Boost Math library.
- [OEIS - Integer Sequences related to Elliptic Integrals](https://oeis.org/) — The On-Line Encyclopedia of Integer Sequences.
