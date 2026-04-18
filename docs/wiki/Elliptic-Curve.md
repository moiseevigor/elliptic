# Elliptic Curves

An elliptic curve is a smooth, projective algebraic curve of genus one with a specified point. Elliptic curves are central objects in number theory, algebraic geometry, and cryptography, and are deeply connected to elliptic integrals and elliptic functions.

---

# Definition

An elliptic curve over a field `K` is defined by the **Weierstrass equation**:

```
y² = x³ + ax + b
```

where `a, b ∈ K` and the discriminant `Δ = -16(4a³ + 27b²) ≠ 0` (ensuring the curve is non-singular).

The set of points `(x, y)` satisfying this equation, together with a "point at infinity" `O`, forms an abelian group under a geometric addition law.

---

# Connection to Elliptic Integrals

The arc length of an ellipse is given by an elliptic integral of the second kind:

```
L = 4a ∫₀^{π/2} √(1 - e² sin²θ) dθ = 4a · E(e²)
```

where `a` is the semi-major axis and `e` is the eccentricity. This library computes such integrals:

```matlab
setup;
a = 5; b = 3;
e2 = 1 - (b/a)^2;
[~, E] = elliptic12(pi/2, e2);
L = 4 * a * E;     % perimeter of the ellipse
```

More generally, the functions `arclength_ellipse` computes the arc length directly:

```matlab
% Arc length of an ellipse with semi-axes a, b
L = arclength_ellipse(a, b);
```

---

# Group Law

Points on an elliptic curve form a group. Given two points `P` and `Q`:

1. Draw the line through `P` and `Q` (or the tangent line if `P = Q`)
2. This line intersects the curve at a third point `R`
3. Reflect `R` across the x-axis to get `P + Q`

The point at infinity `O` serves as the identity element.

---

# Uniformization

Over the complex numbers, every elliptic curve is isomorphic to a complex torus `C / Λ`, where `Λ` is a lattice. The isomorphism is given by the Weierstrass elliptic function:

```
φ: C/Λ → E
z ↦ (℘(z), ℘'(z))
```

where `℘(z)` is the [Weierstrass elliptic function](Weierstrass-Elliptic-Functions). This connects elliptic curves to the doubly-periodic functions in this library.

The Jacobi elliptic functions provide an alternative uniformization via the parameter `m`:

```matlab
% Points on the curve y² = (1-x²)(1-m·x²)
u = linspace(0, 4, 100);
[Sn, Cn, Dn] = ellipj(u, 0.5);
% (Sn, Cn·Dn) parametrizes the curve
```

---

# Applications

## Cryptography

Elliptic curve cryptography (ECC) uses the group structure of elliptic curves over finite fields. The difficulty of the elliptic curve discrete logarithm problem provides security equivalent to much larger RSA keys.

## Number Theory

Elliptic curves are central to:
- The proof of Fermat's Last Theorem (Wiles, 1995)
- The Birch and Swinnerton-Dyer conjecture (one of the Millennium Prize Problems)
- Modular forms and the Langlands program

---

# References

1. J. H. Silverman, "[The Arithmetic of Elliptic Curves](https://www.springer.com/gp/book/9780387094939)", Springer, 2nd ed., 2009.
2. [NIST Digital Library of Mathematical Functions, Ch. 23](https://dlmf.nist.gov/23) — Weierstrass Elliptic and Modular Functions.
3. [Wikipedia: Elliptic curve](https://en.wikipedia.org/wiki/Elliptic_curve)
4. M. Abramowitz and I.A. Stegun, "[Handbook of Mathematical Functions](https://personal.math.ubc.ca/~cbm/aands/)" Dover Publications, 1965.
