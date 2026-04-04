# Weierstrass Elliptic Functions

> **Note:** Weierstrass elliptic function support in this library is currently in development.

The Weierstrass elliptic functions provide a canonical way to parametrize elliptic curves and are fundamental objects in the theory of doubly-periodic functions.

---

# Definition

The **Weierstrass ℘-function** (pe-function) is defined for a lattice `Λ = {2mω₁ + 2nω₃ : m, n ∈ Z}` as:

```
℘(z; Λ) = 1/z² + Σ_{(m,n)≠(0,0)} [1/(z - ω_{m,n})² - 1/ω_{m,n}²]
```

where `ω_{m,n} = 2mω₁ + 2nω₃` ranges over all nonzero lattice points.

The ℘-function is:
- **Even:** `℘(-z) = ℘(z)`
- **Doubly periodic:** `℘(z + 2ω₁) = ℘(z + 2ω₃) = ℘(z)`
- Has a **double pole** at each lattice point

---

# Differential Equation

The ℘-function satisfies the cubic differential equation:

```
(℘')² = 4℘³ - g₂℘ - g₃
```

where the **invariants** `g₂` and `g₃` are determined by the lattice:

```
g₂ = 60 Σ' 1/ω_{m,n}⁴
g₃ = 140 Σ' 1/ω_{m,n}⁶
```

The discriminant `Δ = g₂³ - 27g₃² ≠ 0` ensures the curve `y² = 4x³ - g₂x - g₃` is non-singular.

---

# Related Functions

## Weierstrass Zeta Function

The Weierstrass zeta function `ζ(z)` is defined by:

```
ζ'(z) = -℘(z)
```

It is **not** doubly periodic, but satisfies:

```
ζ(z + 2ωⱼ) = ζ(z) + 2ηⱼ
```

where `ηⱼ = ζ(ωⱼ)` are the quasi-periods.

## Weierstrass Sigma Function

The Weierstrass sigma function `σ(z)` is defined by:

```
σ'(z)/σ(z) = ζ(z)
```

It is an entire, odd function with simple zeros at lattice points.

---

# Connection to Jacobi Functions

The Weierstrass and Jacobi elliptic functions are related by:

```
℘(z) = e₃ + (e₁ - e₃) / sn²(u, m)
```

where `u = z√(e₁ - e₃)`, `m = (e₂ - e₃)/(e₁ - e₃)`, and `e₁, e₂, e₃` are the roots of `4t³ - g₂t - g₃ = 0`.

Using this library's Jacobi elliptic functions:

```matlab
setup;

% Given lattice invariants, compute e1, e2, e3
% e1, e2, e3 are roots of 4t³ - g2*t - g3
e = roots([4 0 -g2 -g3]);
e = sort(real(e), 'descend');  % e1 > e2 > e3

% Compute ℘ via Jacobi functions
m = (e(2) - e(3)) / (e(1) - e(3));
scale = sqrt(e(1) - e(3));
[Sn, ~, ~] = ellipj(z * scale, m);
wp = e(3) + (e(1) - e(3)) ./ Sn.^2;
```

---

# Connection to Theta Functions

The ℘-function can also be expressed using theta functions:

```
℘(z) - e₁ = (π θ₃(0) θ₄(0) θ₂(v) / (2ω₁ θ₁(v)))²
```

where `v = πz/(2ω₁)`. The library's theta functions can be used for these computations:

```matlab
Th1 = theta(1, v, m);
Th2 = theta(2, v, m);
```

---

# Applications

- **Elliptic curve cryptography:** The group law on `y² = 4x³ - g₂x - g₃` is expressed naturally using ℘
- **Integrable systems:** Many soliton equations (KdV, KP) have solutions expressed in terms of ℘
- **String theory:** Weierstrass functions appear in the computation of string amplitudes on tori
- **Classical mechanics:** The spinning top, geodesics on an ellipsoid

---

# References

1. M. Abramowitz and I.A. Stegun, "[Handbook of Mathematical Functions](https://personal.math.ubc.ca/~cbm/aands/)" Dover Publications, 1965, Ch. 18.
2. D. F. Lawden, "[Elliptic Functions and Applications](https://www.amazon.com/Elliptic-Functions-Applications-Mathematical-Sciences/dp/0387969659)" Springer-Verlag, vol. 80, 1989.
3. [NIST Digital Library of Mathematical Functions, Ch. 23](https://dlmf.nist.gov/23) — Weierstrass Elliptic and Modular Functions.
4. [Wikipedia: Weierstrass elliptic function](https://en.wikipedia.org/wiki/Weierstrass_elliptic_function)
5. E. T. Whittaker and G. N. Watson, "A Course of Modern Analysis", 4th ed., Cambridge University Press, 1990, Ch. 20.
