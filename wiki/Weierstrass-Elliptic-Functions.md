# Weierstrass Elliptic Functions

The Weierstrass elliptic functions provide a canonical way to parametrize elliptic curves and are fundamental objects in the theory of doubly-periodic functions. This library implements four functions — `weierstrassP`, `weierstrassPPrime`, `weierstrassZeta`, and `weierstrassSigma` — plus a utility `weierstrassInvariants`, all following the conventions of Abramowitz & Stegun Chapter 18.

All functions accept real inputs, support **scalar broadcasting**, and respect the `elliptic_config` parallel/GPU settings.

---

# Input Convention

Every Weierstrass function takes the half-period roots `(e1, e2, e3)` directly, with the ordering constraint `e1 > e2 > e3`. These are the three real roots of the cubic `4t³ − g₂t − g₃ = 0`.

To work from invariants `(g2, g3)`, use `weierstrassInvariants` in reverse via the companion `weierstrassFromInvariants`:

```matlab
% from invariants to roots
e = sort(real(roots([4, 0, -g2, -g3])), 'descend');
e1 = e(1);  e2 = e(2);  e3 = e(3);
```

---

# Definition

The **Weierstrass ℘-function** is defined for a lattice `Λ = {2mω₁ + 2nω₃ : m, n ∈ Z}` as:

```
℘(z; Λ) = 1/z² + Σ_{(m,n)≠(0,0)} [1/(z − ω_{m,n})² − 1/ω_{m,n}²]
```

The ℘-function is:
- **Even:** `℘(−z) = ℘(z)`
- **Doubly periodic:** `℘(z + 2ω₁) = ℘(z + 2ω₃) = ℘(z)`
- Has a **double pole** at each lattice point

---

# Differential Equation

The ℘-function satisfies the cubic ODE:

```
(℘')² = 4℘³ − g₂℘ − g₃
```

The discriminant `Δ = g₂³ − 27g₃² ≠ 0` ensures `y² = 4x³ − g₂x − g₃` is non-singular.

---

# Functions

## weierstrassInvariants

Converts `(e1, e2, e3)` to the Weierstrass lattice invariants (A&S 18.1.3).

```matlab
[g2, g3, Delta] = weierstrassInvariants(e1, e2, e3)
```

```
g2    = −4·(e1·e2 + e1·e3 + e2·e3)
g3    =  4·e1·e2·e3
Delta = g2³ − 27·g3²
```

**Example:**

```matlab
[g2, g3, D] = weierstrassInvariants(1, 0, -1);
% g2 = 4,  g3 = 0,  D = 64
```

---

## weierstrassP

Evaluates the Weierstrass ℘-function using the connection formula to Jacobi elliptic functions (A&S 18.9.1):

```matlab
P = weierstrassP(z, e1, e2, e3)
```

**Algorithm:** `℘(z) = e3 + (e1 − e3) / sn²(u, m)` where `u = z·√(e1−e3)` and `m = (e2−e3)/(e1−e3)`.

- Double pole at `z = 0` (returns `Inf`)
- Even function: `℘(−z) = ℘(z)`
- `℘(ω₁) = e1` exactly

**Example:**

```matlab
e1 = 1;  e2 = 0;  e3 = -1;
z  = linspace(0.05, 1.2, 300);
P  = weierstrassP(z, e1, e2, e3);

% Verify ODE identity
[g2, g3, ~] = weierstrassInvariants(e1, e2, e3);
dP = weierstrassPPrime(z, e1, e2, e3);
max(abs(dP.^2 - (4.*P.^3 - g2.*P - g3)))  % ~1e-8
```

---

## weierstrassPPrime

Evaluates the derivative `℘'(z)` via the chain rule (A&S 18.9.8):

```matlab
dP = weierstrassPPrime(z, e1, e2, e3)
```

**Algorithm:** `℘'(z) = −2·(e1−e3)^(3/2) · cn·dn / sn³`

- Triple pole at `z = 0` (returns `Inf`)
- Odd function: `℘'(−z) = −℘'(z)`
- `℘'(ω₁) = 0` exactly (simple zero at the half-period)

**Example:**

```matlab
e1 = 1;  e2 = 0;  e3 = -1;
z0 = 0.7;  h = 1e-7;
fd  = (weierstrassP(z0+h, e1,e2,e3) - weierstrassP(z0-h, e1,e2,e3)) / (2*h);
dP  = weierstrassPPrime(z0, e1, e2, e3);
abs(fd - dP)  % ~1e-6  (finite-difference agreement)
```

---

## weierstrassZeta

Evaluates the Weierstrass zeta function `ζ(z)`:

```matlab
Z = weierstrassZeta(z, e1, e2, e3)
```

**Key properties:**

```
ζ'(z) = −℘(z)
ζ(z + 2ω₁) = ζ(z) + 2η₁,   η₁ = ζ(ω₁)   (quasi-periodicity)
ζ(−z) = −ζ(z)               (odd)
```

**Algorithm (A&S 18.10.1):**

1. Compute half-period `ω₁ = K(m)/√(e1−e3)` and `η₁ = ζ(ω₁)` via regularised GL quadrature:
   `η₁ = 1/ω₁ + ∫₀^{ω₁} [−℘(t) + 1/t²] dt`
   The Laurent correction on `[0, ε₀]` avoids catastrophic cancellation.
2. Period-reduce `z` to `z_red ∈ (−ω₁, ω₁]` using `z_red = z − k·2ω₁`.
3. For `z_red > 0`: integrate `ζ(z_red) = η₁ + ∫_{ω₁}^{z_red} [−℘(t)] dt`.
   For `z_red < 0`: use odd symmetry `ζ(z_red) = −ζ(−z_red)` to avoid integrating through the ℘ pole at `t = 0`.

Returns `Inf` at lattice poles; accurate for `|z| < ~4ω₁`.

**Example:**

```matlab
e1 = 1;  e2 = 0;  e3 = -1;
z  = linspace(0.05, 1.1, 200);
Z  = weierstrassZeta(z, e1, e2, e3);

% Verify derivative identity  ζ'(z) = −℘(z)
h  = 1e-6;  z0 = 0.6;
Zp = (weierstrassZeta(z0+h, e1,e2,e3) - weierstrassZeta(z0-h, e1,e2,e3)) / (2*h);
P0 = weierstrassP(z0, e1, e2, e3);
abs(Zp + P0)   % ~1e-10
```

---

## weierstrassSigma

Evaluates the Weierstrass sigma function `σ(z)`:

```matlab
S = weierstrassSigma(z, e1, e2, e3)
```

**Key properties:**

```
σ'(z)/σ(z) = ζ(z)   (logarithmic derivative)
σ(0)  = 0,  σ'(0) = 1
σ(−z) = −σ(z)        (odd, entire)
```

**Algorithm:**

Integrates `log σ(z) = ∫₀^{|z|} [ζ(t) − 1/t] dt`, split at `t_split = min(0.25·ω₁, 0.99·|z|)`:

- `[0, t_split]`: Laurent series `∫ [ζ(t)−1/t]dt = −g₂·t⁴/240 − g₃·t⁶/840`  
  (uses A&S 18.5.5; avoids cancellation since `ζ(t)−1/t = O(t³)`)
- `[t_split, |z|]`: 10-point Gauss-Legendre on `ζ(t) − 1/t`  
  (`t ≥ t_split ≥ 0.25·ω₁` ensures ζ is accurately computed)

Odd symmetry `σ(−z) = −σ(z)` is applied for `z < 0`. Overflows for `|z| > ~4ω₁`.

**Example:**

```matlab
e1 = 1;  e2 = 0;  e3 = -1;
z  = linspace(0.05, 1.0, 200);
S  = weierstrassSigma(z, e1, e2, e3);

% Verify log-derivative: σ'(z)/σ(z) = ζ(z)
h  = 1e-7;  z0 = 0.6;
S0 = weierstrassSigma(z0, e1, e2, e3);
Sp = (weierstrassSigma(z0+h, e1,e2,e3) - weierstrassSigma(z0-h, e1,e2,e3)) / (2*h);
Z0 = weierstrassZeta(z0, e1, e2, e3);
abs(Sp/S0 - Z0)   % ~1e-10
```

---

# Parallel and GPU Support

All four functions follow the same dispatch pattern as the rest of the library:

```matlab
% Enable multi-core parallel execution
elliptic_config('parallel', true);

% Enable GPU acceleration (requires MATLAB Parallel Computing Toolbox + CUDA)
elliptic_config('gpu', true);
```

GPU acceleration cascades through the call chain: `weierstrassSigma` → `weierstrassZeta` → `weierstrassP` → `ellipj` (GPU path). No changes to calling code are needed.

---

# Connection to Jacobi Functions

The fundamental connection formula (A&S 18.9.1):

```
℘(z; e1, e2, e3) = e3 + (e1 − e3) / sn²(u, m)
```

where `u = z·√(e1−e3)` and `m = (e2−e3)/(e1−e3)`.

---

# Connection to Theta Functions

The ℘-function can also be expressed using theta functions:

```
℘(z) − e₁ = (π θ₃(0) θ₄(0) θ₂(v) / (2ω₁ θ₁(v)))²
```

where `v = πz/(2ω₁)`.

---

# Applications

- **Elliptic curve cryptography:** The group law on `y² = 4x³ − g₂x − g₃` is expressed naturally using ℘
- **Integrable systems:** Many soliton equations (KdV, KP) have solutions in terms of ℘
- **Classical mechanics:** The spinning top, geodesics on an ellipsoid
- **String theory:** Weierstrass functions appear in string amplitudes on tori

---

# References

1. M. Abramowitz and I.A. Stegun, "[Handbook of Mathematical Functions](https://personal.math.ubc.ca/~cbm/aands/)" Dover Publications, 1965, Ch. 18.
2. [NIST Digital Library of Mathematical Functions, Ch. 23](https://dlmf.nist.gov/23) — Weierstrass Elliptic and Modular Functions.
3. D. F. Lawden, "[Elliptic Functions and Applications](https://www.amazon.com/Elliptic-Functions-Applications-Mathematical-Sciences/dp/0387969659)" Springer-Verlag, vol. 80, 1989.
4. E. T. Whittaker and G. N. Watson, "A Course of Modern Analysis", 4th ed., Cambridge University Press, 1990, Ch. 20.
5. [Wikipedia: Weierstrass elliptic function](https://en.wikipedia.org/wiki/Weierstrass_elliptic_function)
