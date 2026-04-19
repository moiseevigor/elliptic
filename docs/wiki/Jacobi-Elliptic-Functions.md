# Jacobi Elliptic Functions

The Jacobi elliptic functions are a set of twelve functions that generalize the trigonometric functions to the elliptic case. They arise naturally as solutions to the pendulum equation and other nonlinear differential equations.

---

# Definition

The Jacobi elliptic functions `sn`, `cn`, `dn` are defined through the incomplete elliptic integral of the first kind. Given:

```
u = F(φ | m) = ∫₀^φ dθ / √(1 - m sin²θ)
```

the **amplitude** function is `φ = am(u, m)`, and the three principal Jacobi elliptic functions are:

```
sn(u | m) = sin(am(u, m))
cn(u | m) = cos(am(u, m))
dn(u | m) = √(1 - m sin²(am(u, m)))
```

These satisfy the fundamental identities:

```
sn² + cn² = 1
m·sn² + dn² = 1
```

## Parameter Convention

The parameter `m` (sometimes called the modulus squared) is related to the modular angle `α` and the elliptic modulus `k` by:

```
m = k² = sin²α
```

The complementary parameter is `m₁ = 1 - m`.

---

# Special Cases

| Parameter | `sn(u, m)` | `cn(u, m)` | `dn(u, m)` |
|-----------|------------|------------|------------|
| `m = 0`   | `sin(u)`   | `cos(u)`   | `1`        |
| `m = 1`   | `tanh(u)`  | `sech(u)`  | `sech(u)`  |

When `m = 0`, the elliptic functions reduce to ordinary trigonometric functions. When `m = 1`, they reduce to hyperbolic functions.

---

# Periodicity

The Jacobi elliptic functions are doubly periodic in the complex plane. Their real and imaginary periods are:

| Function | Real period | Imaginary period |
|----------|------------|-----------------|
| `sn`     | `4K`       | `2iK'`          |
| `cn`     | `4K`       | `2K + 2iK'`     |
| `dn`     | `2K`       | `4iK'`          |

where `K = K(m)` is the complete elliptic integral of the first kind and `K' = K(1 - m)`.

---

# The Remaining Nine Functions

Beyond `sn`, `cn`, `dn`, there are nine additional Jacobi elliptic functions defined as ratios:

```
ns = 1/sn    nc = 1/cn    nd = 1/dn
sc = sn/cn   cs = cn/sn   ds = dn/sn
sd = sn/dn   dc = dn/cn   cd = cn/dn
```

The naming convention uses a two-letter code `pq` where the function equals `p(u) / q(u)`, with `s = sn`, `c = cn`, `d = dn`, `n = 1`.

---

# Library Usage

## Real Arguments

```matlab
setup;  % add src/ to path

% Compute sn, cn, dn for real arguments
[Sn, Cn, Dn] = ellipj(u, m);

% Example: generate a table
u = linspace(0, 4, 100);
m = 0.5;
[Sn, Cn, Dn] = ellipj(u, m);
plot(u, Sn, u, Cn, u, Dn);
legend('sn', 'cn', 'dn');
```

## Complex Arguments

```matlab
% Compute Jacobi elliptic functions for complex arguments
[Sni, Cni, Dni] = ellipji(u_complex, m);
```

## Relation to Theta Functions

The Jacobi elliptic functions can be expressed in terms of theta functions:

```matlab
% Using Jacobi theta and eta functions
[Th, H] = jacobiThetaEta(u, m);

% sn(u, m) = H / Th  (up to scaling)
```

---

# Differential Equations

The Jacobi elliptic functions satisfy the following coupled ODEs:

```
d/du sn(u) = cn(u) · dn(u)
d/du cn(u) = -sn(u) · dn(u)
d/du dn(u) = -m · sn(u) · cn(u)
```

These are the equations of motion for a pendulum in terms of elliptic functions.

---

# Examples

## Simple Pendulum Period

The period of a simple pendulum with amplitude `θ₀` is:

```matlab
theta0 = pi/4;                    % amplitude in radians
m = sin(theta0/2)^2;              % parameter
T_ratio = elliptic12(pi/2, m);    % K(m) = F(π/2, m)
% T = 4*sqrt(L/g) * K(m)
```

## Verifying Identities

```matlab
u = 0.7; m = 0.3;
[Sn, Cn, Dn] = ellipj(u, m);

% Check sn^2 + cn^2 = 1
assert(abs(Sn^2 + Cn^2 - 1) < 1e-14);

% Check m*sn^2 + dn^2 = 1
assert(abs(m*Sn^2 + Dn^2 - 1) < 1e-14);
```

---

# References

1. M. Abramowitz and I.A. Stegun, "[Handbook of Mathematical Functions](https://personal.math.ubc.ca/~cbm/aands/)" Dover Publications, 1965, Ch. 16.
2. D. F. Lawden, "[Elliptic Functions and Applications](https://www.amazon.com/Elliptic-Functions-Applications-Mathematical-Sciences/dp/0387969659)" Springer-Verlag, vol. 80, 1989.
3. [NIST Digital Library of Mathematical Functions, Ch. 22](https://dlmf.nist.gov/22) — Jacobi Elliptic Functions.
4. [Wikipedia: Jacobi elliptic functions](https://en.wikipedia.org/wiki/Jacobi_elliptic_functions)
5. N. H. Abel, "[Studies on Elliptic Functions](https://old.maa.org/sites/default/files/images/upload_library/1/abeltranslation.pdf)", English translation from French by Marcus Emmanuel Barnes.
