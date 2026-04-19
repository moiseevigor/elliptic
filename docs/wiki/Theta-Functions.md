# Jacobi Theta Functions

The Jacobi theta functions are the elliptic analogs of the exponential function and may be used to express the Jacobi elliptic functions. They are quasi-doubly periodic and play a central role in the theory of elliptic functions, abelian varieties, moduli spaces, and quadratic forms.

---

# Definitions

There are four Jacobi theta functions, defined by their Fourier series in terms of the nome `q`:

```
θ₁(z, q) = 2 Σ_{n=0}^∞ (-1)ⁿ q^{(n+1/2)²} sin((2n+1)z)

θ₂(z, q) = 2 Σ_{n=0}^∞ q^{(n+1/2)²} cos((2n+1)z)

θ₃(z, q) = 1 + 2 Σ_{n=1}^∞ q^{n²} cos(2nz)

θ₄(z, q) = 1 + 2 Σ_{n=1}^∞ (-1)ⁿ q^{n²} cos(2nz)
```

The nome `q` is related to the parameter `m` by:

```
q = exp(-π K(1-m) / K(m))
```

where `K(m)` is the complete elliptic integral of the first kind.

---

# Properties

## Values at z = 0

```
θ₁(0, q) = 0
θ₂(0, q) = 2 q^{1/4} (1 + q² + q⁶ + q¹² + ...)
θ₃(0, q) = 1 + 2(q + q⁴ + q⁹ + q¹⁶ + ...)
θ₄(0, q) = 1 - 2(q - q⁴ + q⁹ - q¹⁶ + ...)
```

## Symmetries

- `θ₁` is odd: `θ₁(-z) = -θ₁(z)`
- `θ₂`, `θ₃`, `θ₄` are even: `θⱼ(-z) = θⱼ(z)` for `j = 2, 3, 4`

## Jacobi's Identity

The theta functions satisfy the famous identity:

```
θ₂(0)⁴ + θ₄(0)⁴ = θ₃(0)⁴
```

This is closely related to the Pythagorean theorem for elliptic functions.

---

# Relation to Jacobi Elliptic Functions

The Jacobi elliptic functions can be expressed as ratios of theta functions:

```
sn(u, m) = (θ₃(0) / θ₂(0)) · (θ₁(v) / θ₄(v))
cn(u, m) = (θ₄(0) / θ₂(0)) · (θ₂(v) / θ₄(v))
dn(u, m) = (θ₄(0) / θ₃(0)) · (θ₃(v) / θ₄(v))
```

where `v = u / θ₃(0)²` (with appropriate normalization).

---

# Theta and Eta Functions

The library also provides the auxiliary Jacobi theta `Θ(u)` and eta `H(u)` functions via `jacobiThetaEta`. These are related to the standard theta functions by:

```
H(u, m) = θ₁(v, q)    (Jacobi eta)
Θ(u, m) = θ₄(v, q)    (Jacobi theta)
```

where the argument `v` is normalized by `2K/π`.

---

# Nome

The nome `q` is a fundamental quantity relating the parameter `m` to the theta function series. The library provides:

```matlab
q = nomeq(m);           % compute nome from parameter
m = inversenomeq(q);    % recover parameter from nome
```

These are useful for converting between different parameterizations of elliptic functions.

---

# Library Usage

## Theta Functions of Four Types

```matlab
setup;  % add src/ to path

% Compute theta functions
Th1 = theta(1, z, m);    % θ₁
Th2 = theta(2, z, m);    % θ₂
Th3 = theta(3, z, m);    % θ₃
Th4 = theta(4, z, m);    % θ₄

% Table generation
[phi, alpha] = meshgrid(0:5:90, 0:2:90);
Th1 = theta(1, pi/180*phi, sin(pi/180*alpha).^2);
Th2 = theta(2, pi/180*phi, sin(pi/180*alpha).^2);
Th3 = theta(3, pi/180*phi, sin(pi/180*alpha).^2);
Th4 = theta(4, pi/180*phi, sin(pi/180*alpha).^2);
```

## Theta Functions and Derivatives

```matlab
% Compute theta function and its derivative simultaneously
[Th, ThPrime] = theta_prime(1, z, m);    % θ₁ and θ₁'
[Th, ThPrime] = theta_prime(2, z, m);    % θ₂ and θ₂'
[Th, ThPrime] = theta_prime(3, z, m);    % θ₃ and θ₃'
[Th, ThPrime] = theta_prime(4, z, m);    % θ₄ and θ₄'
```

## Jacobi Theta and Eta

```matlab
% Compute Jacobi's Θ and H functions
[Th, H] = jacobiThetaEta(u, m);
```

---

# Examples

## Verifying Jacobi's Identity

```matlab
m = 0.7;
z = 0;
T2 = theta(2, z, m);
T3 = theta(3, z, m);
T4 = theta(4, z, m);
% θ₂⁴ + θ₄⁴ = θ₃⁴
assert(abs(T2^4 + T4^4 - T3^4) < 1e-10);
```

## Plotting Theta Functions

```matlab
z = linspace(0, pi, 200);
m = 0.5;
figure;
plot(z, theta(1, z, m), z, theta(2, z, m), ...
     z, theta(3, z, m), z, theta(4, z, m));
legend('\theta_1', '\theta_2', '\theta_3', '\theta_4');
title('Jacobi Theta Functions (m = 0.5)');
```

---

# References

1. M. Abramowitz and I.A. Stegun, "[Handbook of Mathematical Functions](https://personal.math.ubc.ca/~cbm/aands/)" Dover Publications, 1965, Ch. 16.27-16.36.
2. D. F. Lawden, "[Elliptic Functions and Applications](https://www.amazon.com/Elliptic-Functions-Applications-Mathematical-Sciences/dp/0387969659)" Springer-Verlag, vol. 80, 1989.
3. [NIST Digital Library of Mathematical Functions, Ch. 20](https://dlmf.nist.gov/20) — Theta Functions.
4. [Wikipedia: Theta function](https://en.wikipedia.org/wiki/Theta_function)
5. E. T. Whittaker and G. N. Watson, "A Course of Modern Analysis", 4th ed., Cambridge University Press, 1990, Ch. 21.
