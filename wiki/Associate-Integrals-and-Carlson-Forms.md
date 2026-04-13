# Associate Elliptic Integrals, Carlson Forms, and Bulirsch CEL

This page documents the functions added in the `weierstrass` branch that implement:

1. **Carlson symmetric forms** — R_F, R_D, R_J, R_C (DLMF §19)
2. **Incomplete associate integrals** — B(φ|m), D(φ|m), J(φ,n|m) via Carlson
3. **Complete associate integrals** — B(m), D(m), S(m)
4. **Jacobi-argument forms** — E_u(u|m), D_u(u|m), J_u(u,n|m)
5. **Bulirsch's complete integral** — cel(kc, p, a, b)

---

# Why B, D, J?

The standard integrals F, E, Π are related to B, D, J by:

```
F(φ|m)   = B(φ|m) + D(φ|m)
E(φ|m)   = B(φ|m) + (1−m)·D(φ|m)
Π(n,φ|m) = B(φ|m) + D(φ|m) + n·J(φ,n|m)
```

When F−E is small (m near 0 or φ near 0), computing B and D separately avoids catastrophic cancellation. Similarly, J isolates the third-kind contribution without the pole structure of Π−F.

---

# Carlson Symmetric Forms

## carlsonRC

Computes R_C(x,y) = R_F(x,y,y) via closed-form (DLMF §19.2.17–18):

```
y > x ≥ 0 : R_C = arctan(√((y−x)/x))  / √(y−x)
0 < y < x  : R_C = arctanh(√((x−y)/x)) / √(x−y)
y = x      : R_C = 1/√x
x = 0      : R_C = π/(2√y)
```

```matlab
carlsonRC(0, 0.25)   % = pi
carlsonRC(1, 1)      % = 1.0
```

## carlsonRF

Computes R_F(x,y,z) via Carlson duplication (DLMF §19.36.1). Converges in ~5 iterations for double precision.

```matlab
phi = 0.8;  m = 0.5;
s = sin(phi);  c = cos(phi);  d = sqrt(1-m*s^2);
F = s * carlsonRF(c^2, d^2, 1);   % = F(phi|m)
```

## carlsonRD

Computes R_D(x,y,z) = R_J(x,y,z,z) via Carlson duplication (DLMF §19.36.2).

```matlab
D = (sin(phi)^3 / 3) * carlsonRD(c^2, d^2, 1);   % = D(phi|m)
```

## carlsonRJ

Computes R_J(x,y,z,p) via Carlson duplication (DLMF §19.36.3). At each step requires one R_C sub-call.

```matlab
n = 0.3;
Pi = s*carlsonRF(c^2,d^2,1) + (n*s^3/3)*carlsonRJ(c^2,d^2,1,1-n*s^2);  % = Π(n,phi|m)
```

---

# Associate Incomplete Integrals

## ellipticBDJ

```matlab
[B, D, J] = ellipticBDJ(phi, m, n)   % B, D, J simultaneously
[B, D]    = ellipticBDJ(phi, m)       % only B and D
```

Internally computes:
- `F = s * RF(c², Δ², 1)`
- `D = (s³/3) * RD(c², Δ², 1)`
- `B = F − D`
- `J = (s³/3) * RJ(c², Δ², 1, 1−n·s²)`  [only when n requested]

**Example: verify ODE identities:**
```matlab
phi = 0.8;  m = 0.5;  n = 0.3;
[B, D, J] = ellipticBDJ(phi, m, n);
[F, E]    = elliptic12(phi, m);
Pi        = elliptic3(phi, m, n);
abs(B+D-F)       % ~1e-16
abs(B+(1-m)*D-E) % ~1e-16
abs(B+D+n*J-Pi)  % ~1e-11
```

## ellipticBD

```matlab
[B, D, S] = ellipticBD(m)
```

Complete integrals at φ=π/2. Algorithm: `D = (K−E)/m`, `B = K−D`, `S = (D−B)/m`.

```matlab
m = 0.7;
[B, D, S] = ellipticBD(m);
[K, E] = ellipke(m);
abs(B+D-K)       % ~0 (exact)
abs(B+(1-m)*D-E) % ~0 (exact)
```

---

# Jacobi-Argument Forms

## jacobiEDJ

Converts Jacobi argument u to amplitude φ = am(u|m) via `ellipj`, then calls `ellipticBDJ`.

```matlab
[Eu, Du, Ju] = jacobiEDJ(u, m, n)
```

Key identity: `E_u(u|m) = u − m·D_u(u|m)`.

```matlab
m = 0.5;  n = 0.3;
[K, ~] = ellipke(m);  u = 0.6*K;
[Eu, Du, Ju] = jacobiEDJ(u, m, n);
abs(Eu - (u - m*Du))   % ~1e-16
```

---

# Bulirsch's Complete Integral

## cel, cel1, cel2, cel3

```matlab
C = cel(kc, p, a, b)   % Generalised: cel(kc,p,a,b)
C = cel1(kc)           % K(1-kc^2)
C = cel2(kc, a, b)     % cel(kc, 1, a, b)
C = cel3(kc, p)        % cel(kc, p, 1, 1) = Pi(1-p|m)
```

The function evaluates:
```
cel(kc,p,a,b) = ∫₀^{π/2} (a·cos²φ + b·sin²φ) / ((cos²φ + p·sin²φ)·Δ) dφ
```
where Δ = √(cos²φ + kc²·sin²φ).

**Conversion table** (kc = √(1−m)):

| Legendre | Bulirsch |
|---|---|
| K(m) | `cel1(kc)` |
| E(m) | `cel2(kc, 1, kc²)` |
| B(m) | `cel2(kc, 1, 0)` |
| D(m) | `cel2(kc, 0, 1)` |
| Π(n\|m) | `cel3(kc, 1−n)` |

**Algorithm**: For p=1, uses `a·B(m) + b·D(m)`. For p≠1:
```
cel = a·K + (b − a·p)·(Π(1−p|m) − K) / (1−p)
```
where Π(1−p|m) = K + (1−p)·J_complete(1−p|m).

---

# Parallel and GPU Support

All functions follow the standard library dispatch pattern:

```matlab
elliptic_config('parallel', true);   % multi-core via parcellfun/parfor
elliptic_config('gpu', true);        % GPU via gpuArray
```

`ellipticBD` and `ellipticBDJ` have full parallel/GPU dispatch.  
The Carlson functions (`carlsonRF`, `carlsonRD`, `carlsonRJ`, `carlsonRC`) and `cel` are called element-wise from the dispatch layer and do not dispatch independently.

---

# References

1. NIST DLMF §19 — Elliptic Integrals, https://dlmf.nist.gov/19
2. B.C. Carlson, "Numerical Computation of Real or Complex Elliptic Integrals," *Numer. Algorithms* 10 (1995), 13–26.
3. T. Fukushima, "Elliptic functions and elliptic integrals for celestial mechanics and dynamical astronomy," (2015).
4. R. Bulirsch, "Numerical computation of elliptic integrals and elliptic functions," *Numer. Math.* 7 (1965), 78–90.
5. M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions," Dover, 1965, Ch. 17.
