# Post-0d09740 Regression and Data-Integrity Audit

**Date:** 2026-07-28  
**Baseline:** `0d09740` on `master`  
**Trigger:** [issue #35](https://github.com/moiseevigor/elliptic/issues/35) and its initial [patch](https://github.com/moiseevigor/elliptic/commit/0d09740)  
**Scope:** MATLAB/Octave and Python implementations, numerical special values,
backend dispatch, data-shape preservation, packaging, and CI coverage.

## Assessment

The initial patch added useful special-value tests, but it did not fully close
the regression class. In particular, the issue #35 phase-reduction fix was
missing from the MATLAB GPU implementation, and a private stale copy inside
`elliptic123.m` could still bypass the repaired public implementation. The
audit also found independent numerical cancellation, pole-handling,
backend-dispatch, and CI-coverage defects.

All fixes below are kept atomic at the function boundary and covered by
regression tests. No public function was intentionally removed.

## Confirmed findings and fixes

| Severity | Area | Failure and data impact | Fix |
|---|---|---|---|
| Critical | MATLAB GPU `elliptic12` | The GPU path did not apply the quasi-period correction from issue #35, so phases outside the principal interval could return the wrong branch. | Apply the same exact period reduction and complete-integral correction as the serial path; add GPU regression coverage. |
| High | `elliptic123` dispatch | Stale private `elliptic12i`/`elliptic12ic` copies shadowed the repaired public implementation and preserved old behavior. | Retire the private names so `elliptic123` routes through the maintained public function; add routing parity tests. |
| High | Python runtime dependencies | `theta`, nome, inverse, and complex functions imported SciPy although SciPy was not declared as a runtime dependency. A base installation failed only when those public functions were called. | Replace runtime SciPy calls with backend-native theta series, Carlson/AGM forms, and fixed-iteration inverse solvers; add a SciPy-free install job. |
| High | `elliptic12`, `m=1` | Reducing the phase before checking the singularity hid crossed first-kind poles (`F(pi,1)` became zero) and under-counted `E` after the first quadrant. | Detect pole crossings from the original phase and restore the full quasi-period contribution. |
| High | `ellipticBD`, small `m` | `(K-E)/m` and `(D-B)/m` catastrophically cancelled. At `m=1e-20`, valid finite limits became `D=0` and an order-`1e20` wrong `S`. | Compute `D` directly with Carlson `RD`, derive `B` from `RF-RD/3`, and use a convergent series for `S` near zero in both languages. |
| High | `elliptic3`, near `n=1` | Fixed quadrature lost roughly seven decimal digits near the third-kind endpoint pole. Negative phases crossing an interior pole were not consistently rejected in Python. | Use Carlson `RF/RJ` throughout Python; use a MATLAB/Octave hybrid that keeps quadrature on regular inputs and switches near poles; validate pole crossings using `abs(phi)`. |
| High | `ellipj`, extreme inputs | Clipping valid parameters changed representable data near `m=0` and `m=1`; large phases entered the amplified Landen recursion directly and could produce order-one errors. | Preserve every interior parameter exactly, reduce by `2K` before descent, reconstruct quasi-period signs, and use stable `dn` and `sech` formulas. |
| High | Weierstrass functions | A magnitude heuristic labelled finite near-pole values as infinity (`P(1e-11)` became `Inf`), complex inputs silently lost their imaginary part, and zeta/sigma forced NumPy arrays. | Detect actual lattice points using reduced periods and ULP-scale tolerance, explicitly reject unsupported complex inputs, and keep zeta/sigma backend-native. |
| High | JAX/PyTorch dispatch | Python control flow, NumPy casts, and scalar flattening broke tracing or moved public results off device in complex, inverse, theta, Jacobi EDJ, Weierstrass, and ellipse helpers. | Replace value-dependent Python branches with masked array expressions and fixed iteration counts; preserve the originating array namespace and shape. |
| Medium | Jacobi theta endpoints | `theta(1,v,0)` returned `sin(v)` instead of zero; tiny valid `m` values were collapsed to the endpoint; representable values near `m=1` returned NaN. | Evaluate native theta series with exact endpoint overrides and no arbitrary endpoint clipping. |
| Medium | MATLAB grouping | `uniquetol(...,1e-11)` treated distinct `m` values as identical and silently substituted one result for another. A `5e-12` parameter difference produced a `1.69e-9` vector/scalar discrepancy. | Group exact duplicates only in `elliptic12` and `ellipj`; add a no-substitution regression. |
| Medium | MATLAB complex integral | `elliptic12i` divided by `m` at the exact `m=0` endpoint. | Use a safe internal denominator and restore the analytic `m=0` values exactly. |
| Medium | Carlson `RJ` tracing | Public eager validation converted JAX tracers to NumPy and failed before evaluation; an unselected `RC` branch could still emit invalid arithmetic. | Restrict exception-producing validation to eager NumPy inputs and make internal masked arguments safe. |
| Medium | Ellipse helper | Public inputs were cast to Python `float`, rejecting arrays, Torch tensors, and JAX tracers and losing broadcast shape. | Implement elementwise backend-native broadcasting, explicit NumPy domain errors, and masked invalid values for traced execution. |
| Medium | CI coverage | Octave CI ran only 6 of 15 test files. Torch/JAX jobs installed their backends but no tests passed those arrays. The release path had no base-dependency-only check. | Discover every `test*.m`, add real NumPy/Torch/JAX dispatch tests, add a SciPy-free standalone job, and make publishing depend on all of them. |

## Numerical and differential checks

- Real `F` and `E` were compared over broad random phases, including many
  quasi-periods, against independent reference implementations; observed
  absolute discrepancies were about `1.3e-12` or smaller.
- Real third-kind values were compared to independent Carlson `RF/RJ`
  references; observed relative discrepancies were about `2.5e-15`.
- Two hundred random complex `F`, `E`, and Jacobi `sn` samples were compared
  with 50-digit `mpmath` values; the largest observed discrepancy was about
  `5e-14`.
- The Weierstrass differential equation residual was checked at random regular
  points; the largest observed relative residual was about `4.5e-15`.
- The large-phase `ellipj` regression at
  `u=1,000,000.123`, `m=nextafter(1,0)` was checked against a 70-digit
  `mpmath` reference in both languages.

## Automated verification

- Python default environment: **451 passed, 1 skipped**. The skip is the
  optional JAX-only test when JAX is not installed.
- Python with both JAX and PyTorch installed: **458 passed**.
- Dedicated backend matrix: **10 passed**, covering NumPy, PyTorch, and JAX;
  JAX compiles the core paths and traces the larger fixed-iteration graphs.
- Octave: **220/220 test blocks passed** across all 16 `test*.m` files.
- Python byte-compilation and whitespace/error checks passed.
- A base-dependency-only public API exercise verifies that SciPy is absent and
  not imported at runtime.

## Adversarial review round (2026-08-16)

An independent adversarial review (Codex, cross-checked against mpmath at
40-100 digits) of the four commits above found five material counterexamples
that the original round missed.  All are fixed and regression-tested
(`testEdgeCases.m` block Q, `test_edge_cases.py::TestAdversarialRound`):

| Finding | Failure | Fix |
|---|---|---|
| `elliptic3` negative amplitude with pole in the complete integral | `0*Inf = NaN` for e.g. `Pi(-1\|.5,1)` (MATLAB) | complete-integral correction applied only where a half-period or reflection is present |
| complex `F/E` small `m` (both ports) | A&S 17.4.11 decomposition loses `sqrt(eps/m)` digits: `F(0.2i\|1e-20)` returned 0, error 9.2e-3 at `m=1e-14` | Maclaurin series through `m^2`, switched at `m*max(1, e^(2\|psi\|)) < 1e-4`; crossover error ~2e-12 vs mpmath |
| Weierstrass pole classification (both ports) | tolerance windows (`abs(sn) < eps^(1/3)` MATLAB; `8*eps*max(1,\|z\|)` Python) returned `Inf` for the finite `P(1e-16) = 1e32` | pole only at the exact lattice point (`z_reduced == 0` / `sn == 0`), per the DLMF 23.9.2 Laurent expansion |
| inverse nome small `q` (both ports) | Python bisection had a `2^-64` absolute floor (`m(1e-30)` off by 9 orders); MATLAB tables documented-unreliable outside `[1e-5, 0.76]` | DLMF 20.9.1 closed form `m = (theta2/theta3)^4` with the `q^(1/4)` factor kept outside the ratio; exact at every scale |
| Carlson `RC` scale invariance (Python) | absolute branch tolerance sent small-scale inputs down the degenerate branch: `RC(1e-20,2e-20)` off 27%, contaminating `RJ` | relative branch selection; homogeneity `RF,RC ~ lambda^(-1/2)`, `RD,RJ ~ lambda^(-3/2)` now tested at `lambda = 1e+/-20` |

Also from that round: reversed arc intervals are now signed for circles and
ellipses alike; DLMF citations corrected (19.25.14 for incomplete third kind,
19.25.5/19.25.9 for F/E); nondegenerate `ellipticBD` anchors at
`m = 0.2, 0.7, 0.999`.

## Adversarial review round 2 (self-review, 2026-08-16)

A second, independent adversarial pass fuzzed every public function of both
ports against mpmath at 40 digits over parameter endpoints (`m -> 0`,
`m -> 1`, `m = 1`), extreme argument scales, near-pole and near-lattice
arguments, exact period multiples +/- ulps, and complex arguments across the
branch point.  Each candidate was classified by evaluating scipy at the same
double inputs: where scipy reaches machine precision the loss is ours.
Nine implementation defects resulted, all fixed and pinned
(`testEdgeCases.m` block R, `test_edge_cases.py::TestAdversarialRound2`):

| Defect | Failure | Fix |
|---|---|---|
| `1 - m sin^2` / `1 - n sin^2` formed by subtraction (`elliptic12` py, `elliptic3` both) | `F(pi/2-1e-9 \| 1-eps/2)` off 4e-3; `Pi(pi/2-1e-6 \| m, n=1)` off 3e-5 | form `(1-m) + m cos^2`, `(1-n) + n cos^2` |
| `F(phi \| 1) = log(tan(pi/4+phi/2))` (both) | `F(0\|1) = -1.1e-16`, wrong sign at `1e-16` | `atanh(sin phi)` |
| `(ratio-1)/m` in the A&S 17.4.11 decomposition (both) | `Im F(pi/2 + 1e-9 i)` returned 0; `sqrt(eps/m)` loss for small m | closed cancellation-free `tan^2(mu) = 2 sinh^2 csc^2 / (B' + sqrt(B'^2 + 4C'))`, m cancels analytically |
| Landen back-substitution `asin(c sin/a)` near +/-1 (both) | `cn(9.4 \| 1-eps/2)` off 5e-10 | `atan2(c sin, sqrt(a^2 cos^2 + b^2 sin^2))` via `a^2 - c^2 = b^2` |
| `R_C` arctanh branch for `y << x` and `log(1+tiny)` (both) | `RC(3,1e-10)` off 3e-9; `RC(1+1e-13,1)` off 3e-10 -- contaminated `R_J` and `J` | `log1p(((x-y)/(sqrt x + sqrt y) + sqrt(x-y))/sqrt y)/sqrt(x-y)` |
| `R_J` duplication capped at 30 (both) | `RJ(1e-20,2e-20,3e-20,.5)` 11% off | 60 fixed (py, ratio limit ~3e32) / adaptive to 200 (MATLAB) |
| two zero Carlson arguments (both) | `RF(0,0,1) = 2e6` | `Inf` (DLMF 19.16) |
| inverse `E`: fold of tiny negative z; tol-gated Newton (both) | rel 1e-7 at `z = -1e-9 E1` | oddness first; unconditional / relative-stop Newton |
| `1 - m` from lattice roots by subtraction; MATLAB `nomeq` via `ellipke(1-m)`; MATLAB `weierstrassP` reducing inside `ellipj` | `q(1e-16)` 11% off, `q(1e-17) = 0`; `P(2 omega1 + 1e-9)` off 40% on near-m=1 lattices | `1-m = (e1-e2)/(e1-e3)`; `K' = R_F(0, m, 1)`; reduce by `2 omega1` before `ellipj` |

Reference-construction lesson recorded for future rounds: anchors must be
evaluated at the *exact double* the library receives (`mpf(float(x))`), not
at the decimal the test author typed -- near singularities the two differ at
the 1e-9 level (`F(pi/2 - 1e-9 \| 1-eps/2)`: 19.6599302656 vs 19.6599302792).

## Adversarial review round 3 (dense fuzz + API abuse, 2026-08-16)

Dense random fuzz (600 points/function/seed, m log-uniform to within 1e-14 of
both endpoints, phases over +/-13 periods and at odd multiples of pi/2 with
1e-9 jitter) against mpmath, plus a logic/API abuse probe (shapes, dtypes,
NaN/Inf/-0 propagation, input mutation, vector-vs-scalar consistency,
out-of-domain parameters, complex symmetries).  Everything numerical now sits
at the input's conditioning floor; the logical findings were real:

| Finding | Was | Fix |
|---|---|---|
| Phase reduction `u - k*pi` rounds `k*pi` (both ports, `elliptic12`/`elliptic3`/`ellipticBDJ`) | `eps*\|u\|` in the reduced phase; near pi/2 at m -> 1 amplified 1e5x into Z and Pi (1e-10) | Cody-Waite split `(u - k*PI_HI) - k*PI_LO` |
| python `ellipj` with `m` outside [0,1] or NaN | returned sn(u \| 0.5) -- the interior placeholder leaked | `check_range` (raise on numpy, NaN mask on device backends) + NaN propagation; also `elliptic12`, `nomeq` |
| MATLAB `elliptic12` with NaN `m` | crashed inside `unique()` ("subscripts must be...") | NaN in, NaN out, excluded from the grouping |
| `elliptic12i(-0.0)` / `(0)` (both) | `eps` from the cot(phi) nudge | exact zero (sign preserved) |
| python `R_J` argument ratios beyond 3e32 | 9e-12 at ratio 1.9e44 | 100 duplications (limit ~4e56) |

Verified clean in the same round: no public function mutates its inputs;
vector and scalar calls agree bit-for-bit over 300 random points including
m in {0, 1, 1e-17, nextafter(1,0)}; matrix shapes are preserved by every
function; empty input returns empty; F(conj u) = conj F(u), F(-u) = -F(u) and
sn(conj u) = conj sn(u) hold exactly.

## Adversarial review rounds 4-5 (theta recurrence, parallel and GPU paths, batch independence, 2026-09-02)

Attack surface: the code paths the unit tests never exercised on a
developer machine (parallel chunking, the GPU kernels under identity
stubs) and the invariant "a value must not depend on what else is in the
batch".

| # | Attack | Result | Fix |
|---|--------|--------|-----|
| 4.1 | theta functions at large argument (`v ~ 1e8`) | `sin(k*v)` with `k*v` as a double product loses `eps*|k v|` (1e-12 at `v = 1e8`) in every series of `theta.m`, `theta_prime.m`, `jacobiThetaEta.m`, `weierstrassZeta.m`, `weierstrassSigma.m`, `theta.py`, `weierstrass.py` | one shared `theta_series.m` (all three MATLAB theta callers) and `theta._trig_start` (Python): `sin/cos((2n+1)v)` by angle-addition recurrence from `sin v, cos v`, so every term carries only the rounding of the reduced argument |
| 4.2 | parallel chunk path: `N` an exact multiple of `chunk_size` | `par_worker` re-entered the parallel dispatcher from inside a worker (stub `get_nworkers` ignored the `parallel` flag) and recursed until SIGILL | recursion guard in `par_worker.m`: the worker forces `elliptic_config('parallel', false)` for the duration of its call (restored in `unwind_protect_cleanup`) |
| 4.3 | chunked vs serial `elliptic3` on 1000 random points | 6 of 1000 differed by one ulp (rel 2.4e-16): the Carlson cores stopped on a whole-vector convergence test, so the number of duplication steps applied to an element depended on its batch mates | per-element convergence in `carlsonRF/RD/RJ.m`: an `active` mask freezes each element at its own converged step |
| 4.4 | same attack on `elliptic12` and `ellipj` | the final AGM row used was `max(n)` over the batch (`a(mn,:)`), so K and the Landen back-substitution scale changed with batch composition | first converged row per element (`a(n+1)`) in both CPU and GPU paths |
| 4.5 | GPU `elliptic12` at `u = 1e6`, `m = 1 - eps/2` vs CPU | Z off by 3.6e-11, E by 1.2e-10: the GPU phase reduction lacked the Cody-Waite tail term (`k * 1.22e-16 = 3.9e-11` at `k = 318310`) | tail term added to the GPU reduction, matching the CPU path |

New test: `testParallel.m` block "chunking is exact" builds temporary
`get_nworkers`/`parcellfun` stubs, evaluates every parallel-capable
function serially and chunked for `N` in `{cs-1, cs, 2cs, 3cs, 3cs+7}`
(`cs = chunk_size`) and requires bit-identical results.  The GPU-stub
probe (identity `gpuArray/gather`) now agrees with the CPU path to 0.0 on
4005 points including the large-`u` near-`m = 1` cases.

## Adversarial review round 6 (cross-port parity sweep at extreme m, 2026-09-02)

Attack: 3000 random points with `m` drawn from `U(0,1)`, `1 - 10^U(-16,-3)`
and `10^U(-16,-3)`, `u` from `U(-3,3)` and `U(-1e5,1e5)`, characteristic
`n` from `U(0,0.999)`; both ports evaluated on the identical doubles and
every disagreement above 1e-13 was adjudicated with mpmath at the exact
inputs.  Eight defects, several older than the April refactor:

| # | Function | Defect | Fix |
|---|----------|--------|-----|
| 6.1 | `elliptic12.m` (CPU + GPU) | `m` in `[eps^2, ~5e-16]`: the AGM converges in one step, no Landen step ran, the scale `e` stayed 0 and `F = E = Inf` | scale `2^(n-2)` in closed form (`n <= 1 -> 1/2`) |
| 6.2 | `theta.m`, `theta_prime.m`, `jacobiThetaEta.m` | nome from `ellipke(1-m)`: `1-m` rounds first, `q` was 30% off at `m ~ 1e-16` and `theta1` off by 1e-5 | `K'(m) = R_F(0, m, 1)` from the exact `m`, as `nomeq` already did |
| 6.3 | `elliptic12.m` (CPU + GPU) | `E/K = 1 - sum 2^(j-1) c_j^2` stopped one AGM term early (`2^(n-2) c_{n-1}^2 ~ 1e-14` near `m -> 1`); `E(-1.65|1-4e-15)` off by 1.4e-13 | the C sum takes `n` terms, the descent still `n-1` |
| 6.4 | all `k*pi` reductions (both ports) | the "Cody-Waite" split used `double(pi)` as the head, so `k*pi` itself rounded by `eps*|u|` (2.3e-10 at `u = 1e6`); Jacobi Zeta at `u = 8e4` was off by 1.5e-11 | `sub_kpi` / `_xputils.sub_kpi`: 25-bit `PI_A`, `PI_B` (products exact for `k < 2^28`) plus `PI_C`; verified to 4e-16 over 20000 random `k < 2^27` |
| 6.5 | `elliptic3.m` | the reflection `Pi(pi-u) = 2 Pi(pi/2) - Pi(u)` used `elliptic3(double(pi/2))`: `cos(double(pi/2)) = 6e-17`, and near `m = 1` the sliver `6e-17/(sqrt(1-m)(1-c))` is 2e-7 (relative 3e-10) | complete integral from the exact Carlson form `R_F(0,1-m,1) + (c/3) R_J(0,1-m,1,1-c)` |
| 6.6 | `carlsonRJ.m`, `carlson.py` | series term `E3 = XYZ + 2 E2 P + 3 P^3`; DLMF 19.36.2 has `4 P^3`.  With the 0.0015 stopping tolerance the O(eps^4) residual was 1e-13 relative (Python masked it by running 100 duplications) | coefficient corrected; checked in exact arithmetic: residual 3e-22 |
| 6.7 | `ellipticBDJ` (both ports) | `n > 1` with the phase beyond the pole: MATLAB returned complex `J` silently, Python 1.147 (principal value 0.859); `n = 1` returned NaN from `0 * inf` in the period term | error (MATLAB / NumPy) or NaN (traced backends) beyond the pole; the complete `J(n|m)` is only added where a period was removed |
| 6.8 | `elliptic3.m` | rejected `c < 0` although the integral is standard there and the Python port accepts it | `c <= 1` accepted; `c < 0` routed to the Carlson branch (the 20-node rule loses digits as `1/(1+|c| sin^2)` narrows) |

After the fixes the same 3000-point sweep agrees with mpmath to `< 2e-15`
relative on F, E, Z, Pi, theta1 in both ports (previously up to 3e-10, 1e-5
for theta1).  The remaining cross-port gap is `sn, cn, dn` at `|u| ~ 1e5`
(1e-11): the `4K` period is not a constant, so `u - 4kK` rounds by
`eps*|u|` in both ports; this is the documented limit of `ellipj`.

| 6.9 | `cel` (both ports) | evaluated through `m = 1 - kc^2`, which loses `kc` entirely below ~1e-8: `cel1(1e-9)` was Inf (MATLAB) / 2e6 (Python) against `ln(4/kc) = 22.1`; MATLAB also rejected `kc > 1` (`m < 0`) and both returned Inf for `p < 0` | Bulirsch's own kc-native algorithm (Numer. Math. 13, 1969) in both ports: any real `kc`, `p < 0` is the Cauchy principal value (`= Re Pi(1-p | m)`, checked against mpmath), bit-identical across ports, 1e-16 from `kc = 1e-300` to 100 |
| 6.10 | `elliptic12i.py` (Jacobi Zeta output) | complete `K`, `E` taken as `F(double(pi/2)|m)`: `cos(double(pi/2)) = 6e-17`, K 5.8e-9 relative short at `m = 1-eps/2`, Z off by 1.8e-11 | exact Carlson complete forms `R_F(0,1-m,1)`, `R_F - (m/3) R_D` |
| 6.11 | `ellipticBDJ` (both ports) | `Delta^2 = 1 - m sin^2 phi` cancels near `phi = pi/2` as `m -> 1` (relative 2.5e-9 at `m = 1-1e-8`), which `R_D` turned into 3e-10 in `D(phi|m)` | `Delta^2 = (1-m) + m cos^2 phi` |
| 6.12 | `jacobiEDJ` (both ports) | took `am(u)` at `|u| ~ 1e3` (rounding `eps*|am|`) and only then reduced, where the map `phi -> D` is steep (`1/sqrt(1-m)`): `D_u(1520|1-1e-8)` off by 3e-10 on top of 6.11 | reduce `u` by `2K` first, amplitude of the reduced argument, add `2k` times the complete integrals; now at the `eps*|u|` floor (`dD_u/du = sn^2 <= 1`) |
| 6.13 | `arclength_ellipse.m` | `if (a < b) ... elseif (a > b)` on arrays uses all-elements semantics: any mixed array fell through to the circle formula `a (theta1 - theta0)` for every element | elementwise masks after broadcasting scalars (the Python port was already elementwise) |
| 6.14 | `elliptic12i` (both ports) | the period term `pi*ceil(phi/pi - 0.5 + eps)` (Python: `+ 1e-14`) was counted from a separately rounded quantity: for `phi` a few ulps (Python: 3e-14) below `pi/2` it added a period the sign term `(-1)^floor(2phi/pi)` had not crossed, and `Re F` came out `3K` instead of `K`.  `asin(sqrt(3))` lands 2 ulps below `pi/2`, so `elliptic123` inherited `K(3) = 3.003` (mpmath 1.001) | period `pi*ceil(k/2)` from the same `k = floor(2phi/pi)` |
| 6.15 | `elliptic123.m` (complete `m > 1`) | evaluated `elliptic12i(asin(sqrt(m)), 1/m)`, i.e. exactly on the branch point of `F(.|1/m)`, where the decomposition is `sqrt(eps)`-conditioned (1e-8 after 6.14) | DLMF 19.7.3 closed forms `K(m) = (K(1/m) - i K(1-1/m))/sqrt(m)` and the matching `E`; `elliptic123(pi/2, m)` routes there too |
| 6.16 | `inversenomeq` (both ports), `nome2m.m` | above `q_max = 0.7789534` the 30-term theta series is not converged: MATLAB returned `m > 1` (1.034 at `q = 0.999`), Python raised; `nome2m` captured its whole input array in the `fzero` objective and errored on any array (bracket also covered only `q < 0.62`) | the true `1 - m` is below `eps/2` there, so both ports return exactly 1 (clamped `<= 1` below); `nome2m` is an alias of `inversenomeq` |

Deliberate limits found in this sweep: complex `sn, cn, dn` near the poles
`u = iK'` and Weierstrass functions near lattice points carry the
conditioning `eps*K'/|u - iK'|` of the rounded half-period (1e-11 at a
distance 0.05); `elliptic12i` exactly at its branch point
`pi/2 + i acosh(1/sqrt(m))` is `sqrt(eps)`-conditioned (1e-8) because the
function has a square-root singularity there; `theta` at `|v| ~ 1e3` and
`jacobiThetaEta` at `|u| ~ 1e4` carry `eps*|v|` from the argument itself.

Also in this round: every MATLAB docstring `Example:` block now runs as a
test (`testDocExamples.m`) and the Python docstrings run under
`pytest --doctest-modules` (one example printed a 0-d array and was fixed).

## Deliberate limits and residual risk

- CUDA/OpenCL hardware was not available during this audit. GPU source paths
  and dispatch tests were reviewed, and CPU/Octave tests passed, but the
  corrected MATLAB/Octave GPU kernels still require hardware execution before
  a release claim can include direct GPU validation.
- MATLAB Parallel Computing Toolbox was not available. Octave reported its
  parallel-package skips; serial behavior and parallel dispatch code were
  reviewed, but a real multi-worker run remains release validation work.
- Weierstrass functions currently support real inputs only. They now reject
  complex input explicitly instead of silently discarding data.
- `elliptic3` deliberately rejects real paths that cross a third-kind pole;
  Cauchy principal-value continuation is not implemented.
- Near-pole and near-lattice conditioning: `F(phi|m)` with `m -> 1` at
  `phi -> pi/2` and the Weierstrass functions within `~1e-9 omega1` of a
  lattice point are evaluated to the input's conditioning floor
  (`2 eps |z| / |z - 2k omega1|`, i.e. ~1e-6 relative at `1e-9 omega1`).
  This is a property of the double input, not of the algorithm; the same
  inputs move the true value by that much.
- `R_J` in the python port uses 60 fixed duplications (JAX-traceable), valid
  for max/min argument ratios up to ~3e32; MATLAB iterates adaptively.
- Jacobi phase reduction is double precision: the residual phase carries an
  absolute uncertainty ~`|u|*eps`, so `ellipj` holds full precision to
  `|u| ~ 1e12`, degrades linearly beyond, and has lost the phase entirely by
  `|u| ~ 1e16`.  Every double-precision implementation (MATLAB's and SciPy's
  included) shares this bound; extended-precision reduction would need K(m)
  to ~32 digits.
- `elliptic12i` follows the A&S 17.4.11 real-decomposition branch: on
  `Re u = pi/2` above the branch point, `Re F = K(m)`, which diverges as
  `m -> 1`.  mpmath/Mathematica may return values on a different sheet
  there; the convention is now documented in both ports.
