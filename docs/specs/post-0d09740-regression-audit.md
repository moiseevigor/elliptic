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
