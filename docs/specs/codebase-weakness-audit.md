# Codebase Weakness Audit — Python Elliptic Library

**Date:** 2026-04-21
**Scope:** `/home/igor/Work/elliptic/python/elliptic/`
**Purpose:** Identify weak or breaking points (special values, complex numbers, inf/NaN, memory, exception handling, backend tracing) to guide future implementation work.

File:line refs are approximate; verify before acting.

---

## HIGH severity — tracing failures or silent data loss

- **`complex_elliptic.py:40-41, 124-125`** — `if np.any(m < 0) or np.any(m > 1)` on possibly-traced JAX arrays. Breaks JAX tracing outright.
- **`weierstrass.py:34-37`** — `xp.asarray(..., dtype=xp.float64)` silently drops imaginary parts; no complex-input validation. Also `weierstrassZeta` / `weierstrassSigma` force `.ravel()` via numpy (lines 90-92, 146-148), breaking JAX/torch dispatch entirely.
- **`inverse.py:62`** — `if np.max(np.abs(res)) < tol: break` inside Newton loop. Python-level control flow on a possibly-traced array → JAX trace failure.
- **`theta.py:195-197, 243-246`** — `if np.any(maskN)` gates under tracing.
- **`applications.py:46-47, 56-61`** — `float(a)`, `float(E1)`, `float(E0)` casts reject traced arrays with an opaque error.

## MEDIUM — numerical / correctness

- **`carlson.py:78-89, 113-133, 169-193`** — fixed 20/30 iterations. At m very close to 1 (e.g. `1 − 1e-15`) convergence is unverified; silent accuracy loss. No runtime tolerance check.
- **`carlson.py:34-51` (`_rc_xp`)** — all three branches (`rc_gt`, `rc_lt`, `rc_eq`) are computed then `where`-selected; `x_safe = 1.0` fallback when `x < 1e-300` can still feed non-selected branches. Masked NaNs can propagate in reverse-mode autodiff.
- **`complex_elliptic.py:56-65`** — discriminant `sqrt(max(b²/4 − c, 0))` followed by `1/sqrt(...)`: catastrophic cancellation at extreme m; result can go subnormal and divide blows up.
- **`inverse.py:64-65`** — Newton divides by `max(denom, 1e-15)`. At true pole (m=1, φ=π/2) the guard masks non-convergence; Newton stalls and returns wrong φ without signalling.
- **`nome.py:56-65`** — `q_max` computed from `nextafter(1, 0)` is slightly below the true bound; boundary rejections are brittle. `nome.py:70-73` uses a Python `for` + `brentq` over flattened array — destroys batched JAX/torch perf (not a crash, but severe regression).
- **`ellipj.py:36-37`** — `m_safe = clip(m, 1e-15, 1-1e-15)` is a design hazard: correctness relies on every later branch re-patching with `where`; any new codepath using `m_safe` directly is silently wrong.
- **`elliptic12.py:82`** — `log(tan(π/4 + u/2))` for m=1 produces NaN near u=π/2 (tan overflow before log). `near_pole_m1` guard exists but doesn't cover all subnormal edges.
- **`weierstrass.py:134, 178`** — `result = inf` when `|z_reduced| < 1e-5`: heuristic, not a true pole; returns `inf` for large-but-finite values.

## LOW

- **`theta.py:137-139, 207-208`** — returns NaN at m ≈ 1 silently; inconsistent with `elliptic3.py` which raises. Consider raising for uniformity.
- **`_agm.py:46-49`** — bare `except Exception` around JAX stacking swallows real errors.
- **`ellipticBD.py:41-47`** — redundant `where(m==0, 1, m)` inside already-masked branch; cosmetic only.

---

## Cross-cutting recommendations

1. **JAX safety pass** — ban `np.any(...)` in control flow, ban `float()` / `.item()` on inputs, ban in-place assignment, ban Python `for` over flattened arrays. Validation should live behind a `not _is_traced(...)` guard or use pure `where`.
2. **Complex dispatch** — pick a canonical complex promotion rule. Either accept complex and promote to `complex128`, or reject early with `ValueError` — never silently drop imaginary parts via `float64` cast.
3. **Carlson convergence check** — after the fixed iteration, assert `max(|X|, |Y|, |Z|) < threshold` (or add 5-10 more iters) so users near m=1 get accuracy, not silent ~1e-8 error.
4. **Fused `where` hygiene** — for `where(cond, safe, unsafe)` in autodiff, guard the unsafe branch so it doesn't compute NaN gradients. Use the double-where trick: `safe_x = where(cond, x, 1); y = f(safe_x); out = where(cond, y, fallback)`.
5. **Pole handling policy** — document one uniform rule: raise `ValueError` vs return `inf` / `nan`. Today it's mixed: `elliptic3` raises, `theta` returns NaN, `weierstrass` returns `inf`.
6. **Replace `m_safe` pattern** — prefer computing the unsafe expression inside `where(cond, ..., fallback)` without mutating `m` globally, to remove the refactoring hazard.

---

## Suggested priority order for future PRs

1. JAX tracing fixes (HIGH cluster): `complex_elliptic`, `inverse`, `theta`, `applications` validation paths.
2. Complex-argument policy + `weierstrass.py` backend restore.
3. Carlson convergence guards + where-branch NaN hygiene.
4. Unify pole handling policy across `elliptic3` / `theta` / `weierstrass`.
5. Nome batched path (replace scalar `brentq` loop).
6. Retire `m_safe` clip pattern in `ellipj.py`.
