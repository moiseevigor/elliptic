# Benchmark Comparison: Before vs After Vectorization

**Platform:** Octave 6.4.0 on x86_64-pc-linux-gnu
**Date:** 2026-04-04
**Methodology:** 5 iterations per size, 1 warmup run excluded. Mean and std reported.

## Changes tested

- **elliptic12.m**: `find()` + `uint32` → logical indexing, precomputed `e_vals` in Landen loop
- **elliptic3.m**: 10-iteration Gauss-Legendre while loop → single vectorized matrix operation
- **ellipj.m**: `find()` + `uint32` → logical indexing in AGM descent loop
- **jacobiThetaEta.m**: `find()` + `uint32` → logical indexing in AGM descent loop

## Results

### elliptic12

| Size | Points | Before (s) | After (s) | Speedup |
|------|--------|-----------|----------|---------|
| 100x100 | 10,000 | 0.0971 | 0.0989 | 0.98x |
| 500x500 | 250,000 | 2.5539 | 2.6064 | 0.98x |
| 1000x1000 | 1,000,000 | 10.2694 | 10.1719 | 1.01x |
| 2000x2000 | 4,000,000 | 41.4147 | 40.8450 | 1.01x |

**Verdict: No significant change (~1%).** The logical indexing optimization in the Landen loop and AGM convergence loop provides negligible benefit. The bottleneck is elsewhere (likely the AGM iteration itself and `uniquetol_compat` overhead).

### elliptic3

| Size | Points | Before (s) | After (s) | Speedup |
|------|--------|-----------|----------|---------|
| 100 | 100 | 0.000456 | 0.000233 | **1.96x** |
| 1,000 | 1,000 | 0.000758 | 0.000509 | **1.49x** |
| 10,000 | 10,000 | 0.003333 | 0.003312 | 1.01x |
| 100,000 | 100,000 | 0.032550 | 0.046247 | **0.70x** |

**Verdict: Mixed.** The vectorized Gauss-Legendre shows ~2x speedup at small sizes (reduced loop overhead) but is **slower at large sizes** (100k pts). The matrix broadcast `t(:) * u` creates a 10xN intermediate matrix that causes memory pressure at large N, defeating the purpose.

### ellipj

| Size | Points | Before (s) | After (s) | Speedup |
|------|--------|-----------|----------|---------|
| 100x100 | 10,000 | 0.0955 | 0.0976 | 0.98x |
| 500x500 | 250,000 | 2.4927 | 2.5236 | 0.99x |
| 1000x1000 | 1,000,000 | 10.2564 | 10.1287 | 1.01x |
| 2000x2000 | 4,000,000 | 40.4401 | 40.4683 | 1.00x |

**Verdict: No significant change.** Same as elliptic12 — the logical indexing swap is neutral.

### jacobiThetaEta

| Size | Points | Before (s) | After (s) | Speedup |
|------|--------|-----------|----------|---------|
| 100x100 | 10,000 | 0.1018 | 0.1003 | 1.01x |
| 500x500 | 250,000 | 2.5334 | 2.5581 | 0.99x |
| 1000x1000 | 1,000,000 | 10.4582 | 10.3573 | 1.01x |
| 2000x2000 | 4,000,000 | 42.3392 | 41.7325 | 1.01x |

**Verdict: No significant change.**

## Conclusion

The `find()` → logical indexing changes in **elliptic12, ellipj, and jacobiThetaEta** provide **no measurable benefit** on Octave. The inner AGM/Landen loops only iterate ~7 times regardless of input size, so the `find()` overhead is trivial compared to the vectorized arithmetic on large arrays.

The **elliptic3 Gauss-Legendre vectorization** helps at small sizes (~2x) but hurts at large sizes due to 10xN memory allocation. A chunked approach or keeping the original loop may be preferable.

### Bottleneck analysis

All four functions show near-identical ~10s/1M points scaling, suggesting the bottleneck is:
1. **The `uniquetol_compat` fallback** (O(n log n) sort-based, runs on every call)
2. **The vectorized arithmetic itself** (element-wise ops on million-element arrays)
3. **Single-threaded execution** — Octave does not parallelize element-wise operations across cores

### Recommendations

1. **Revert the elliptic3 vectorization** for large inputs, or add a size threshold
2. **Profile `uniquetol_compat`** — it may dominate runtime for large inputs with many unique m values
3. **Keep the logical indexing changes** — they're cleaner code even if not faster
4. **For real speedup**: consider MEX/C implementation of hot loops, or parallel Octave packages
