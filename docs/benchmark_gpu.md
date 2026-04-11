# Benchmark: Serial CPU / Parallel CPU / GPU

**Platform:** Intel i7-7700K @ 4.20 GHz (8 cores) · 32 GB RAM  
**GPU:** NVIDIA GeForce GTX 1080 Ti (11264 MiB VRAM, compute 6.1)  
**Software:** Octave 6.4.0 · ocl 1.2.4 (OpenCL 3.0) · parallel 4.0.2  
**Driver:** NVIDIA 535.288.01  
**Date:** 2026-04-11  
**Methodology:** 3 repetitions per cell, minimum wall-clock time reported.
GPU time includes host↔device memory transfer. Parallel uses 8 workers
(`parcellfun`). Parallel chunk size = `ceil(N / nWorkers)`.

---

## elliptic12

`[F,E,Z] = elliptic12(u, m)` — AGM / Descending Landen transformation

| N | Serial (s) | Parallel (s) | par/serial | GPU (s) | gpu/serial | GPU Mpts/s |
|---|---|---|---|---|---|---|
| 10 000 | 0.0086 | 0.0382 | 0.2× | 0.0083 | 1.0× | 1.2 |
| 50 000 | 0.0483 | 0.0456 | 1.1× | 0.0139 | 3.5× | 3.6 |
| 100 000 | 0.1004 | 0.0662 | 1.5× | 0.0253 | 4.0× | 4.0 |
| 500 000 | 0.5640 | 0.2643 | 2.1× | 0.1442 | 3.9× | 3.5 |
| 1 000 000 | 1.2963 | 0.5293 | **2.4×** | 0.2798 | **4.6×** | 3.6 |
| 2 000 000 | 2.5610 | 1.0383 | 2.5× | 0.5569 | 4.6× | 3.6 |
| 4 000 000 | 5.1140 | 2.0912 | 2.4× | 1.1256 | 4.5× | 3.6 |

**Hardware utilisation @ N = 1 M:**  
Serial CPU 16% (1 core) · Parallel CPU 90% (8 cores) · GPU compute 38% · GPU mem 596 MiB

---

## ellipj

`[sn,cn,dn,am] = ellipj(u, m)` — AGM / Ascending Landen back-substitution

| N | Serial (s) | Parallel (s) | par/serial | GPU (s) | gpu/serial | GPU Mpts/s |
|---|---|---|---|---|---|---|
| 10 000 | 0.0066 | 0.0364 | 0.2× | 0.0060 | 1.1× | 1.7 |
| 50 000 | 0.0398 | 0.0431 | 0.9× | 0.0128 | 3.1× | 3.9 |
| 100 000 | 0.0812 | 0.0624 | 1.3× | 0.0236 | 3.4× | 4.2 |
| 500 000 | 0.4828 | 0.2348 | 2.1× | 0.1430 | 3.4× | 3.5 |
| 1 000 000 | 1.0886 | 0.4653 | **2.3×** | 0.2852 | **3.8×** | 3.5 |
| 2 000 000 | 2.1787 | 0.9278 | 2.3× | 0.5733 | 3.8× | 3.5 |
| 4 000 000 | 4.4240 | 1.8575 | 2.4× | 1.1441 | 3.9× | 3.5 |

**Hardware utilisation @ N = 1 M:**  
Serial CPU 16% · Parallel CPU 85% (8 cores) · GPU compute 27% · GPU mem 472 MiB

---

## elliptic3

`Pi = elliptic3(u, m, c)` — 10-point Gauss-Legendre quadrature (no AGM)

| N | Serial (s) | Parallel (s) | par/serial | GPU (s) | gpu/serial | GPU Mpts/s |
|---|---|---|---|---|---|---|
| 10 000 | 0.0037 | 0.0405 | 0.1× | 0.0038 | 1.0× | 2.7 |
| 50 000 | 0.0150 | 0.0395 | 0.4× | 0.0055 | 2.7× | 9.1 |
| 100 000 | 0.0319 | 0.0475 | 0.7× | 0.0078 | 4.1× | 12.8 |
| 500 000 | 0.2122 | 0.1393 | 1.5× | 0.0306 | 6.9× | 16.3 |
| 1 000 000 | 0.4451 | 0.2887 | **1.5×** | 0.0564 | **7.9×** | 17.7 |
| 2 000 000 | 1.2056 | 0.5827 | 2.1× | 0.0995 | 12.1× | 20.1 |
| 4 000 000 | 2.6679 | 1.1673 | 2.3× | 0.1988 | **13.4×** | 20.1 |

**Hardware utilisation @ N = 1 M:**  
Serial CPU 16% · Parallel CPU 89% (8 cores) · GPU compute 41% · GPU mem 244 MiB

`elliptic3` is the clear GPU winner: all 10 quadrature iterations are
embarrassingly parallel (no `gather` synchronisation between iterations),
so the GPU stays busy. At N = 4 M it sustains **20 Mpts/s** — 13× faster
than serial.

---

## jacobiThetaEta

`[Th,H] = jacobiThetaEta(u, m)` — AGM + theta product

| N | Serial (s) | Parallel (s) | par/serial | GPU (s) | gpu/serial | GPU Mpts/s |
|---|---|---|---|---|---|---|
| 10 000 | 0.0099 | 0.0378 | 0.3× | 0.0085 | 1.2× | 1.2 |
| 50 000 | 0.0503 | 0.0480 | 1.0× | 0.0193 | 2.6× | 2.6 |
| 100 000 | 0.1116 | 0.0704 | 1.6× | 0.0384 | 2.9× | 2.6 |
| 500 000 | 0.7014 | 0.2944 | 2.4× | 0.2444 | 2.9× | 2.0 |
| 1 000 000 | 1.4598 | 0.5832 | **2.5×** | 0.4935 | **3.0×** | 2.0 |
| 2 000 000 | 2.9953 | 1.2110 | 2.5× | 0.9996 | 3.0× | 2.0 |
| 4 000 000 | 6.1122 | 2.4190 | 2.5× | 1.9750 | 3.1× | 2.0 |

**Hardware utilisation @ N = 1 M:**  
Serial CPU 15% · Parallel CPU 85% (8 cores) · GPU compute 28% · GPU mem 580 MiB

`jacobiThetaEta` carries the most AGM state per element (the `prodth`
product array) and thus uses the most GPU memory.

---

## Summary

| Function | Best parallel speedup | Best GPU speedup | GPU Mpts/s (peak) |
|---|---|---|---|
| `elliptic12` | 2.5× | **4.6×** | 3.6 |
| `ellipj` | 2.4× | **3.9×** | 3.5 |
| `elliptic3` | 2.3× | **13.4×** | **20.1** |
| `jacobiThetaEta` | 2.5× | **3.1×** | 2.0 |

### Why GPU < 8× for AGM functions

AGM-based functions (`elliptic12`, `ellipj`, `jacobiThetaEta`) require
`gather()` inside the while loop to check convergence on the CPU:

```matlab
while any(gather(abs(c(:,ii))) > tol)
```

Each iteration synchronises the GPU with the CPU (~5–7 iterations total).
This synchronisation overhead, combined with the relatively small
arithmetic intensity, limits GPU occupancy to 27–41%.

`elliptic3` has no convergence loop — the 10 Gauss-Legendre iterations
execute entirely on the GPU without synchronisation, reaching 41% and
delivering 13× speedup.

### Parallel crossover point

The parallel path incurs ~35 ms dispatch overhead (parcellfun spawn cost
in Octave).  It only beats serial when the computation exceeds this:

- `elliptic12`, `ellipj`, `jacobiThetaEta`: crossover ≈ N = 50 000–100 000
- `elliptic3`: crossover ≈ N = 300 000 (faster serial path raises bar)

Below the crossover, always use serial (default) or GPU.

---

## Reproducing these results

```octave
pkg load parallel ocl
addpath src
bench_gpu
```

Raw data: [`bench_gpu_results.csv`](../bench_gpu_results.csv)  
Full GPU setup guide: [`GPU.md`](GPU.md)
