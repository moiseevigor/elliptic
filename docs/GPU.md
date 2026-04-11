# GPU Acceleration

This library supports GPU-accelerated evaluation of all four core functions
(`elliptic12`, `ellipj`, `elliptic3`, `jacobiThetaEta`) through an opt-in
configuration flag.  The same source code paths work in both **MATLAB**
(Parallel Computing Toolbox / CUDA) and **Octave** (ocl Forge package /
OpenCL).

---

## Benchmark results

Measured on: Intel i7-7700K (8 cores, 4.2 GHz) · GTX 1080 Ti (11 GB) ·
Octave 6.4.0 · ocl 1.2.4 · CUDA driver 535, OpenCL 3.0.

| Function | N | Serial | Parallel (8 cores) | GPU (OpenCL) | GPU speedup |
|---|---|---|---|---|---|
| `elliptic12` | 1 M | 1.296 s | 0.529 s (2.4×) | 0.280 s | **4.6×** |
| `ellipj` | 1 M | 1.089 s | 0.465 s (2.3×) | 0.285 s | **3.8×** |
| `elliptic3` | 4 M | 2.668 s | 1.167 s (2.3×) | 0.199 s | **13.4×** |
| `jacobiThetaEta` | 1 M | 1.460 s | 0.583 s (2.5×) | 0.494 s | **3.1×** |

`elliptic3` achieves 13× because it is a pure Gauss-Legendre quadrature
(no AGM, no sequential dependencies) — every element is completely
independent and maps trivially to GPU threads.

### Hardware utilisation at N = 1 M

| Mode | CPU load | GPU compute | GPU memory |
|---|---|---|---|
| Serial | 16% (1 core) | — | — |
| Parallel | 85–90% (8 cores) | — | — |
| GPU (OpenCL) | ~15% (host overhead) | 27–41% | 244–596 MiB / 11264 MiB |

GPU compute utilisation is bounded by the AGM convergence loop, which
involves repeated `gather()` calls to check CPU-side convergence — these
synchronisation points limit occupancy.  `elliptic3` reaches 41% because
it has no such synchronisation.

---

## Requirements

### MATLAB

| Requirement | Version |
|---|---|
| MATLAB | R2019b or later |
| Parallel Computing Toolbox | any |
| CUDA-capable GPU | compute capability ≥ 3.5 |
| CUDA Toolkit | compatible with your MATLAB release |

No extra packages needed — MATLAB's `gpuArray` / `gather` are part of
the Parallel Computing Toolbox.

Verify your GPU:
```matlab
gpuDeviceCount          % must be > 0
gpuDevice               % shows device name, compute capability
```

### Octave

| Requirement | Version |
|---|---|
| Octave | 6.0 or later |
| ocl package | 1.2.4 or later |
| OpenCL runtime | any vendor (NVIDIA, AMD, Intel) |
| NVIDIA GPU | any with OpenCL support |

---

## Installation

### Octave — NVIDIA GPU via OpenCL

#### 1. Install the OpenCL ICD loader and NVIDIA ICD

```bash
# Debian / Ubuntu
sudo apt-get install ocl-icd-libopencl1 ocl-icd-opencl-dev clinfo

# Verify OpenCL sees the GPU
clinfo | grep -E "Platform|Device Name|Max compute"
```

If `clinfo` shows no platforms, the NVIDIA OpenCL ICD is missing.  It is
typically installed as part of the NVIDIA driver:

```bash
# Check the ICD file location
ls /etc/OpenCL/vendors/
# Expected: nvidia.icd  (created by the NVIDIA driver installer)
```

If the file is missing, reinstall the driver with:
```bash
sudo apt-get install --reinstall nvidia-opencl-icd-<version>
# or re-run the NVIDIA .run installer and tick the OpenCL option
```

#### 2. Install the Octave `ocl` package

```octave
pkg install -forge ocl
```

If the Forge version does not compile (OpenCL header path issues):
```bash
# Install headers first
sudo apt-get install opencl-headers

# Then retry
octave --eval "pkg install -forge ocl"
```

Manual installation from source (latest):
```bash
git clone https://github.com/alexej-schelle/ocl.git
cd ocl
octave --eval "pkg install ."
```

#### 3. Load and verify

```octave
pkg load ocl
a = gpuArray(rand(1, 1000));
b = sqrt(a);
assert(max(abs(gather(b) - sqrt(gather(a)))) < 1e-12)
disp('OCL working')
```

---

### MATLAB — NVIDIA GPU via CUDA

MATLAB handles everything through the Parallel Computing Toolbox.  No
separate OpenCL setup is needed.

```matlab
% Check that a GPU is visible
gpuDeviceCount              % → 1 (or more)
d = gpuDevice(1);           % select GPU 1
fprintf('Using %s (compute %g)\n', d.Name, d.ComputeCapability)

% Quick smoke test
a = gpuArray(rand(1, 1000));
gather(sqrt(a));
disp('CUDA gpuArray working')
```

---

## Enabling GPU acceleration

GPU mode is **disabled by default**.  Enable it once per session before
calling any elliptic function:

```matlab
% Octave: load the package first
pkg load ocl          % Octave only

% Enable GPU
addpath src           % or wherever the library lives
elliptic_config('gpu', true)

% Now all function calls automatically use the GPU
[F, E, Z] = elliptic12(u, m);
Pi        = elliptic3(u, m, c);
[sn, cn]  = ellipj(u, m);
[Th, H]   = jacobiThetaEta(u, m);

% Disable again
elliptic_config('gpu', false)
```

### Check GPU status

```matlab
has_gpu()             % returns true when config=true AND hardware present
elliptic_config('gpu')  % returns current setting (true/false)
```

### Combining GPU and parallel

GPU and parallel modes are mutually exclusive — GPU takes precedence:

```matlab
elliptic_config('gpu',      true)
elliptic_config('parallel', true)   % ignored when gpu=true
```

---

## Performance tips

### When to use GPU

| Scenario | Recommendation |
|---|---|
| N < 50 000 | Serial CPU — GPU transfer overhead dominates |
| 50 000 ≤ N < 200 000 | GPU starts to break even (2–4×) |
| N ≥ 500 000 | GPU recommended (3–13×) |
| `elliptic3` any large N | GPU strongly recommended (10–13×) |

### Memory budget

The GPU path allocates `3 × N × 12 × 8` bytes (three `N × 12` double
arrays for the AGM).  For a GTX 1080 Ti with 11 GB VRAM:

| N | GPU memory |
|---|---|
| 1 M | ~300 MiB |
| 4 M | ~1.2 GiB |
| 10 M | ~3 GiB |
| 30 M | ~8.6 GiB (near limit) |

`elliptic3` uses much less (no AGM arrays, just three N-vectors).

### Warm-up calls

The first GPU call compiles OpenCL kernels and incurs a one-time overhead
of 0.5–2 seconds.  For benchmarking, always run one warm-up call before
measuring:

```matlab
elliptic_config('gpu', true)
elliptic12(u, m);          % warm-up (kernel compilation)
tic; elliptic12(u, m); toc % actual measurement
```

---

## Troubleshooting

### `has_gpu()` returns false in Octave

1. Check the `gpu` config is set: `elliptic_config('gpu', true)`
2. Check the ocl package is loaded: `pkg load ocl`
3. Test OCL directly:
   ```octave
   pkg load ocl; gpuArray(1.0)
   ```
4. Check OpenCL sees the device: run `clinfo` in a terminal

### `OclArray: dimensions of both arrays must match exactly`

This error should not occur in the current codebase.  If you see it, make
sure you are running the latest version of the library (`git pull`).  It
was a known issue fixed by switching to column-vector layout throughout
the GPU subfunctions.

### NVIDIA driver shows OpenCL version < 3.0

The ocl package requires at least OpenCL 1.2.  Drivers ≥ 384 support
OpenCL 3.0.  Update via:
```bash
sudo apt-get install nvidia-driver-<latest>
```

### Octave `pkg install -forge ocl` fails to compile

```bash
# Check that opencl headers are installed
dpkg -l | grep opencl-headers
# Install if missing
sudo apt-get install opencl-headers ocl-icd-opencl-dev
# Retry
octave --eval "pkg install -forge ocl"
```

### MATLAB: `gpuArray` not defined

The Parallel Computing Toolbox is not installed or not licensed.  Run
`ver` in MATLAB to confirm.  Without the toolbox, the library falls back
to CPU automatically (`has_gpu()` returns false).

---

## Running the benchmark

```matlab
% Octave
pkg load parallel ocl
addpath src
bench_gpu

% MATLAB (with parpool already started)
addpath src
bench_gpu
```

Results are written to `bench_gpu_results.csv` with columns:
`function, N, t_serial_s, t_parallel_s, t_gpu_s, speedup_par, speedup_gpu, mpts_per_s_gpu`
