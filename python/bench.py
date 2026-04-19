"""Comprehensive benchmark for the elliptic Python package.

Measures wall-clock time, throughput (Melem/s), and memory for every
major function across NumPy, PyTorch CPU/CUDA, and JAX backends.

Usage
-----
    python bench.py                   # NumPy only, default sizes
    python bench.py --torch           # + PyTorch CPU (+ CUDA if available)
    python bench.py --jax             # + JAX
    python bench.py --all             # all backends
    python bench.py --sizes 1e4 1e5   # custom input sizes
    python bench.py --reps 5          # repetitions per timing
    python bench.py --csv results.csv # save results
    python bench.py --all --csv out.csv --reps 5

Output columns
--------------
    Function   Backend   N       ms/call   Melem/s   peak_MB
"""
from __future__ import annotations

import argparse
import csv
import gc
import io
import os
import platform
import sys
import time
import tracemalloc
from dataclasses import dataclass, field
from typing import Any, Callable, List, Optional, Tuple

import numpy as np

# ── optional imports ────────────────────────────────────────────────────────
try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False

try:
    import torch
    _HAS_TORCH = True
    _CUDA = torch.cuda.is_available()
except ImportError:
    _HAS_TORCH = False
    _CUDA = False

try:
    import jax
    import jax.numpy as jnp
    _HAS_JAX = True
except ImportError:
    _HAS_JAX = False

try:
    import pynvml
    pynvml.nvmlInit()
    _HAS_NVML = True
except Exception:
    _HAS_NVML = False

import elliptic

# ── data types ───────────────────────────────────────────────────────────────

@dataclass
class Row:
    function:  str
    backend:   str
    N:         int
    ms:        float          # wall time per call (ms)
    melem_s:   float          # million elements / second
    peak_mb:   float          # peak RAM delta (MB); GPU-MB when applicable
    gpu_mb:    float = 0.0    # GPU VRAM delta (MB) — 0 when not applicable


# ── timing helpers ───────────────────────────────────────────────────────────

def _rss_mb() -> float:
    if _HAS_PSUTIL:
        return psutil.Process(os.getpid()).memory_info().rss / 1024**2
    return 0.0


def _gpu_mb() -> float:
    if _CUDA:
        return torch.cuda.memory_allocated() / 1024**2
    return 0.0


def _nvml_gpu_util() -> Optional[int]:
    """GPU utilisation % via pynvml (first device), or None."""
    if not _HAS_NVML:
        return None
    try:
        h = pynvml.nvmlDeviceGetHandleByIndex(0)
        return pynvml.nvmlDeviceGetUtilizationRates(h).gpu
    except Exception:
        return None


def _time_fn(fn: Callable, args: tuple, reps: int) -> Tuple[float, float]:
    """Return (mean_ms, peak_ram_mb) for *fn(*args)* over *reps* calls."""
    # warm-up (also JIT compile for JAX/Torch)
    fn(*args)
    if _HAS_JAX:
        import jax
        jax.block_until_ready(fn(*args))

    gc.collect()
    rss0 = _rss_mb()

    tracemalloc.start()
    t0 = time.perf_counter()
    for _ in range(reps):
        result = fn(*args)
        # block JAX async dispatch
        if _HAS_JAX:
            import jax
            jax.block_until_ready(result)
        # block CUDA async kernels
        if _CUDA and isinstance(result, tuple):
            for r in result:
                if hasattr(r, 'device') and str(getattr(r, 'device', '')) != 'cpu':
                    torch.cuda.synchronize()
                    break
        elif _CUDA and hasattr(result, 'device'):
            torch.cuda.synchronize()
    elapsed = (time.perf_counter() - t0) / reps * 1000  # ms

    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_mb = max(peak / 1024**2, _rss_mb() - rss0)

    return elapsed, peak_mb


def _time_gpu(fn: Callable, args: tuple, reps: int) -> Tuple[float, float, float]:
    """Return (mean_ms, cpu_mb, gpu_mb) using CUDA events for accurate GPU timing.

    *fn* is a zero-argument callable (args are captured in a closure).
    """
    if not _CUDA:
        raise RuntimeError("CUDA not available")

    torch.cuda.reset_peak_memory_stats()
    gc.collect()
    cpu0 = _rss_mb()

    start = torch.cuda.Event(enable_timing=True)
    end   = torch.cuda.Event(enable_timing=True)
    start.record()
    for _ in range(reps):
        fn()
    end.record()
    torch.cuda.synchronize()
    elapsed = start.elapsed_time(end) / reps   # ms

    gpu_peak = torch.cuda.max_memory_allocated() / 1024**2
    cpu_peak = max(0.0, _rss_mb() - cpu0)
    return elapsed, cpu_peak, gpu_peak


# ── input factories ──────────────────────────────────────────────────────────

def _make_numpy(N: int) -> dict:
    rng = np.random.default_rng(42)
    return {
        "phi": rng.uniform(0.01, np.pi / 2, N),
        "m":   rng.uniform(0.01, 0.99,      N),
        "n":   rng.uniform(0.01, 0.99,      N),
        "u":   rng.uniform(0.01, 1.5,       N),
        "e1":  np.full(N, 2.0),
        "e2":  np.zeros(N),
        "e3":  np.full(N, -2.0),
        "kc":  np.sqrt(rng.uniform(0.01, 0.99, N)),
        "a0":  rng.uniform(0.5, 2.0, N),
        "b0":  rng.uniform(0.1, 0.9, N),
    }


def _make_torch(d: dict, device: str = "cpu") -> dict:
    return {k: torch.tensor(v, dtype=torch.float64, device=device)
            for k, v in d.items()}


def _make_jax(d: dict) -> dict:
    return {k: jnp.array(v) for k, v in d.items()}


# ── function registry ────────────────────────────────────────────────────────
# Each entry: (label, fn, lambda d: args_tuple, n_outputs_for_throughput)

BENCHMARKS = [
    ("elliptic12",   elliptic.elliptic12,   lambda d: (d["phi"], d["m"])),
    ("elliptic3",    elliptic.elliptic3,    lambda d: (d["phi"], d["m"], d["n"])),
    ("ellipj",       elliptic.ellipj,       lambda d: (d["u"],   d["m"])),
    ("ellipticBDJ",  elliptic.ellipticBDJ,  lambda d: (d["phi"], d["m"], d["n"])),
    ("ellipticBD",   elliptic.ellipticBD,   lambda d: (d["m"],)),
    ("jacobiEDJ",    elliptic.jacobiEDJ,    lambda d: (d["u"],   d["m"])),
    ("carlsonRF",    elliptic.carlsonRF,    lambda d: (d["m"],   1.0 - d["m"], d["m"] * 0 + 1.0)),
    ("carlsonRD",    elliptic.carlsonRD,    lambda d: (d["m"],   1.0 - d["m"], d["m"] * 0 + 1.0)),
    ("cel1",         elliptic.cel1,         lambda d: (d["kc"],)),
    ("jacobiThetaEta", elliptic.jacobiThetaEta, lambda d: (d["u"], d["m"])),
    ("theta3",       lambda *a: elliptic.theta(3, *a), lambda d: (d["phi"], d["m"])),
    ("weierstrassP", elliptic.weierstrassP, lambda d: (d["u"], d["e1"], d["e2"], d["e3"])),
    ("agm",          elliptic.agm,          lambda d: (d["a0"], d["b0"])),
]

# Subset of functions that work cleanly with Torch/JAX tensors
# (theta and jacobiThetaEta use scipy K internally — numpy path only)
_TORCH_JAX_SKIP = {"jacobiThetaEta", "theta3", "inversenomeq"}


# ── pretty printing ──────────────────────────────────────────────────────────

_COL = {
    "function":  30,
    "backend":   12,
    "N":          8,
    "ms":        10,
    "Melem/s":   10,
    "peak_MB":    9,
    "gpu_MB":     9,
}

def _hdr():
    h = (f"{'Function':<{_COL['function']}}"
         f"{'Backend':<{_COL['backend']}}"
         f"{'N':>{_COL['N']}}"
         f"{'ms/call':>{_COL['ms']}}"
         f"{'Melem/s':>{_COL['Melem/s']}}"
         f"{'peak_MB':>{_COL['peak_MB']}}"
         f"{'gpu_MB':>{_COL['gpu_MB']}}")
    return h + "\n" + "─" * len(h)


def _row_str(r: Row) -> str:
    return (f"{r.function:<{_COL['function']}}"
            f"{r.backend:<{_COL['backend']}}"
            f"{r.N:>{_COL['N']},}"
            f"{r.ms:>{_COL['ms']}.2f}"
            f"{r.melem_s:>{_COL['Melem/s']}.1f}"
            f"{r.peak_mb:>{_COL['peak_MB']}.1f}"
            f"{r.gpu_mb:>{_COL['gpu_MB']}.1f}")


def _size_tag(N: int) -> str:
    if N >= 1_000_000: return f"{N//1_000_000}M"
    if N >= 1_000:     return f"{N//1_000}k"
    return str(N)


# ── system info ──────────────────────────────────────────────────────────────

def _sysinfo() -> str:
    lines = [
        "=" * 70,
        "  elliptic Python benchmark",
        "=" * 70,
        f"  Python      : {sys.version.split()[0]}",
        f"  Platform    : {platform.platform()}",
        f"  NumPy       : {np.__version__}",
    ]
    if _HAS_TORCH:
        lines.append(f"  PyTorch     : {torch.__version__}")
        if _CUDA:
            lines.append(f"  CUDA device : {torch.cuda.get_device_name(0)}")
            lines.append(f"  CUDA mem    : {torch.cuda.get_device_properties(0).total_memory / 1024**3:.1f} GB")
        else:
            lines.append("  CUDA        : not available")
    if _HAS_JAX:
        lines.append(f"  JAX         : {jax.__version__}")
        lines.append(f"  JAX devices : {jax.devices()}")
    if _HAS_PSUTIL:
        cpu = platform.processor() or "unknown"
        lines.append(f"  CPU         : {cpu}")
        lines.append(f"  RAM         : {psutil.virtual_memory().total / 1024**3:.1f} GB")
    lines.append("=" * 70)
    return "\n".join(lines)


# ── main benchmark loops ─────────────────────────────────────────────────────

def run_numpy(sizes: List[int], reps: int) -> List[Row]:
    rows = []
    for N in sizes:
        d = _make_numpy(N)
        for label, fn, make_args in BENCHMARKS:
            try:
                args = make_args(d)
                ms, peak = _time_fn(fn, args, reps)
                rows.append(Row(label, "numpy", N, ms, N / ms / 1000, peak))
            except Exception as e:
                rows.append(Row(label, "numpy", N, float("nan"), 0.0, 0.0))
    return rows


def run_torch_cpu(sizes: List[int], reps: int) -> List[Row]:
    rows = []
    for N in sizes:
        d_np  = _make_numpy(N)
        d_tor = _make_torch(d_np, "cpu")
        for label, fn, make_args in BENCHMARKS:
            if label in _TORCH_JAX_SKIP:
                continue
            try:
                args = make_args(d_tor)
                ms, peak = _time_fn(fn, args, reps)
                rows.append(Row(label, "torch-cpu", N, ms, N / ms / 1000, peak))
            except Exception:
                pass
    return rows


def run_torch_cuda(sizes: List[int], reps: int) -> List[Row]:
    """Run functions natively on CUDA tensors; optionally with torch.compile."""
    _has_compile = hasattr(torch, "compile")

    def _compiled(fn):
        if _has_compile:
            try:
                return torch.compile(fn, fullgraph=False, dynamic=True)
            except Exception:
                return fn
        return fn

    rows = []
    for N in sizes:
        d_np   = _make_numpy(N)
        d_cuda = _make_torch(d_np, "cuda")
        for label, fn, make_args in BENCHMARKS:
            if label in _TORCH_JAX_SKIP:
                continue
            try:
                args = make_args(d_cuda)
                cfn  = _compiled(fn)
                # warm-up: let torch.compile trace the graph
                cfn(*args); torch.cuda.synchronize()
                cfn(*args); torch.cuda.synchronize()
                ms, cpu_mb, gpu_mb = _time_gpu(lambda: cfn(*args), (), reps)
                rows.append(Row(label, "torch-cuda", N, ms, N / ms / 1000, cpu_mb, gpu_mb))
            except Exception as e:
                # fall back: direct call without compile
                try:
                    cfn2 = fn
                    cfn2(*args); torch.cuda.synchronize()
                    ms, cpu_mb, gpu_mb = _time_gpu(lambda: cfn2(*args), (), reps)
                    rows.append(Row(label, "torch-cuda", N, ms, N / ms / 1000, cpu_mb, gpu_mb))
                except Exception:
                    pass
    return rows


def run_jax(sizes: List[int], reps: int) -> List[Row]:
    # Functions use np.asarray() internally, so JAX JIT tracing fails.
    # We pass JAX arrays which convert to numpy — measures dispatch overhead.
    # JAX GPU (XLA) init is attempted; falls back to CPU on CuDNN mismatch.
    rows = []
    for N in sizes:
        d_np = _make_numpy(N)
        try:
            d_jax = _make_jax(d_np)   # triggers XLA GPU init on first call
        except Exception as e:
            print(f"  JAX array creation failed (N={N}): {e}", file=sys.stderr)
            continue
        for label, fn, make_args in BENCHMARKS:
            if label in _TORCH_JAX_SKIP:
                continue
            try:
                args = make_args(d_jax)
                ms, peak = _time_fn(fn, args, reps)
                rows.append(Row(label, "jax", N, ms, N / ms / 1000, peak))
            except Exception:
                pass
    return rows


# ── summary table ────────────────────────────────────────────────────────────

def _speedup_summary(rows: List[Row], sizes: List[int]) -> str:
    """Print per-function speedup of each backend vs numpy at each size."""
    from collections import defaultdict
    # index: (function, N, backend) -> ms
    idx = {(r.function, r.N, r.backend): r.ms for r in rows}

    backends = []
    for r in rows:
        if r.backend != "numpy" and r.backend not in backends:
            backends.append(r.backend)
    if not backends:
        return ""

    lines = ["\n── Speedup vs NumPy ──────────────────────────────────────────────────"]
    hdr = f"{'Function':<22}" + f"{'N':>8}" + "".join(f"{b:>14}" for b in backends)
    lines.append(hdr)
    lines.append("─" * len(hdr))

    fns = list(dict.fromkeys(r.function for r in rows))
    for fn in fns:
        for N in sizes:
            np_ms = idx.get((fn, N, "numpy"))
            if np_ms is None or np_ms != np_ms:  # nan
                continue
            parts = [f"{fn:<22}{N:>8,}"]
            any_other = False
            for b in backends:
                ms = idx.get((fn, N, b))
                if ms and ms == ms and np_ms > 0:
                    parts.append(f"{np_ms/ms:>13.1f}x")
                    any_other = True
                else:
                    parts.append(f"{'—':>14}")
            if any_other:
                lines.append("".join(parts))
    return "\n".join(lines)


# ── CSV output ───────────────────────────────────────────────────────────────

def _write_csv(rows: List[Row], path: str):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["function", "backend", "N", "ms_per_call", "melem_s", "peak_mb", "gpu_mb"])
        for r in rows:
            w.writerow([r.function, r.backend, r.N,
                        f"{r.ms:.4f}", f"{r.melem_s:.2f}",
                        f"{r.peak_mb:.2f}", f"{r.gpu_mb:.2f}"])
    print(f"\n  Results saved → {path}")


# ── entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Benchmark elliptic functions across backends",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--torch",  action="store_true", help="enable PyTorch CPU (+ CUDA if available)")
    parser.add_argument("--jax",    action="store_true", help="enable JAX")
    parser.add_argument("--all",    action="store_true", help="enable all backends")
    parser.add_argument("--sizes",  type=float, nargs="+", default=[1e4, 1e5, 1e6],
                        help="input sizes (e.g. 1e4 1e5 1e6)")
    parser.add_argument("--reps",   type=int, default=3,
                        help="timing repetitions per call")
    parser.add_argument("--csv",    type=str, default="",
                        help="save results to CSV file")
    parser.add_argument("--no-gpu-util", action="store_true",
                        help="skip GPU utilisation polling")
    args = parser.parse_args()

    sizes = [int(s) for s in args.sizes]
    do_torch = args.torch or args.all
    do_jax   = args.jax   or args.all
    do_cuda  = do_torch and _CUDA

    print(_sysinfo())
    print()

    all_rows: List[Row] = []

    # ── NumPy ────────────────────────────────────────────────────────────────
    print("Running NumPy …")
    np_rows = run_numpy(sizes, args.reps)
    all_rows.extend(np_rows)

    # ── PyTorch CPU ──────────────────────────────────────────────────────────
    if do_torch and _HAS_TORCH:
        print("Running PyTorch CPU …")
        all_rows.extend(run_torch_cpu(sizes, args.reps))

    # ── PyTorch CUDA ─────────────────────────────────────────────────────────
    if do_cuda:
        print("Running PyTorch CUDA …")
        all_rows.extend(run_torch_cuda(sizes, args.reps))

    # ── JAX ──────────────────────────────────────────────────────────────────
    if do_jax and _HAS_JAX:
        print("Running JAX …")
        all_rows.extend(run_jax(sizes, args.reps))

    # ── Print results ─────────────────────────────────────────────────────────
    print()
    print(_hdr())

    last_fn = None
    for r in all_rows:
        if r.function != last_fn and last_fn is not None:
            print()  # blank line between function groups
        last_fn = r.function
        print(_row_str(r))

    # ── Speedup table ─────────────────────────────────────────────────────────
    summary = _speedup_summary(all_rows, sizes)
    if summary:
        print(summary)

    # ── CSV ───────────────────────────────────────────────────────────────────
    if args.csv:
        _write_csv(all_rows, args.csv)

    # ── GPU util note ─────────────────────────────────────────────────────────
    if do_cuda and not args.no_gpu_util:
        util = _nvml_gpu_util()
        if util is not None:
            print(f"\n  GPU utilisation after run: {util}%")


if __name__ == "__main__":
    main()
