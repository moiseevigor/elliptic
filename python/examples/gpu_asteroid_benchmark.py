#!/usr/bin/env python3
"""GPU benchmark: massive Keplerian propagation + elliptic arc lengths.

Demonstrates throughput of the elliptic package on realistic solar-system
workloads (hundreds of thousands to millions of asteroids × many epochs).

Pipeline per (body, epoch):
    1. mean anomaly  M(t) = M0 + n (t - t0)
    2. Kepler solve  M = E - e sin E
    3. heliocentric state vector (r, v) in inertial frame
    4. elliptic arc length from perihelion:  L(E) = a E(E | e²)
       (incomplete elliptic integral of the second kind, from this package)

Backends benchmarked:
    numpy  — CPU baseline
    torch  — CPU tensor
    torch_cuda — GPU via CUDA
    jax    — JAX if installed (CPU or GPU depending on device)

Honesty caveat:
    This is a Keplerian / osculating-ellipse throughput benchmark.
    It is NOT a full n-body integrator and it is NOT a replacement for
    JPL Horizons for precision ephemerides.

Usage:
    python3 python/examples/gpu_asteroid_benchmark.py \\
        [--bodies 10000,100000,1000000] [--epochs 128,512,2048] \\
        [--dtype float64] [--out FILE]
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import time
from pathlib import Path

import sys

import numpy as np

# Allow running from a checkout without installing the package.
_PKG_PARENT = Path(__file__).resolve().parents[1]
if str(_PKG_PARENT) not in sys.path:
    sys.path.insert(0, str(_PKG_PARENT))

from elliptic import elliptic12, ellipticBD  # noqa: E402

GM_SUN = 0.0002959122082855911  # AU^3 / day^2


# ---------------------------------------------------------------------------
# Backend selection
# ---------------------------------------------------------------------------

def available_backends():
    backends = [("numpy", "cpu")]
    try:
        import torch
        backends.append(("torch", "cpu"))
        if torch.cuda.is_available():
            backends.append(("torch", "cuda"))
    except ImportError:
        pass
    try:
        import jax
        backends.append(("jax", jax.default_backend()))
    except ImportError:
        pass
    return backends


def as_tensor(xp_name, device, arr, dtype):
    if xp_name == "numpy":
        return arr.astype(dtype)
    if xp_name == "torch":
        import torch
        t = torch.as_tensor(arr)
        t = t.to(dtype=getattr(torch, dtype), device=device)
        return t
    if xp_name == "jax":
        import jax.numpy as jnp
        return jnp.asarray(arr, dtype=dtype)
    raise ValueError(xp_name)


def sync(xp_name, device):
    if xp_name == "torch" and device == "cuda":
        import torch
        torch.cuda.synchronize()
    elif xp_name == "jax":
        import jax
        jax.block_until_ready(jax.numpy.zeros(1))


def get_xp(xp_name):
    if xp_name == "numpy":
        return np
    if xp_name == "torch":
        import torch
        return torch
    if xp_name == "jax":
        import jax.numpy as jnp
        return jnp


# ---------------------------------------------------------------------------
# Kepler solve (vectorized Newton–Raphson, fixed iteration count for tracing)
# ---------------------------------------------------------------------------

def kepler_solve(xp, M, e, iters: int = 12):
    E = M + e * xp.sin(M)
    for _ in range(iters):
        E = E - (E - e * xp.sin(E) - M) / (1.0 - e * xp.cos(E))
    return E


# ---------------------------------------------------------------------------
# Full propagation kernel
# ---------------------------------------------------------------------------

def propagate(xp_name, device, a, e, inc, Om, w, M0, epoch_jd, t_grid):
    """Return position tensor of shape (N_bodies, N_epochs, 3), plus E."""
    xp = get_xp(xp_name)
    N, T = a.shape[0], t_grid.shape[0]

    # broadcast (N, 1) x (1, T)
    a2 = a[:, None]; e2 = e[:, None]
    n = xp.sqrt(GM_SUN / a2**3)
    M = M0[:, None] + n * (t_grid[None, :] - epoch_jd[:, None])

    # wrap to [-pi, pi] to keep Newton iteration stable
    pi = math.pi
    M = xp.remainder(M + pi, 2 * pi) - pi

    E = kepler_solve(xp, M, e2)
    ca = xp.cos(E); sa = xp.sin(E)
    b = a2 * xp.sqrt(1.0 - e2 * e2)
    x_op = a2 * (ca - e2)
    y_op = b * sa

    cw = xp.cos(w)[:, None];  sw = xp.sin(w)[:, None]
    cO = xp.cos(Om)[:, None]; sO = xp.sin(Om)[:, None]
    ci = xp.cos(inc)[:, None]; si = xp.sin(inc)[:, None]
    R11 = cO*cw - sO*sw*ci
    R12 = -cO*sw - sO*cw*ci
    R21 = sO*cw + cO*sw*ci
    R22 = -sO*sw + cO*cw*ci
    R31 = sw*si
    R32 = cw*si
    x = R11*x_op + R12*y_op
    y = R21*x_op + R22*y_op
    z = R31*x_op + R32*y_op
    return x, y, z, E


# ---------------------------------------------------------------------------
# Elliptic arc length workload (uses the elliptic package on the backend)
# ---------------------------------------------------------------------------

def arc_lengths(xp_name, device, a, e, E):
    """L(E) = a · E(E | e²) — shape (N, T). Perimeter P = 4a·E(e²)."""
    xp = get_xp(xp_name)
    m = (e * e)[:, None]
    m_broad = xp.broadcast_to(m, E.shape)
    _, E_inc, _ = elliptic12(E, m_broad)
    L = a[:, None] * E_inc

    # complete perimeter
    Bc, Dc, _ = ellipticBD(e * e)
    perimeter = 4.0 * a * (Bc + (1.0 - e * e) * Dc)
    return L, perimeter


# ---------------------------------------------------------------------------
# Dataset generators
# ---------------------------------------------------------------------------

def load_or_generate(path: Path | None, n_bodies: int):
    if path and path.exists():
        z = np.load(path, allow_pickle=True)
        take = min(n_bodies, len(z["a"]))
        idx = np.random.default_rng(0).choice(len(z["a"]), size=take, replace=False)
        return (z["a"][idx], z["e"][idx], z["i"][idx],
                z["om"][idx], z["w"][idx], z["ma"][idx], z["epoch"][idx])
    # synthetic main-belt-flavoured sample
    rng = np.random.default_rng(20260421)
    a = rng.uniform(1.5, 4.5, n_bodies)
    e = np.clip(rng.beta(2, 15, n_bodies) * 1.2, 0.0, 0.95)
    inc = np.radians(np.abs(rng.normal(0, 10, n_bodies)))
    Om = rng.uniform(0, 2 * math.pi, n_bodies)
    w = rng.uniform(0, 2 * math.pi, n_bodies)
    M0 = rng.uniform(0, 2 * math.pi, n_bodies)
    epoch = np.full(n_bodies, 2460000.5)
    return a, e, inc, Om, w, M0, epoch


# ---------------------------------------------------------------------------
# Timing
# ---------------------------------------------------------------------------

def time_kernel(fn, xp_name, device, warmup: int = 2, reps: int = 3) -> float:
    for _ in range(warmup):
        out = fn()
    sync(xp_name, device)
    t0 = time.perf_counter()
    for _ in range(reps):
        out = fn()
    sync(xp_name, device)
    return (time.perf_counter() - t0) / reps


def mem_mb(xp_name: str, device: str) -> float:
    if xp_name == "torch" and device == "cuda":
        import torch
        return torch.cuda.max_memory_allocated() / (1024 ** 2)
    return 0.0


def run_one(xp_name, device, dtype, data_np, t_grid_np):
    a_np, e_np, inc_np, Om_np, w_np, M0_np, epoch_np = data_np
    to = lambda x: as_tensor(xp_name, device, x, dtype)
    a, e, inc = to(a_np), to(e_np), to(inc_np)
    Om, w, M0 = to(Om_np), to(w_np), to(M0_np)
    epoch = to(epoch_np); t_grid = to(t_grid_np)

    if xp_name == "torch" and device == "cuda":
        import torch; torch.cuda.reset_peak_memory_stats()

    def kernel_prop():
        return propagate(xp_name, device, a, e, inc, Om, w, M0, epoch, t_grid)

    t_prop = time_kernel(kernel_prop, xp_name, device)
    x, y, z, E = kernel_prop()

    def kernel_arc():
        return arc_lengths(xp_name, device, a, e, E)

    t_arc = time_kernel(kernel_arc, xp_name, device)

    N, T = a_np.shape[0], t_grid_np.shape[0]
    body_epochs = N * T
    return {
        "backend": f"{xp_name}_{device}",
        "dtype": dtype,
        "bodies": N,
        "epochs": T,
        "body_epochs": body_epochs,
        "t_propagate_ms": 1e3 * t_prop,
        "t_arc_length_ms": 1e3 * t_arc,
        "t_total_ms": 1e3 * (t_prop + t_arc),
        "throughput_meps": body_epochs / (t_prop + t_arc) / 1e6,
        "gpu_mem_mb": mem_mb(xp_name, device),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bodies", default="10000,100000,1000000",
                    help="comma-separated body counts")
    ap.add_argument("--epochs", default="128,512,2048",
                    help="comma-separated epoch counts")
    ap.add_argument("--dtype", default="float64", choices=["float32", "float64"])
    ap.add_argument("--npz", type=Path,
                    default=Path(__file__).resolve().parents[2]
                    / "examples/gpu-asteroid-swarm/data/asteroids_full.npz")
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parents[2]
                    / "examples/gpu-asteroid-swarm/data/benchmark_results.json")
    args = ap.parse_args()

    body_list = [int(s) for s in args.bodies.split(",")]
    epoch_list = [int(s) for s in args.epochs.split(",")]
    backends = available_backends()

    device_info = {"platform": platform.platform(),
                   "python": platform.python_version(),
                   "numpy": np.__version__}
    try:
        import torch
        device_info["torch"] = torch.__version__
        if torch.cuda.is_available():
            device_info["cuda_device"] = torch.cuda.get_device_name(0)
            device_info["cuda_version"] = torch.version.cuda
    except ImportError:
        pass

    print(f"Backends: {backends}")
    print(f"Sizes: bodies={body_list}  epochs={epoch_list}  dtype={args.dtype}")

    t_grid_template = np.linspace(0.0, 10 * 365.25, max(epoch_list))
    rows = []
    for N in body_list:
        data = load_or_generate(args.npz if args.npz.exists() else None, N)
        for T in epoch_list:
            t_grid = t_grid_template[:T]
            for xp_name, device in backends:
                try:
                    row = run_one(xp_name, device, args.dtype, data, t_grid)
                    rows.append(row)
                    print(f"  {row['backend']:>14s}  N={N:<8d} T={T:<5d} "
                          f"total={row['t_total_ms']:8.2f} ms  "
                          f"{row['throughput_meps']:6.2f} M body-epochs/s  "
                          f"gpu={row['gpu_mem_mb']:.0f} MB")
                except Exception as err:
                    print(f"  {xp_name}_{device}  N={N} T={T}: {err}")

    # speedup columns relative to numpy at each (N,T)
    base = {(r["bodies"], r["epochs"]): r["t_total_ms"]
            for r in rows if r["backend"] == "numpy_cpu"}
    for r in rows:
        b = base.get((r["bodies"], r["epochs"]))
        r["speedup_vs_numpy"] = round(b / r["t_total_ms"], 3) if b else None

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({
        "device_info": device_info,
        "results": rows,
    }, indent=2))
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
