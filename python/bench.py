"""Benchmark: serial vs multicore vs GPU across key functions.

Usage:
    python3 bench.py             # numpy only
    python3 bench.py --torch     # add torch CPU and CUDA (if available)
    python3 bench.py --all       # numpy + torch + jax
"""
from __future__ import annotations

import argparse
import time
import numpy as np
import elliptic

N_SMALL  = 10_000
N_MEDIUM = 100_000
N_LARGE  = 1_000_000


def _time(fn, *args, reps=3):
    # warm-up
    fn(*args)
    t0 = time.perf_counter()
    for _ in range(reps):
        fn(*args)
    return (time.perf_counter() - t0) / reps


def bench_numpy(N):
    phi = np.random.uniform(0.01, np.pi / 2, N)
    m   = np.random.uniform(0.01, 0.99, N)
    n   = np.random.uniform(0.01, 0.99, N)
    u   = np.random.uniform(0.01, 1.5, N)
    e1  = np.full(N, 2.0); e2 = np.zeros(N); e3 = np.full(N, -2.0)

    rows = []
    for label, fn, args in [
        ("elliptic12",   elliptic.elliptic12,   (phi, m)),
        ("elliptic3",    elliptic.elliptic3,     (phi, m, n)),
        ("ellipj",       elliptic.ellipj,        (u, m)),
        ("ellipticBDJ",  elliptic.ellipticBDJ,   (phi, m, n)),
        ("ellipticBD",   elliptic.ellipticBD,    (m,)),
        ("carlsonRF",    elliptic.carlsonRF,     (m, 1-m, np.ones(N))),
        ("weierstrassP", elliptic.weierstrassP,  (u, e1, e2, e3)),
    ]:
        t = _time(fn, *args)
        rows.append((label, t * 1000))

    return rows


def bench_torch(N):
    import torch
    phi = torch.from_numpy(np.random.uniform(0.01, np.pi / 2, N))
    m   = torch.from_numpy(np.random.uniform(0.01, 0.99, N))
    n   = torch.from_numpy(np.random.uniform(0.01, 0.99, N))

    rows = []
    for label, fn, args in [
        ("elliptic12 [torch-cpu]", elliptic.elliptic12, (phi, m)),
        ("ellipticBDJ [torch-cpu]", elliptic.ellipticBDJ, (phi, m, n)),
    ]:
        t = _time(fn, *args)
        rows.append((label, t * 1000))

    if torch.cuda.is_available():
        phi_c = phi.cuda(); m_c = m.cuda(); n_c = n.cuda()
        # warm-up GPU
        elliptic.elliptic12(phi_c, m_c)
        for label, fn, args in [
            ("elliptic12 [cuda]",   elliptic.elliptic12,  (phi_c, m_c)),
            ("ellipticBDJ [cuda]",  elliptic.ellipticBDJ, (phi_c, m_c, n_c)),
        ]:
            t = _time(fn, *args)
            rows.append((label, t * 1000))
    else:
        rows.append(("CUDA", float("nan")))
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--torch", action="store_true")
    parser.add_argument("--all",   action="store_true")
    args = parser.parse_args()

    print(f"{'Function':<30} {'N':>10}  {'ms/call':>10}")
    print("-" * 56)

    for N, tag in [(N_SMALL, "10k"), (N_MEDIUM, "100k"), (N_LARGE, "1M")]:
        rows = bench_numpy(N)
        for label, ms in rows:
            print(f"{label:<30} {tag:>10}  {ms:>10.2f}")
        print()

    if args.torch or args.all:
        print("--- PyTorch ---")
        rows = bench_torch(N_MEDIUM)
        for label, ms in rows:
            print(f"{label:<45}  {ms:>10.2f} ms")


if __name__ == "__main__":
    main()
