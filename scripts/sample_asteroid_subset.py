#!/usr/bin/env python3
"""Build a lightweight browser-friendly asteroid subset.

Reads data/asteroids_full.npz and emits:
    data/asteroids_subset.json       (plain; loaded by the browser example)
    data/asteroids_subset.json.gz    (optional; --gzip)

The subset is stratified across orbit classes and covers the full
semi-major-axis range so the 3D visualization looks representative.

Default size: 10 000 bodies — small enough to render smoothly with
instanced Three.js meshes while preserving the visible structure of
the main belt, Hildas, Trojans, etc.

Usage:
    python3 scripts/sample_asteroid_subset.py [--n 10000] [--in FILE] [--out FILE]
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path

import numpy as np


def stratified_sample(data: dict, n_target: int, seed: int = 20260421) -> np.ndarray:
    rng = np.random.default_rng(seed)
    N = len(data["a"])
    if n_target >= N:
        return np.arange(N)
    classes, inv = np.unique(data["orbit_class"], return_inverse=True)
    per_class = max(1, n_target // len(classes))
    idx = []
    for k in range(len(classes)):
        pool = np.where(inv == k)[0]
        take = min(len(pool), per_class)
        idx.append(rng.choice(pool, size=take, replace=False))
    idx = np.concatenate(idx)
    if len(idx) > n_target:
        idx = rng.choice(idx, size=n_target, replace=False)
    elif len(idx) < n_target:
        pool = np.setdiff1d(np.arange(N), idx, assume_unique=False)
        extra = rng.choice(pool, size=n_target - len(idx), replace=False)
        idx = np.concatenate([idx, extra])
    return idx


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10_000)
    ap.add_argument("--in", dest="infile", type=Path,
                    default=Path(__file__).resolve().parent.parent
                    / "examples/gpu-asteroid-swarm/data/asteroids_full.npz")
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parent.parent
                    / "examples/gpu-asteroid-swarm/data/asteroids_subset.json")
    ap.add_argument("--gzip", action="store_true", help="also write .gz next to output")
    args = ap.parse_args()

    z = np.load(args.infile, allow_pickle=True)
    data = {k: z[k] for k in z.files}
    idx = stratified_sample(data, args.n)

    payload = {
        "count": int(len(idx)),
        "epoch_jd": float(np.median(data["epoch"][idx])),
        "units": {"a": "AU", "i": "rad", "Om": "rad", "w": "rad", "M0": "rad"},
        "fields": ["a", "e", "i", "Om", "w", "M0", "orbit_class"],
        "bodies": [
            {
                "a":   round(float(data["a"][k]),   6),
                "e":   round(float(data["e"][k]),   6),
                "i":   round(float(data["i"][k]),   6),
                "Om":  round(float(data["om"][k]),  6),
                "w":   round(float(data["w"][k]),   6),
                "M0":  round(float(data["ma"][k]),  6),
                "c":   str(data["orbit_class"][k]),
            } for k in idx
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(payload, f, separators=(",", ":"))
    size_kb = args.out.stat().st_size / 1024
    print(f"Wrote {args.out} — {len(idx)} bodies, {size_kb:.1f} KB")
    if args.gzip:
        gz_path = args.out.with_suffix(args.out.suffix + ".gz")
        with gzip.open(gz_path, "wt") as f:
            json.dump(payload, f, separators=(",", ":"))
        print(f"Wrote {gz_path} — {gz_path.stat().st_size/1024:.1f} KB")


if __name__ == "__main__":
    main()
