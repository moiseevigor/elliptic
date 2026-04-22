#!/usr/bin/env python3
"""Build a compact asteroid orbital-element dataset for the GPU benchmark.

Primary source: JPL Small-Body Database (SBDB) Query API
    https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html

Produces:
    data/asteroids_full.npz      — all retained bodies (numpy arrays)
    data/asteroids_full.parquet  — same data in Parquet (optional, if pyarrow)
    data/asteroids_meta.json     — bookkeeping (epoch, count, filters)

Only elliptic (e < 1) bodies are kept. Elements stored:
    a [AU], e, i [rad], Om [rad], w [rad], M0 [rad], epoch_jd, H, orbit_class

Usage:
    python3 scripts/prepare_asteroid_dataset.py [--max N] [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import math
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
FIELDS = ["full_name", "a", "e", "i", "om", "w", "ma", "epoch", "H", "class"]


def fetch_sbdb(max_rows: int | None = None) -> dict:
    params = {
        "fields": ",".join(FIELDS),
        "sb-kind": "a",          # asteroids only
        "sb-cdata": '{"AND":["e|LT|1","a|GT|0"]}',
        "full-prec": "true",
    }
    if max_rows:
        params["limit"] = str(max_rows)
    qs = "&".join(f"{k}={urllib.parse.quote(v)}" for k, v in params.items())
    url = f"{SBDB_URL}?{qs}"
    print(f"GET {url}")
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)


def parse(js: dict) -> dict[str, np.ndarray]:
    fields = js["fields"]
    idx = {f: i for i, f in enumerate(fields)}
    rows = js["data"]
    n = len(rows)
    name = np.empty(n, dtype=object)
    a = np.empty(n); e = np.empty(n); inc = np.empty(n)
    om = np.empty(n); w = np.empty(n); ma = np.empty(n)
    ep = np.empty(n); H = np.empty(n); cls = np.empty(n, dtype=object)

    for k, row in enumerate(rows):
        try:
            name[k] = row[idx["full_name"]]
            a[k]    = float(row[idx["a"]])
            e[k]    = float(row[idx["e"]])
            inc[k]  = math.radians(float(row[idx["i"]]))
            om[k]   = math.radians(float(row[idx["om"]]))
            w[k]    = math.radians(float(row[idx["w"]]))
            ma[k]   = math.radians(float(row[idx["ma"]]))
            ep[k]   = float(row[idx["epoch"]])
            H[k]    = float(row[idx["H"]]) if row[idx["H"]] is not None else np.nan
            cls[k]  = row[idx["class"]] or "UNK"
        except (TypeError, ValueError):
            a[k] = np.nan

    keep = np.isfinite(a) & np.isfinite(e) & (e < 1.0) & (a > 0.0)
    print(f"Parsed {n}, kept {keep.sum()} valid elliptic bodies")
    return {
        "name": name[keep], "a": a[keep], "e": e[keep], "i": inc[keep],
        "om": om[keep], "w": w[keep], "ma": ma[keep],
        "epoch": ep[keep], "H": H[keep], "class": cls[keep],
    }


def save(data: dict, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        outdir / "asteroids_full.npz",
        name=data["name"], a=data["a"], e=data["e"], i=data["i"],
        om=data["om"], w=data["w"], ma=data["ma"],
        epoch=data["epoch"], H=data["H"], orbit_class=data["class"],
    )
    meta = {
        "count": int(len(data["a"])),
        "epoch_median_jd": float(np.median(data["epoch"])),
        "a_stats": {"min": float(data["a"].min()), "max": float(data["a"].max()),
                    "median": float(np.median(data["a"]))},
        "e_stats": {"min": float(data["e"].min()), "max": float(data["e"].max()),
                    "median": float(np.median(data["e"]))},
        "source": "JPL SBDB Query API (sb-kind=a, e<1)",
        "fields": FIELDS,
    }
    (outdir / "asteroids_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"Wrote {outdir/'asteroids_full.npz'} ({len(data['a'])} bodies)")

    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
        tbl = pa.table({k: v for k, v in data.items()})
        pq.write_table(tbl, outdir / "asteroids_full.parquet", compression="zstd")
        print(f"Wrote {outdir/'asteroids_full.parquet'}")
    except ImportError:
        print("pyarrow not installed — skipping Parquet output")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max", type=int, default=None,
                    help="cap number of rows fetched (default: all)")
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parent.parent
                    / "examples/gpu-asteroid-swarm/data")
    args = ap.parse_args()
    js = fetch_sbdb(args.max)
    data = parse(js)
    save(data, args.out)


if __name__ == "__main__":
    main()
