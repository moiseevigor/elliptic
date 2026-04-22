#!/usr/bin/env python3
"""Validate Keplerian propagation on a sampled subset against JPL Horizons.

For a small set of named asteroids, this queries Horizons for the
heliocentric ecliptic state vector on a given date, propagates the same
objects with our pure-Kepler integrator, and reports position errors.

This is NOT a claim that Keplerian propagation is as accurate as
Horizons — the point is to characterize the drift honestly so the
reader can see how much planetary perturbations matter over the chosen
horizon.

Output: examples/gpu-asteroid-swarm/data/validation_sample.json

Usage:
    python3 scripts/validate_against_horizons.py [--date YYYY-MM-DD] [--targets FILE]
"""
from __future__ import annotations

import argparse
import json
import math
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"

DEFAULT_TARGETS = [
    ("1",    "(1) Ceres"),
    ("2",    "(2) Pallas"),
    ("3",    "(3) Juno"),
    ("4",    "(4) Vesta"),
    ("433",  "(433) Eros"),
    ("1566", "(1566) Icarus"),
    ("4179", "(4179) Toutatis"),
    ("25143","(25143) Itokawa"),
    ("99942","(99942) Apophis"),
    ("101955","(101955) Bennu"),
]


def solve_kepler(M: np.ndarray, e: float, tol: float = 1e-13) -> np.ndarray:
    E = M.copy() if hasattr(M, "copy") else np.asarray(M).copy()
    for _ in range(50):
        dE = (M - E + e * np.sin(E)) / (1 - e * np.cos(E))
        E = E + dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def kepler_state(a, e, i, Om, w, M0, epoch_jd, jd):
    GM_SUN = 0.0002959122082855911  # AU^3 / day^2
    n = math.sqrt(GM_SUN / a**3)
    M = M0 + n * (jd - epoch_jd)
    M = (M + math.pi) % (2 * math.pi) - math.pi
    E = solve_kepler(np.array([M]), e)[0]
    ca, sa = math.cos(E), math.sin(E)
    x_op = a * (ca - e)
    y_op = a * math.sqrt(1 - e * e) * sa
    cw, sw = math.cos(w), math.sin(w)
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(i), math.sin(i)
    R11 = cO * cw - sO * sw * ci
    R12 = -cO * sw - sO * cw * ci
    R21 = sO * cw + cO * sw * ci
    R22 = -sO * sw + cO * cw * ci
    R31 = sw * si
    R32 = cw * si
    x = R11 * x_op + R12 * y_op
    y = R21 * x_op + R22 * y_op
    z = R31 * x_op + R32 * y_op
    return np.array([x, y, z])


def horizons_state(cmd: str, date: str) -> np.ndarray:
    """Return Sun-centred ecliptic J2000 position [AU] on given date."""
    jd_start = date
    # Horizons wants step>=1 — request two days and use the first
    params = {
        "format": "json",
        "COMMAND": f"'{cmd};'",
        "OBJ_DATA": "NO",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "VECTORS",
        "CENTER": "'@sun'",
        "START_TIME": f"'{jd_start}'",
        "STOP_TIME": f"'{jd_start} 00:01'",
        "STEP_SIZE": "'1'",
        "REF_PLANE": "'ECLIPTIC'",
        "REF_SYSTEM": "'J2000'",
        "VEC_TABLE": "'1'",
        "OUT_UNITS": "'AU-D'",
    }
    qs = "&".join(f"{k}={urllib.parse.quote(v)}" for k, v in params.items())
    with urllib.request.urlopen(f"{HORIZONS}?{qs}", timeout=60) as r:
        js = json.load(r)
    txt = js["result"]
    # result block sits between $$SOE / $$EOE markers
    try:
        body = txt.split("$$SOE")[1].split("$$EOE")[0].strip()
    except IndexError:
        raise RuntimeError(f"Horizons returned no ephemeris for {cmd}")
    lines = [ln.strip() for ln in body.splitlines() if ln.strip()]
    # expected: "JD = ..." then "X = ... Y = ... Z = ..."
    xyz_line = next((ln for ln in lines if ln.startswith("X")), None)
    if xyz_line is None:
        raise RuntimeError(f"Could not parse Horizons output for {cmd}")
    parts = xyz_line.replace("=", " ").split()
    x = float(parts[1]); y = float(parts[3]); z = float(parts[5])
    return np.array([x, y, z])


def jd_from_iso(date: str) -> float:
    y, mo, d = (int(p) for p in date.split("-"))
    a = (14 - mo) // 12
    y2 = y + 4800 - a
    m2 = mo + 12 * a - 3
    return (d + (153 * m2 + 2) // 5 + 365 * y2 + y2 // 4
            - y2 // 100 + y2 // 400 - 32045.5)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="2026-01-01")
    ap.add_argument("--npz", type=Path,
                    default=Path(__file__).resolve().parent.parent
                    / "examples/gpu-asteroid-swarm/data/asteroids_full.npz")
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parent.parent
                    / "examples/gpu-asteroid-swarm/data/validation_sample.json")
    args = ap.parse_args()

    z = np.load(args.npz, allow_pickle=True)
    names = np.array([str(x) for x in z["name"]])
    # Leading integer of each full_name, e.g. '     1 Ceres (A801 AA)' -> 1.
    def _leading_int(s: str):
        tok = s.strip().split(None, 1)[0] if s.strip() else ""
        return int(tok) if tok.lstrip("-").isdigit() else -1
    leading_ids = np.array([_leading_int(n) for n in names])
    jd_target = jd_from_iso(args.date)
    print(f"Validation date: {args.date}  (JD {jd_target:.1f})")

    results = []
    for cmd, label in DEFAULT_TARGETS:
        want = int(cmd)
        matches = np.where(leading_ids == want)[0]
        if matches.size == 0:
            print(f"  skip {label}: not in local dataset")
            continue
        k = int(matches[0])
        r_pred = kepler_state(z["a"][k], z["e"][k], z["i"][k], z["om"][k],
                              z["w"][k], z["ma"][k], z["epoch"][k], jd_target)
        try:
            r_ref = horizons_state(cmd, args.date)
        except Exception as err:
            print(f"  skip {label}: Horizons query failed ({err})")
            continue
        err_au = float(np.linalg.norm(r_pred - r_ref))
        dt_years = (jd_target - float(z["epoch"][k])) / 365.25
        print(f"  {label:22s}  |Δr| = {err_au:.3e} AU   ΔT = {dt_years:+.2f} yr")
        results.append({
            "name": label, "cmd": cmd,
            "dt_years": round(dt_years, 3),
            "pred_xyz_au": [round(v, 9) for v in r_pred.tolist()],
            "horizons_xyz_au": [round(v, 9) for v in r_ref.tolist()],
            "err_au": round(err_au, 9),
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({
        "date": args.date, "jd": jd_target,
        "note": "Position error between pure-Kepler propagation (this "
                "repo) and JPL Horizons ephemerides. Horizons includes "
                "full planetary perturbations; differences grow with |ΔT|.",
        "samples": results,
    }, indent=2))
    print(f"Wrote {args.out} ({len(results)} targets)")


if __name__ == "__main__":
    main()
