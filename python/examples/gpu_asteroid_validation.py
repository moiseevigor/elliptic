#!/usr/bin/env python3
"""Short-horizon validation: pure-Kepler vs JPL Horizons on named asteroids.

This is a thin wrapper around scripts/validate_against_horizons.py that
reads the resulting validation_sample.json and prints a summary table.
Useful for quick sanity-checking after regenerating the data.
"""
from __future__ import annotations

import json
from pathlib import Path

DEFAULT = (Path(__file__).resolve().parents[2]
           / "examples/gpu-asteroid-swarm/data/validation_sample.json")


def main(path: Path = DEFAULT) -> None:
    if not path.exists():
        print(f"{path} not found — run scripts/validate_against_horizons.py first")
        return
    payload = json.loads(path.read_text())
    print(f"Validation date: {payload['date']}  (JD {payload['jd']:.1f})")
    print(f"{'Object':26s} {'ΔT [yr]':>9s} {'|Δr| [AU]':>12s}")
    print("-" * 50)
    for s in payload["samples"]:
        print(f"{s['name']:26s} {s['dt_years']:+9.2f} {s['err_au']:12.4e}")


if __name__ == "__main__":
    main()
