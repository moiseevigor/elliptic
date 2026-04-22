# GPU Asteroid Swarm

Flagship GPU example: massive batched Keplerian propagation of real
small-body orbital elements, with elliptic-integral arc lengths from the
`elliptic` package. Published under
<https://moiseevigor.github.io/elliptic/examples/gpu-asteroid-swarm/>.

## Pipeline

```
raw MPC / SBDB catalog
         │
         ▼
scripts/prepare_asteroid_dataset.py  →  data/asteroids_full.npz (+ meta)
         │
         ▼
scripts/sample_asteroid_subset.py    →  data/asteroids_subset.json
         │
         ▼
python/examples/gpu_asteroid_benchmark.py   →  data/benchmark_results.json
scripts/validate_against_horizons.py        →  data/validation_sample.json
matlab/examples/gpu_asteroid_benchmark.m    →  data/benchmark_results_matlab.json
         │
         ▼
index.html + app.js + style.css  (loads the JSON artefacts at runtime)
```

## Regenerate everything

```bash
# 1. full catalog from JPL SBDB (optionally --max to cap)
python3 scripts/prepare_asteroid_dataset.py

# 2. compact browser-friendly subset
python3 scripts/sample_asteroid_subset.py --n 10000

# 3. GPU benchmark (NumPy baseline + PyTorch CPU/CUDA)
python3 python/examples/gpu_asteroid_benchmark.py

# 4. Horizons validation (small sampled set of named asteroids)
python3 scripts/validate_against_horizons.py --date 2026-01-01

# 5. MATLAB benchmark (in MATLAB / Octave)
#    >> run('matlab/examples/gpu_asteroid_benchmark.m')
```

## Preview locally

```bash
cd examples/gpu-asteroid-swarm
python3 -m http.server 8080
# then open http://localhost:8080/
```

The page degrades gracefully:

- If `data/asteroids_subset.json` is missing, a synthetic main-belt
  subset is generated in-browser.
- `benchmark_results.json` reflects whatever hardware the benchmark
  script was last run on — re-run it locally to overwrite with your
  own numbers.

## Scientific scope — honest version

- This is a **Keplerian / osculating-ellipse** throughput demo.
  Gravity between planets is ignored.
- Elliptic integrals are used for exact geometric quantities on the
  orbital ellipses (arc length, perimeter).
- For precision ephemerides use
  [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/) or SPICE.
- The Validation section in the page quantifies how much Horizons
  differs from pure-Kepler across various objects and time offsets.

## Data attributions

- Orbital elements:
  [JPL Small-Body Database Query API](https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html),
  [Minor Planet Center MPCORB](https://www.minorplanetcenter.net/iau/MPCORB.html).
- Ephemeris reference:
  [JPL Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html).
- Planetary mean elements: JPL J2000 approximation
  (Standish & Williams 1992).
