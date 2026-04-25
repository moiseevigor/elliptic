// ═══════════════════════════════════════════════════════════════════════
// Resonance Stability of the Asteroid Belt
// Semi-analytical Poincaré stability index via K(k²)
//
// Theory: Murray & Dermott (1999) Chapters 7–8; Wisdom (1982)
// ═══════════════════════════════════════════════════════════════════════

// ── Constants ──────────────────────────────────────────────────────────
const GM_SUN = 4 * Math.PI * Math.PI;   // [AU³ yr⁻²]
const A_J    = 5.2026;                   // Jupiter semi-major axis [AU]
const MU_J   = 9.5474e-4;               // M_Jupiter / M_Sun

// ── Elliptic integral K(m) via AGM ─────────────────────────────────────
// Complete elliptic integral of the first kind K(m) = F(π/2 | m).
// Used to compute the exact libration period T_lib = 4K(k²)/ω₀.
function ellipticK(m) {
  if (m >= 1) return Infinity;
  if (m <= 0) return Math.PI / 2;
  let a = 1, b = Math.sqrt(1 - m);
  for (let i = 0; i < 32; i++) {
    const an = (a + b) / 2;
    b = Math.sqrt(a * b);
    a = an;
    if (Math.abs(a - b) < 1e-15 * a) break;
  }
  return Math.PI / (2 * a);
}

// ── Laplace coefficient b_{1/2}^{(j)}(α) by trapezoidal quadrature ────
// (Murray & Dermott, Appendix A, eq. A.1)
function lapB(j, alpha, n = 1024) {
  let sum = 0;
  const h = Math.PI / n;
  for (let k = 0; k <= n; k++) {
    const psi = k * h;
    const d   = 1 - 2 * alpha * Math.cos(psi) + alpha * alpha;
    const f   = Math.cos(j * psi) / Math.sqrt(d);
    sum += (k === 0 || k === n) ? 0.5 * f : f;
  }
  return (2 / Math.PI) * h * sum;
}

function dLapB(j, alpha) {
  const h = 1e-5;
  return (lapB(j, alpha + h) - lapB(j, alpha - h)) / (2 * h);
}

// ── Resonant disturbing-function coefficient |f_d| ─────────────────────
// For a first-order p:q resonance (p − q = 1), asteroid inside Jupiter.
// Dominant direct term: Murray & Dermott Table B.6 form.
// Resonant angle: σ = pλ − qλ_J − ϖ
function resonantFd(p, q) {
  const alpha = Math.pow(q / p, 2 / 3);   // a_res / a_J
  const j = q;                              // harmonic index in disturbing function
  const b  = lapB(j, alpha);
  const db = dLapB(j, alpha);
  let f = -((2 * j + 1) * b + alpha * db) / 2;
  if (j === 1) f += alpha / 2;             // indirect term for j=1 (2:1 and 3:2)
  return Math.abs(f);
}

// ── Resonance geometry ────────────────────────────────────────────────
function aRes(p, q) { return A_J * Math.pow(q / p, 2 / 3); }
function nRes(p, q) { return Math.sqrt(GM_SUN / Math.pow(aRes(p, q), 3)); }

// Separatrix half-width in a [AU] (Murray & Dermott eq. 8.76)
function halfwidth(p, q, fd, e) {
  const ar = aRes(p, q);
  return ar * Math.sqrt(16 * MU_J * fd * Math.max(e, 1e-4) / (3 * p * p));
}

// Small-oscillation frequency ω₀ [rad yr⁻¹]
function omega0(p, q, fd, e) {
  return nRes(p, q) * Math.sqrt(1.5 * p * p * MU_J * fd * Math.max(e, 1e-4));
}

// ── Resonance table (first-order only) ───────────────────────────────
// Includes the main Kirkwood gaps AND the stable Hilda 3:2 resonance.
// Second-order resonances (5:2, 7:3) require a separate treatment (width ∝ e
// rather than √e) and are excluded from this first-order survey.
const RESONANCES = [
  { name: '4:1', p: 4, q: 1 },
  { name: '3:1', p: 3, q: 1 },
  { name: '2:1', p: 2, q: 1 },
  { name: '3:2', p: 3, q: 2 },   // Hilda cluster — stable resonance
  { name: '4:3', p: 4, q: 3 },
].map(r => ({ ...r, fd: resonantFd(r.p, r.q), a: aRes(r.p, r.q) }));

// ── Per-body stability computation ────────────────────────────────────
function stability(a, e) {
  let best_k2   = Infinity;
  let best_om0  = 1;
  let best_name = '—';

  for (const r of RESONANCES) {
    const dw = halfwidth(r.p, r.q, r.fd, e);
    const k2 = Math.pow((a - r.a) / (dw + 1e-12), 2);
    if (k2 < best_k2) {
      best_k2   = k2;
      best_om0  = omega0(r.p, r.q, r.fd, e);
      best_name = r.name;
    }
  }

  // Stability index: +1 = resonance centre, 0 = separatrix, <0 = circulating
  const S = 1 - best_k2;
  // Exact libration period via K(k²); infinite at and beyond separatrix
  const K_val  = best_k2 < 1 ? ellipticK(best_k2) : Infinity;
  const T_lib  = best_k2 < 1 ? 4 * K_val / best_om0 : Infinity;
  const T_orb  = 2 * Math.PI * Math.sqrt(Math.pow(a, 3) / GM_SUN);

  return { S, k_sq: best_k2, K: K_val, T_lib, T_ratio: T_lib / T_orb, res: best_name };
}

// ── Figure helpers ─────────────────────────────────────────────────────
function figW(el) {
  return el.getBoundingClientRect().width || el.parentElement.getBoundingClientRect().width || 600;
}

function clipGroup(g, m, W, H) {
  const id = 'clip-' + Math.random().toString(36).slice(2, 8);
  g.append('clipPath').attr('id', id)
    .append('rect').attr('x', m.l).attr('y', m.t)
    .attr('width', W - m.l - m.r).attr('height', H - m.t - m.b);
  return g.append('g').attr('clip-path', `url(#${id})`);
}

// ══════════════════════════════════════════════════════════════════════
// §1  Pendulum phase portrait
// ══════════════════════════════════════════════════════════════════════
function drawPhasePlot() {
  const svg = document.getElementById('fig-phase');
  const W = figW(svg), H = 320;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const g = d3.select(svg); g.selectAll('*').remove();
  const m = { t: 18, b: 40, l: 52, r: 16 };

  const xS = d3.scaleLinear().domain([-Math.PI, Math.PI]).range([m.l, W - m.r]);
  const yS = d3.scaleLinear().domain([-2.6, 2.6]).range([H - m.b, m.t]);

  // Grid lines
  const plot = clipGroup(g, m, W, H);
  [-Math.PI, -Math.PI/2, 0, Math.PI/2, Math.PI].forEach(x => {
    plot.append('line').attr('x1', xS(x)).attr('x2', xS(x)).attr('y1', m.t).attr('y2', H - m.b)
      .attr('stroke', '#f0f0f0').attr('stroke-width', 1);
  });

  // Energy contours H = ½ṡ² − cos(σ) = const
  // At k²: H = 2k² − 1  ⟹  ṡ = ±√(2(H + cos σ)) = ±√(2(2k² − 1 + cos σ))
  const colors = d3.schemeTableau10;
  const k_vals = [0.0, 0.25, 0.5, 0.75, 0.90, 0.97];   // k (not k²)
  k_vals.forEach((k, ci) => {
    const H_energy = 2 * k * k - 1;
    const color = ci < 5 ? d3.interpolateBlues(0.3 + ci * 0.14) : '#e65100';
    for (const sign of [1, -1]) {
      const pts = [];
      for (let i = 0; i <= 240; i++) {
        const sigma = -Math.PI + (i / 240) * 2 * Math.PI;
        const sdot_sq = 2 * (H_energy + Math.cos(sigma));
        if (sdot_sq < 0) continue;
        pts.push([sigma, sign * Math.sqrt(sdot_sq)]);
      }
      if (pts.length > 1) {
        plot.append('path')
          .attr('d', d3.line().x(p => xS(p[0])).y(p => yS(p[1]))(pts))
          .attr('fill', 'none').attr('stroke', color).attr('stroke-width', k >= 0.9 ? 2 : 1.5)
          .attr('stroke-opacity', 0.85);
      }
    }
  });

  // Separatrix (k = 1, H = 1): ṡ = ±√(2(1+cos σ)) = ±2|cos(σ/2)|
  for (const sign of [1, -1]) {
    const pts = [];
    for (let i = 0; i <= 400; i++) {
      const sigma = -Math.PI + (i / 400) * 2 * Math.PI;
      pts.push([sigma, sign * 2 * Math.abs(Math.cos(sigma / 2))]);
    }
    plot.append('path')
      .attr('d', d3.line().x(p => xS(p[0])).y(p => yS(p[1]))(pts))
      .attr('fill', 'none').attr('stroke', '#b71c1c').attr('stroke-width', 2.2)
      .attr('stroke-dasharray', '6,3');
  }

  // Circulating orbits (k > 1)
  [1.1, 1.3].forEach(k => {
    const H_energy = 2 * k * k - 1;
    for (const sign of [1, -1]) {
      const pts = [];
      for (let i = 0; i <= 240; i++) {
        const sigma = -Math.PI + (i / 240) * 2 * Math.PI;
        const v = 2 * (H_energy + Math.cos(sigma));
        pts.push([sigma, sign * Math.sqrt(Math.max(v, 0))]);
      }
      plot.append('path')
        .attr('d', d3.line().x(p => xS(p[0])).y(p => yS(p[1]))(pts))
        .attr('fill', 'none').attr('stroke', '#9e9e9e').attr('stroke-width', 1.5)
        .attr('stroke-dasharray', '3,3').attr('stroke-opacity', 0.6);
    }
  });

  // Axes
  g.append('g').attr('transform', `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xS).tickValues([-Math.PI, -Math.PI / 2, 0, Math.PI / 2, Math.PI])
      .tickFormat(d => d === 0 ? '0' : (d / Math.PI).toFixed(1) + 'π'))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  g.append('g').attr('transform', `translate(${m.l},0)`)
    .call(d3.axisLeft(yS).ticks(5))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  g.append('text').attr('x', (W + m.l) / 2).attr('y', H - 10)
    .attr('font-family', 'Source Sans 3').attr('font-size', 12).attr('text-anchor', 'middle').attr('fill', '#555')
    .text('resonance angle σ');
  g.append('text').attr('transform', `translate(16,${(H + m.t) / 2})rotate(-90)`)
    .attr('font-family', 'Source Sans 3').attr('font-size', 12).attr('text-anchor', 'middle').attr('fill', '#555')
    .text('dσ/dt  [normalised]');

  // Labels
  g.append('text').attr('x', xS(0) + 4).attr('y', yS(0.18))
    .attr('font-family', 'Source Sans 3').attr('font-size', 11).attr('fill', '#1565c0')
    .text('libration  (k < 1, S > 0)');
  g.append('text').attr('x', xS(-Math.PI) + 6).attr('y', yS(2.4))
    .attr('font-family', 'Source Sans 3').attr('font-size', 11).attr('fill', '#888')
    .text('circulation  (k > 1)');
  g.append('text').attr('x', xS(Math.PI / 2) - 2).attr('y', yS(1.05))
    .attr('font-family', 'Source Sans 3').attr('font-size', 11).attr('fill', '#b71c1c')
    .text('separatrix  k = 1, S = 0');
}

// ══════════════════════════════════════════════════════════════════════
// §2  K(k²) divergence — libration period vs amplitude
// ══════════════════════════════════════════════════════════════════════
function drawKPlot() {
  const svg = document.getElementById('fig-kplot');
  const W = figW(svg), H = 260;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const g = d3.select(svg); g.selectAll('*').remove();
  const m = { t: 16, b: 40, l: 56, r: 20 };

  const kPts = d3.range(0, 0.998, 0.002);
  const T_norm = kPts.map(k => ({
    k,
    ratio: (2 / Math.PI) * ellipticK(k * k)   // T_lib / T_0 = (2/π)K(k²)
  }));

  const xK = d3.scaleLinear().domain([0, 1]).range([m.l, W - m.r]);
  const yK = d3.scaleLog().domain([1, 16]).range([H - m.b, m.t]);

  const plot = clipGroup(g, m, W, H);

  // Horizontal reference at T/T_0 = 1
  plot.append('line').attr('x1', m.l).attr('x2', W - m.r)
    .attr('y1', yK(1)).attr('y2', yK(1))
    .attr('stroke', '#e0e0e0').attr('stroke-width', 1);

  // K(k²) curve
  plot.append('path')
    .attr('d', d3.line().x(d => xK(d.k)).y(d => yK(d.ratio))(T_norm))
    .attr('fill', 'none').attr('stroke', '#1565c0').attr('stroke-width', 2.2);

  // Shade the chaotic zone near k = 1
  plot.append('rect')
    .attr('x', xK(0.9)).attr('width', xK(1) - xK(0.9))
    .attr('y', m.t).attr('height', H - m.t - m.b)
    .attr('fill', '#b71c1c').attr('fill-opacity', 0.07);

  g.append('text').attr('x', xK(0.95)).attr('y', m.t + 14)
    .attr('text-anchor', 'middle').attr('font-family', 'Source Sans 3')
    .attr('font-size', 10).attr('fill', '#b71c1c').text('chaotic layer');

  // Axes
  g.append('g').attr('transform', `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xK).ticks(5).tickFormat(d => d === 1 ? '1 (sep.)' : d.toFixed(1)))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  g.append('g').attr('transform', `translate(${m.l},0)`)
    .call(d3.axisLeft(yK).tickValues([1, 2, 4, 8]).tickFormat(d => d + '×'))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  g.append('text').attr('x', (W + m.l) / 2).attr('y', H - 10)
    .attr('font-family', 'Source Sans 3').attr('font-size', 12).attr('text-anchor', 'middle').attr('fill', '#555')
    .text('amplitude parameter k');
  g.append('text').attr('transform', `translate(16,${(H + m.t) / 2})rotate(-90)`)
    .attr('font-family', 'Source Sans 3').attr('font-size', 12).attr('text-anchor', 'middle').attr('fill', '#555')
    .text('T_lib / T_0  =  (2/π) K(k²)');
}

// ══════════════════════════════════════════════════════════════════════
// §3  Stability map — (a, e) coloured by S
// ══════════════════════════════════════════════════════════════════════
function drawStabilityMap(bodies) {
  const svg = document.getElementById('fig-stability');
  const W = figW(svg), H = 380;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const g = d3.select(svg); g.selectAll('*').remove();
  const m = { t: 20, b: 40, l: 52, r: 20 };

  const a_arr = bodies.map(b => b.a);
  const e_arr = bodies.map(b => b.e);

  // Focus on main belt + inner solar system
  const aMin = 1.6, aMax = 5.5;
  const mask = bodies.map((b, i) => b.a >= aMin && b.a <= aMax && b.e < 0.5);
  const bShow = bodies.filter((_, i) => mask[i]);

  const xA = d3.scaleLinear().domain([aMin, aMax]).range([m.l, W - m.r]);
  const yE = d3.scaleLinear().domain([0, 0.48]).range([H - m.b, m.t]);

  // Compute stability for each body
  const stab = bShow.map(b => stability(b.a, b.e));

  // Color scale: S ∈ (−∞, 1]
  // circulating far  (S < −1):  steel blue
  // circulating near (−1 ≤ S < 0):  light orange → orange
  // separatrix       (S ≈ 0):    deep red
  // librating        (S > 0):    green
  function stabColor(S, k_sq) {
    if (k_sq >= 1) {
      // circulating — map |S| to a grey-blue scale
      const t = Math.max(0, Math.min(1, 1 - Math.exp(-k_sq + 1)));
      return d3.interpolate('#4fc3f7', '#37474f')(t);
    }
    // librating: map S from 0 (red) to 1 (green)
    return d3.interpolateRdYlGn((S + 0.1) / 1.1);
  }

  const plot = clipGroup(g, m, W, H);

  // Resonance gap lines
  RESONANCES.forEach(r => {
    if (r.a >= aMin && r.a <= aMax) {
      plot.append('line')
        .attr('x1', xA(r.a)).attr('x2', xA(r.a))
        .attr('y1', m.t).attr('y2', H - m.b)
        .attr('stroke', '#ddd').attr('stroke-width', 1).attr('stroke-dasharray', '3,3');
      g.append('text').attr('x', xA(r.a)).attr('y', m.t - 5)
        .attr('text-anchor', 'middle').attr('font-family', 'Source Sans 3')
        .attr('font-size', 9.5).attr('fill', '#aaa').text(r.name);
    }
  });

  // Points
  plot.append('g').selectAll('circle').data(bShow).join('circle')
    .attr('cx', b => xA(b.a))
    .attr('cy', b => yE(b.e))
    .attr('r', 1.6)
    .attr('fill', (b, i) => stabColor(stab[i].S, stab[i].k_sq))
    .attr('fill-opacity', 0.8);

  // Axes
  g.append('g').attr('transform', `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xA).ticks(8).tickFormat(d => d + ' AU'))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  g.append('g').attr('transform', `translate(${m.l},0)`)
    .call(d3.axisLeft(yE).ticks(5))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  g.append('text').attr('x', (W + m.l) / 2).attr('y', H - 10)
    .attr('font-family', 'Source Sans 3').attr('font-size', 12).attr('text-anchor', 'middle').attr('fill', '#555')
    .text('semi-major axis a [AU]');
  g.append('text').attr('transform', `translate(16,${(H + m.t) / 2})rotate(-90)`)
    .attr('font-family', 'Source Sans 3').attr('font-size', 12).attr('text-anchor', 'middle').attr('fill', '#555')
    .text('eccentricity e');

  // Colour legend
  const legX = W - m.r - 130, legY = H - m.b - 60;
  const legW = 120, legH = 10;
  const defs = g.append('defs');
  const grad = defs.append('linearGradient').attr('id', 'stab-grad');
  for (let i = 0; i <= 10; i++) {
    grad.append('stop').attr('offset', `${i * 10}%`)
      .attr('stop-color', d3.interpolateRdYlGn(i / 10));
  }
  g.append('rect').attr('x', legX).attr('y', legY).attr('width', legW).attr('height', legH)
    .attr('fill', 'url(#stab-grad)').attr('rx', 2);
  g.append('text').attr('x', legX).attr('y', legY - 4)
    .attr('font-family', 'Source Sans 3').attr('font-size', 10).attr('fill', '#888')
    .text('separatrix');
  g.append('text').attr('x', legX + legW).attr('y', legY - 4)
    .attr('text-anchor', 'end').attr('font-family', 'Source Sans 3').attr('font-size', 10).attr('fill', '#888')
    .text('libration centre');
}

// ══════════════════════════════════════════════════════════════════════
// §4  Libration period map — T_lib / T_orb coloured scatter
// ══════════════════════════════════════════════════════════════════════
function drawPeriodMap(bodies) {
  const svg = document.getElementById('fig-tlib');
  const W = figW(svg), H = 300;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const g = d3.select(svg); g.selectAll('*').remove();
  const m = { t: 16, b: 40, l: 52, r: 20 };

  const aMin = 1.8, aMax = 4.4;
  const bLib = bodies.filter(b => {
    const s = stability(b.a, b.e);
    return b.a >= aMin && b.a <= aMax && b.e > 0.02 && b.e < 0.45 && s.k_sq < 1;
  });
  const stab = bLib.map(b => stability(b.a, b.e));

  const xA = d3.scaleLinear().domain([aMin, aMax]).range([m.l, W - m.r]);
  const yR = d3.scaleLog().domain([1, 500]).range([H - m.b, m.t]);

  const plot = clipGroup(g, m, W, H);

  RESONANCES.forEach(r => {
    if (r.a >= aMin && r.a <= aMax) {
      plot.append('line')
        .attr('x1', xA(r.a)).attr('x2', xA(r.a))
        .attr('y1', m.t).attr('y2', H - m.b)
        .attr('stroke', '#f0f0f0').attr('stroke-width', 1);
    }
  });

  const colorT = d3.scaleSequentialLog(d3.interpolatePlasma).domain([1, 500]);
  plot.append('g').selectAll('circle').data(bLib).join('circle')
    .attr('cx', (b, i) => xA(b.a))
    .attr('cy', (b, i) => yR(Math.min(stab[i].T_ratio, 499)))
    .attr('r', 2).attr('fill', (b, i) => colorT(Math.min(stab[i].T_ratio, 499)))
    .attr('fill-opacity', 0.75);

  g.append('g').attr('transform', `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xA).ticks(7).tickFormat(d => d + ' AU'))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  g.append('g').attr('transform', `translate(${m.l},0)`)
    .call(d3.axisLeft(yR).tickValues([1, 5, 20, 100, 500]).tickFormat(d => d + '×'))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  g.append('text').attr('x', (W + m.l) / 2).attr('y', H - 10)
    .attr('font-family', 'Source Sans 3').attr('font-size', 12).attr('text-anchor', 'middle').attr('fill', '#555')
    .text('semi-major axis a [AU]');
  g.append('text').attr('transform', `translate(16,${(H + m.t) / 2})rotate(-90)`)
    .attr('font-family', 'Source Sans 3').attr('font-size', 12).attr('text-anchor', 'middle').attr('fill', '#555')
    .text('T_lib / T_orb');
}

// ── Bootstrap ──────────────────────────────────────────────────────────
function loadJSON(url) {
  return fetch(url).then(r => r.json());
}

window.addEventListener('load', () => {
  // Static figures that don't need data
  drawPhasePlot();
  drawKPlot();
  window.addEventListener('resize', () => { drawPhasePlot(); drawKPlot(); });

  // Data-driven figures
  loadJSON('../gpu-asteroid-swarm/data/asteroids_subset.json').then(sub => {
    const bodies = sub.bodies;
    drawStabilityMap(bodies);
    drawPeriodMap(bodies);
    window.addEventListener('resize', () => {
      drawStabilityMap(bodies);
      drawPeriodMap(bodies);
    });
  }).catch(err => {
    console.warn('Could not load subset data:', err);
  });
});
