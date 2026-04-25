// ═══════════════════════════════════════════════════════════════════════
// Dubins Car, Euler's Elastica & the Visual Cortex
//
// The curvature κ(s) of an optimal curve in SE(2) satisfies the
// nonlinear pendulum equation.  Solutions are Jacobi elliptic functions;
// their spatial period is governed by K(k²).
//
// References:
//   Moiseev & Sachkov (2010) ESAIM:COCV 16:380–399  arXiv:0807.4731
//   Petitot (2003) J. Physiology (Paris) 97:265–309
//   Citti & Sarti (2006) J. Math. Imaging & Vision 24:307–326
// ═══════════════════════════════════════════════════════════════════════

// ── K(m) via arithmetic-geometric mean ─────────────────────────────────
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

// ── Jacobi elliptic functions via descending Landen transform ───────────
// Returns {sn, cn, dn, am} for u and m ∈ [0, 1].
function ellipj(u, m) {
  if (m < 0) m = 0;
  if (m > 1) m = 1;
  if (m < 1e-9) return { sn: Math.sin(u), cn: Math.cos(u), dn: 1, am: u };
  if (m > 1 - 1e-9) {
    const t = Math.tanh(u);
    const s = 1 / Math.cosh(u);
    return { sn: t, cn: s, dn: s, am: 2 * Math.atan(Math.exp(u)) - Math.PI / 2 };
  }
  // Build Landen sequence
  const MAX = 24;
  const a = [1.0], b = [Math.sqrt(1 - m)], c = [Math.sqrt(m)];
  let N = 0;
  for (let n = 0; n < MAX; n++) {
    a.push((a[n] + b[n]) / 2);
    b.push(Math.sqrt(a[n] * b[n]));
    c.push((a[n] - b[n]) / 2);
    N = n + 1;
    if (Math.abs(c[N]) < 1e-15) break;
  }
  let phi = Math.pow(2, N) * a[N] * u;
  for (let n = N; n > 0; n--) {
    phi = (Math.asin(c[n] / a[n] * Math.sin(phi)) + phi) / 2;
  }
  const sn = Math.sin(phi), cn = Math.cos(phi);
  const dn = Math.sqrt(1 - m * sn * sn);
  return { sn, cn, dn, am: phi };
}

// ── Elastica curvature κ(s) ─────────────────────────────────────────────
// Inflectional (k ∈ (0,1)): κ(s) = 2k · sn(s | k²)   [changes sign]
// Separatrix   (k = 1):     κ(s) = 2 · sech(s)         [Euler spiral]
// Non-infl.    (m ∈ (0,1)): κ(s) = 2 · dn(s | m)       [stays positive]
function kappaInflectional(s, k) {
  return 2 * k * ellipj(s, k * k).sn;
}
function kappaSeparatrix(s) {
  return 2 / Math.cosh(s);
}
function kappaNonInflectional(s, m) {
  return 2 * ellipj(s, m).dn;
}

// ── Integrate the elastica ODE ──────────────────────────────────────────
// Given κ(s), integrate θ(s) = ∫κds and (x,y)(s) = ∫(cosθ, sinθ)ds.
// Uses midpoint rule with N steps over [sMin, sMax].
function integrateElastica(kappaFn, sMin, sMax, N, theta0) {
  const ds = (sMax - sMin) / N;
  let x = 0, y = 0, theta = theta0 !== undefined ? theta0 : 0;
  const pts = [{ x, y, theta, kappa: kappaFn(sMin), s: sMin }];

  for (let i = 0; i < N; i++) {
    const s = sMin + i * ds;
    // Midpoint rule for better accuracy
    const kMid = kappaFn(s + ds / 2);
    const thetaMid = theta + kMid * (ds / 2);
    x += Math.cos(thetaMid) * ds;
    y += Math.sin(thetaMid) * ds;
    theta += kMid * ds;
    pts.push({ x, y, theta, kappa: kappaFn(s + ds), s: s + ds });
  }
  return pts;
}

// ── Figure width helper ─────────────────────────────────────────────────
function figW(el) {
  return el.getBoundingClientRect().width
    || (el.parentElement && el.parentElement.getBoundingClientRect().width)
    || 680;
}

// ══════════════════════════════════════════════════════════════════════
// §1  Elastica explorer (interactive dual-panel: κ(s) | curve)
// ══════════════════════════════════════════════════════════════════════
function drawElastica() {
  const container = document.getElementById('elastica-container');
  const svg = document.getElementById('fig-elastica');
  if (!container || !svg) return;

  // Read slider
  const kInput = document.getElementById('elastica-k');
  const typeInput = document.getElementById('elastica-type');
  if (!kInput || !typeInput) return;

  const k = parseFloat(kInput.value);
  const type = typeInput.value;  // "inflectional" | "noninflectional"

  // Update displayed value
  const kLabel = document.getElementById('k-label');
  if (kLabel) kLabel.textContent = k.toFixed(2);

  // Compute curve
  let kappaFn, sMin, sMax, curveLabel, kLabel2;
  const N = 800;

  if (type === 'inflectional') {
    if (k < 0.001) {
      kappaFn = () => 0;
      sMin = -6; sMax = 6;
      curveLabel = 'straight line (k → 0)';
    } else {
      const m = k * k;
      const Km = ellipticK(m);
      kappaFn = s => kappaInflectional(s, k);
      sMin = -2 * Km; sMax = 2 * Km;  // one full period
      curveLabel = `inflectional  k = ${k.toFixed(2)},  period = 4K(k²) = ${(4*Km).toFixed(2)}`;
    }
  } else {
    // non-inflectional: κ = 2·dn(s, k²)
    const m = k * k;
    const Km = ellipticK(m);
    kappaFn = s => kappaNonInflectional(s, m);
    sMin = -Km; sMax = Km;  // one period
    curveLabel = `non-inflectional  m = k² = ${m.toFixed(2)},  half-period = K(m) = ${Km.toFixed(2)}`;
  }

  const pts = integrateElastica(kappaFn, sMin, sMax, N);

  // Centre the curve
  const xMid = (d3.min(pts, p => p.x) + d3.max(pts, p => p.x)) / 2;
  const yMid = (d3.min(pts, p => p.y) + d3.max(pts, p => p.y)) / 2;
  pts.forEach(p => { p.xc = p.x - xMid; p.yc = p.y - yMid; });

  const W = figW(svg);
  const H = 360;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const g = d3.select(svg); g.selectAll('*').remove();

  const half = W / 2 - 4;
  const mLeft  = { t: 24, b: 40, l: 52, r: 8 };
  const mRight = { t: 24, b: 40, l: 8,  r: 16 };

  // ── Left panel: κ(s) ──────────────────────────────────────────────
  const LW = half, LH = H;
  const sArr = pts.map(p => p.s);
  const kArr = pts.map(p => p.kappa);

  const xS = d3.scaleLinear().domain([sMin, sMax]).range([mLeft.l, LW - mLeft.r]);
  const yK = d3.scaleLinear()
    .domain([d3.min(kArr) * 1.2, d3.max(kArr) * 1.2])
    .nice()
    .range([LH - mLeft.b, mLeft.t]);

  const leftG = g.append('g');
  // Zero line
  leftG.append('line').attr('x1', mLeft.l).attr('x2', LW - mLeft.r)
    .attr('y1', yK(0)).attr('y2', yK(0))
    .attr('stroke', '#e8e8e8').attr('stroke-width', 1);
  // κ(s) curve
  leftG.append('path')
    .attr('d', d3.line().x(p => xS(p.s)).y(p => yK(p.kappa))(pts))
    .attr('fill', 'none')
    .attr('stroke', type === 'inflectional' ? '#1565c0' : '#2e7d32')
    .attr('stroke-width', 2);

  leftG.append('g').attr('transform', `translate(0,${LH - mLeft.b})`)
    .call(d3.axisBottom(xS).ticks(5))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  leftG.append('g').attr('transform', `translate(${mLeft.l},0)`)
    .call(d3.axisLeft(yK).ticks(5))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));

  // κ formula label
  const kappaStr = type === 'inflectional'
    ? `κ(s) = 2k · sn(s | k²)`
    : `κ(s) = 2 · dn(s | k²)`;
  leftG.append('text').attr('x', (LW + mLeft.l) / 2).attr('y', mLeft.t - 6)
    .attr('text-anchor', 'middle')
    .attr('font-family', 'Source Sans 3').attr('font-size', 11).attr('fill', '#444')
    .text(kappaStr);

  leftG.append('text').attr('x', (LW + mLeft.l) / 2).attr('y', LH - 10)
    .attr('text-anchor', 'middle')
    .attr('font-family', 'Source Sans 3').attr('font-size', 11).attr('fill', '#777')
    .text('arc length  s');

  leftG.append('text')
    .attr('transform', `translate(14,${(LH + mLeft.t) / 2})rotate(-90)`)
    .attr('text-anchor', 'middle')
    .attr('font-family', 'Source Sans 3').attr('font-size', 11).attr('fill', '#777')
    .text('curvature  κ(s)');

  // ── Divider ─────────────────────────────────────────────────────────
  g.append('line')
    .attr('x1', half).attr('x2', half)
    .attr('y1', mLeft.t).attr('y2', H - mLeft.b)
    .attr('stroke', '#e8e8e8').attr('stroke-width', 1);

  // ── Right panel: (x, y) curve ─────────────────────────────────────
  const RW = W - half, RH = H;
  const rOff = half;

  const xExt = [d3.min(pts, p => p.xc), d3.max(pts, p => p.xc)];
  const yExt = [d3.min(pts, p => p.yc), d3.max(pts, p => p.yc)];
  const range = Math.max(xExt[1] - xExt[0], yExt[1] - yExt[0]) / 2 + 0.5;
  const cx = (xExt[0] + xExt[1]) / 2, cy = (yExt[0] + yExt[1]) / 2;

  const xC = d3.scaleLinear().domain([cx - range, cx + range]).range([rOff + mRight.l, W - mRight.r]);
  const yC = d3.scaleLinear().domain([cy - range, cy + range]).range([RH - mRight.b, mRight.t]);

  const rightG = g.append('g');

  // Axes
  if (yC(0) > mRight.t && yC(0) < RH - mRight.b) {
    rightG.append('line').attr('x1', rOff + mRight.l).attr('x2', W - mRight.r)
      .attr('y1', yC(0)).attr('y2', yC(0)).attr('stroke', '#f0f0f0').attr('stroke-width', 1);
  }
  if (xC(0) > rOff + mRight.l && xC(0) < W - mRight.r) {
    rightG.append('line').attr('x1', xC(0)).attr('x2', xC(0))
      .attr('y1', mRight.t).attr('y2', RH - mRight.b).attr('stroke', '#f0f0f0').attr('stroke-width', 1);
  }

  // Colour gradient by arc length
  const defs = g.append('defs');
  const gradId = 'elast-grad';
  const grad = defs.append('linearGradient').attr('id', gradId)
    .attr('gradientUnits', 'userSpaceOnUse')
    .attr('x1', xC(pts[0].xc)).attr('y1', yC(pts[0].yc))
    .attr('x2', xC(pts[pts.length-1].xc)).attr('y2', yC(pts[pts.length-1].yc));
  const baseColor = type === 'inflectional' ? '#1565c0' : '#2e7d32';
  grad.append('stop').attr('offset', '0%').attr('stop-color', '#90caf9');
  grad.append('stop').attr('offset', '50%').attr('stop-color', baseColor);
  grad.append('stop').attr('offset', '100%').attr('stop-color', '#90caf9');

  rightG.append('path')
    .attr('d', d3.line().x(p => xC(p.xc)).y(p => yC(p.yc))(pts))
    .attr('fill', 'none')
    .attr('stroke', `url(#${gradId})`)
    .attr('stroke-width', 2.2);

  // Start/end points
  rightG.append('circle').attr('cx', xC(pts[0].xc)).attr('cy', yC(pts[0].yc))
    .attr('r', 3.5).attr('fill', baseColor);
  rightG.append('circle').attr('cx', xC(pts[pts.length-1].xc)).attr('cy', yC(pts[pts.length-1].yc))
    .attr('r', 3.5).attr('fill', baseColor);

  rightG.append('g').attr('transform', `translate(0,${RH - mRight.b})`)
    .call(d3.axisBottom(xC).ticks(4))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));

  // Type label top-right
  rightG.append('text').attr('x', W - mRight.r - 4).attr('y', mRight.t - 6)
    .attr('text-anchor', 'end')
    .attr('font-family', 'Source Sans 3').attr('font-size', 10.5).attr('fill', '#888')
    .text(type === 'inflectional' ? 'inflectional elastica' : 'non-inflectional elastica');
}

// ══════════════════════════════════════════════════════════════════════
// §2  Euler / Cornu spiral (separatrix, k=1)
// ══════════════════════════════════════════════════════════════════════
function drawCornuSpiral() {
  const svg = document.getElementById('fig-cornu');
  if (!svg) return;

  const W = figW(svg), H = 320;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const g = d3.select(svg); g.selectAll('*').remove();

  // κ(s) = 2·sech(s), integrated numerically
  const N = 2000, S = 12;
  const pts = integrateElastica(kappaSeparatrix, -S, S, N);

  // Centre
  const xMid = (d3.min(pts, p=>p.x) + d3.max(pts, p=>p.x)) / 2;
  const yMid = (d3.min(pts, p=>p.y) + d3.max(pts, p=>p.y)) / 2;
  pts.forEach(p => { p.xc = p.x - xMid; p.yc = p.y - yMid; });

  const m = { t: 20, b: 36, l: 16, r: 16 };
  const range = Math.max(
    d3.max(pts, p=>Math.abs(p.xc)),
    d3.max(pts, p=>Math.abs(p.yc))
  ) + 0.3;

  const xC = d3.scaleLinear().domain([-range, range]).range([m.l, W - m.r]);
  const yC = d3.scaleLinear().domain([-range, range]).range([H - m.b, m.t]);

  // Axes
  g.append('line').attr('x1', m.l).attr('x2', W - m.r)
    .attr('y1', yC(0)).attr('y2', yC(0)).attr('stroke', '#eee').attr('stroke-width', 1);
  g.append('line').attr('x1', xC(0)).attr('x2', xC(0))
    .attr('y1', m.t).attr('y2', H - m.b).attr('stroke', '#eee').attr('stroke-width', 1);

  // Colour by arc length: purple → teal → gold
  const colorScale = d3.scaleSequential(d3.interpolatePlasma).domain([-S, S]);

  // Draw as coloured segments
  for (let i = 0; i < pts.length - 1; i += 2) {
    const p0 = pts[i], p1 = pts[Math.min(i+2, pts.length-1)];
    g.append('line')
      .attr('x1', xC(p0.xc)).attr('y1', yC(p0.yc))
      .attr('x2', xC(p1.xc)).attr('y2', yC(p1.yc))
      .attr('stroke', colorScale(p0.s)).attr('stroke-width', 2.2)
      .attr('stroke-linecap', 'round');
  }

  // Origin dot (s=0, inflection point, κ=2 max)
  const origin = pts.find(p => Math.abs(p.s) < 0.01) || pts[N];
  g.append('circle').attr('cx', xC(origin.xc)).attr('cy', yC(origin.yc))
    .attr('r', 5).attr('fill', '#b71c1c').attr('fill-opacity', 0.85);
  g.append('text').attr('x', xC(origin.xc) + 8).attr('y', yC(origin.yc) + 4)
    .attr('font-family', 'Source Sans 3').attr('font-size', 11).attr('fill', '#b71c1c')
    .text('s = 0,  κ_max = 2');

  // Asymptote label
  const pEnd = pts[pts.length - 1];
  g.append('text').attr('x', xC(pEnd.xc) - 4).attr('y', yC(pEnd.yc) - 6)
    .attr('text-anchor', 'end')
    .attr('font-family', 'Source Sans 3').attr('font-size', 10).attr('fill', '#888')
    .text('s → +∞');

  g.append('text').attr('x', W - m.r).attr('y', H - m.b + 14)
    .attr('text-anchor', 'end')
    .attr('font-family', 'Source Sans 3').attr('font-size', 10).attr('fill', '#888')
    .text('x');
  g.append('text').attr('x', m.l + 2).attr('y', m.t + 10)
    .attr('font-family', 'Source Sans 3').attr('font-size', 10).attr('fill', '#888')
    .text('y');
}

// ══════════════════════════════════════════════════════════════════════
// §3  K(k²) gives the spatial period of curvature oscillation
// ══════════════════════════════════════════════════════════════════════
function drawPeriodPlot() {
  const svg = document.getElementById('fig-period');
  if (!svg) return;

  const W = figW(svg), H = 270;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const g = d3.select(svg); g.selectAll('*').remove();
  const m = { t: 18, b: 40, l: 56, r: 24 };

  const kPts = d3.range(0.001, 0.999, 0.002);
  // Inflectional period = 4K(k²)
  const datInfl = kPts.map(k => ({ k, T: 4 * ellipticK(k * k) }));
  // Non-inflectional half-period = 2K(k²) [dn has period 2K]
  const datNon = kPts.map(k => ({ k, T: 2 * ellipticK(k * k) }));

  const xK = d3.scaleLinear().domain([0, 1]).range([m.l, W - m.r]);
  const yT = d3.scaleLog().domain([Math.PI, 100]).range([H - m.b, m.t]);

  const plot = g;

  // Chaotic-layer shading
  plot.append('rect')
    .attr('x', xK(0.9)).attr('width', xK(1) - xK(0.9))
    .attr('y', m.t).attr('height', H - m.t - m.b)
    .attr('fill', '#b71c1c').attr('fill-opacity', 0.06);

  // Reference: circular period = π
  plot.append('line').attr('x1', m.l).attr('x2', W - m.r)
    .attr('y1', yT(Math.PI)).attr('y2', yT(Math.PI))
    .attr('stroke', '#e0e0e0').attr('stroke-width', 1).attr('stroke-dasharray', '4,3');

  // Inflectional
  plot.append('path')
    .attr('d', d3.line().defined(d => isFinite(yT(d.T)) && yT(d.T) > m.t)
      .x(d => xK(d.k)).y(d => yT(d.T))(datInfl))
    .attr('fill', 'none').attr('stroke', '#1565c0').attr('stroke-width', 2);

  // Non-inflectional
  plot.append('path')
    .attr('d', d3.line().defined(d => isFinite(yT(d.T)) && yT(d.T) > m.t)
      .x(d => xK(d.k)).y(d => yT(d.T))(datNon))
    .attr('fill', 'none').attr('stroke', '#2e7d32').attr('stroke-width', 2)
    .attr('stroke-dasharray', '6,3');

  // Labels
  plot.append('text').attr('x', xK(0.3)).attr('y', yT(datInfl[150].T) - 7)
    .attr('font-family', 'Source Sans 3').attr('font-size', 11).attr('fill', '#1565c0')
    .text('inflectional  T = 4K(k²)');
  plot.append('text').attr('x', xK(0.3)).attr('y', yT(datNon[150].T) + 14)
    .attr('font-family', 'Source Sans 3').attr('font-size', 11).attr('fill', '#2e7d32')
    .text('non-inflectional  T = 2K(k²)');
  plot.append('text').attr('x', xK(0.05)).attr('y', yT(Math.PI) - 5)
    .attr('font-family', 'Source Sans 3').attr('font-size', 10).attr('fill', '#aaa')
    .text('π (k → 0 limit)');
  plot.append('text').attr('x', xK(0.95)).attr('y', m.t + 14)
    .attr('text-anchor', 'middle')
    .attr('font-family', 'Source Sans 3').attr('font-size', 10).attr('fill', '#b71c1c')
    .text('T → ∞');

  // Euler spiral annotation at k=1
  plot.append('line').attr('x1', xK(0.998)).attr('x2', xK(0.998))
    .attr('y1', m.t + 5).attr('y2', H - m.b)
    .attr('stroke', '#b71c1c').attr('stroke-width', 1.5).attr('stroke-dasharray', '4,2');

  // Axes
  g.append('g').attr('transform', `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xK).ticks(5).tickFormat(d => d === 1 ? '1 (sep.)' : d.toFixed(1)))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  g.append('g').attr('transform', `translate(${m.l},0)`)
    .call(d3.axisLeft(yT).tickValues([Math.PI, 2*Math.PI, 4*Math.PI, 8*Math.PI, 16*Math.PI, 32*Math.PI, 64*Math.PI])
      .tickFormat(d => (d / Math.PI).toFixed(0) + 'π'))
    .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  g.append('text').attr('x', (W + m.l) / 2).attr('y', H - 10)
    .attr('text-anchor', 'middle')
    .attr('font-family', 'Source Sans 3').attr('font-size', 12).attr('fill', '#555')
    .text('modulus  k');
  g.append('text').attr('transform', `translate(16,${(H + m.t) / 2})rotate(-90)`)
    .attr('text-anchor', 'middle')
    .attr('font-family', 'Source Sans 3').attr('font-size', 12).attr('fill', '#555')
    .text('spatial period  T(k)');
}

// ══════════════════════════════════════════════════════════════════════
// §4  Geodesic family from a shared base point
//     Shows the exponential map of SE(2): multiple elastica from origin
// ══════════════════════════════════════════════════════════════════════
function drawGeodesicFamily() {
  const svg = document.getElementById('fig-family');
  if (!svg) return;

  const W = figW(svg), H = 380;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const g = d3.select(svg); g.selectAll('*').remove();
  const m = { t: 20, b: 36, l: 20, r: 20 };

  // Draw 10 inflectional curves (k = 0.1 … 0.95) + Euler spiral
  const kVals = [0.1, 0.2, 0.35, 0.5, 0.62, 0.74, 0.84, 0.92, 0.97, 1.0];
  const colors = kVals.map((k, i) =>
    k === 1 ? '#b71c1c' : d3.interpolateViridis(i / (kVals.length - 1)));

  const allPts = kVals.map((k, i) => {
    if (k === 1) {
      return { k, pts: integrateElastica(kappaSeparatrix, -6, 6, 600) };
    }
    const Km = ellipticK(k * k);
    return { k, pts: integrateElastica(s => kappaInflectional(s, k), -2 * Km, 2 * Km, 600) };
  });

  // Find global extent
  let xAll = [], yAll = [];
  allPts.forEach(({ pts }) => pts.forEach(p => { xAll.push(p.x); yAll.push(p.y); }));
  const xExt = d3.extent(xAll), yExt = d3.extent(yAll);
  const cx = (xExt[0] + xExt[1]) / 2, cy = (yExt[0] + yExt[1]) / 2;
  const range = Math.max(xExt[1] - xExt[0], yExt[1] - yExt[0]) / 2 + 0.5;

  const xC = d3.scaleLinear().domain([cx - range, cx + range]).range([m.l, W - m.r]);
  const yC = d3.scaleLinear().domain([cy - range, cy + range]).range([H - m.b, m.t]);

  // Axes
  if (yC(0) > m.t && yC(0) < H - m.b) {
    g.append('line').attr('x1', m.l).attr('x2', W - m.r)
      .attr('y1', yC(0)).attr('y2', yC(0)).attr('stroke', '#f0f0f0').attr('stroke-width', 1);
  }

  // Draw curves
  allPts.forEach(({ k, pts }, i) => {
    const xOff = pts[0].x, yOff = pts[0].y;
    g.append('path')
      .attr('d', d3.line().x(p => xC(p.x - xOff)).y(p => yC(p.y - yOff))(pts))
      .attr('fill', 'none')
      .attr('stroke', colors[i])
      .attr('stroke-width', k === 1 ? 2.4 : 1.8)
      .attr('stroke-opacity', 0.88);
  });

  // Origin dot
  const o = allPts[0].pts[0];
  g.append('circle').attr('cx', xC(0)).attr('cy', yC(0)).attr('r', 4).attr('fill', '#555');

  // Legend
  const legX = W - m.r - 130, legY = H - m.b - 10;
  g.append('text').attr('x', legX).attr('y', legY)
    .attr('font-family', 'Source Sans 3').attr('font-size', 10.5).attr('fill', '#888')
    .text('k = 0.1 (blue) → 1.0 (red)');
}

// ── Bootstrap ──────────────────────────────────────────────────────────
window.addEventListener('load', () => {
  // Wire up slider
  const kInput = document.getElementById('elastica-k');
  const typeInput = document.getElementById('elastica-type');
  if (kInput) kInput.addEventListener('input', drawElastica);
  if (typeInput) typeInput.addEventListener('change', drawElastica);

  // Draw all figures
  drawElastica();
  drawCornuSpiral();
  drawPeriodPlot();
  drawGeodesicFamily();

  window.addEventListener('resize', () => {
    drawElastica();
    drawCornuSpiral();
    drawPeriodPlot();
    drawGeodesicFamily();
  });
});
