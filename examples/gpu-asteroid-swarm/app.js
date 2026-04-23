// GPU Asteroid Swarm — main interactive module
// ──────────────────────────────────────────────────────────────────────────
// Depends on: global d3 (classic script), Three.js module, KaTeX auto-render
// Data files (graceful degradation if missing):
//   data/asteroids_subset.json
//   data/benchmark_results.json
//   data/validation_sample.json

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

// ───────────── small utilities ──────────────────────────────────────────
const DEG = Math.PI/180;
const GM_SUN = 0.0002959122082855911; // AU^3 / day^2

const CLASS_COLOR = {
  MBA: '#1565c0', IMB: '#0d47a1', OMB: '#1976d2',
  HIL: '#6a1b9a', TJN: '#00838f', HUN: '#8e24aa',
  MCA: '#ef6c00', NEA: '#c62828', JFC: '#5d4037', UNK: '#757575',
  AMO: '#d84315', APO: '#b71c1c', ATE: '#ad1457', IEO: '#880e4f',
  AST: '#455a64', CEN: '#2e7d32', TNO: '#1b5e20',
};
const CLASS_NAME = {
  MBA:'Main belt', IMB:'Inner main belt', OMB:'Outer main belt',
  HIL:'Hilda', TJN:'Trojan', HUN:'Hungaria',
  MCA:'Mars crosser', NEA:'NEA', JFC:'Jupiter-family',
  UNK:'unknown', AMO:'Amor', APO:'Apollo', ATE:'Aten',
  IEO:'Atira/IEO', AST:'asteroid', CEN:'Centaur', TNO:'TNO',
};
const MAIN_BELT = new Set(['MBA', 'IMB', 'OMB']);
const NEA_CLASSES = new Set(['NEA', 'AMO', 'APO', 'ATE', 'IEO']);

// ───────────── Kepler ───────────────────────────────────────────────────
function solveKepler(M, e) {
  let E = M + e*Math.sin(M);
  for (let k = 0; k < 16; k++) {
    const dE = (M - E + e*Math.sin(E)) / (1 - e*Math.cos(E));
    E += dE;
    if (Math.abs(dE) < 1e-13) break;
  }
  return E;
}

function bodyXYZ(b, t_days) {
  const n = Math.sqrt(GM_SUN / (b.a*b.a*b.a));
  let M = b.M0 + n * t_days;
  M = ((M + Math.PI) % (2*Math.PI) + 2*Math.PI) % (2*Math.PI) - Math.PI;
  const E = solveKepler(M, b.e);
  const ca = Math.cos(E), sa = Math.sin(E);
  const x_op = b.a*(ca - b.e);
  const y_op = b.a*Math.sqrt(1-b.e*b.e)*sa;
  const cw=Math.cos(b.w), sw=Math.sin(b.w);
  const cO=Math.cos(b.Om), sO=Math.sin(b.Om);
  const ci=Math.cos(b.i), si=Math.sin(b.i);
  return {
    x: (cO*cw - sO*sw*ci)*x_op + (-cO*sw - sO*cw*ci)*y_op,
    y: (sO*cw + cO*sw*ci)*x_op + (-sO*sw + cO*cw*ci)*y_op,
    z: (sw*si)*x_op + (cw*si)*y_op,
    E, r: b.a*(1 - b.e*ca),
  };
}

// ───────────── elliptic integrals (Carlson; mirrors the repo's JS) ─────
function carlsonRF(x, y, z) {
  const cr = 0.0027;
  let xn=x, yn=y, zn=z;
  for (let i = 0; i < 30; i++) {
    const A = (xn + yn + zn)/3;
    const X = 1 - xn/A, Y = 1 - yn/A, Z = -(X+Y);
    if (Math.max(Math.abs(X), Math.abs(Y), Math.abs(Z)) < cr) break;
    const sx=Math.sqrt(xn), sy=Math.sqrt(yn), sz=Math.sqrt(zn);
    const lam = sx*sy + sy*sz + sz*sx;
    xn = (xn+lam)/4; yn = (yn+lam)/4; zn = (zn+lam)/4;
  }
  const A = (xn+yn+zn)/3;
  const X = 1-xn/A, Y = 1-yn/A, Z = -(X+Y);
  const E2 = X*Y - Z*Z, E3 = X*Y*Z;
  return (1/Math.sqrt(A)) * (1 - E2/10 + E3/14 + E2*E2/24 - 3*E2*E3/44);
}
function carlsonRD(x, y, z) {
  const cr = 0.0015;
  let xn=x, yn=y, zn=z, S=0, pow4=1;
  for (let i = 0; i < 30; i++) {
    const A = (xn+yn+3*zn)/5;
    const X = 1-xn/A, Y = 1-yn/A, Z = 1-zn/A;
    if (Math.max(Math.abs(X),Math.abs(Y),Math.abs(Z)) < cr) break;
    const sx=Math.sqrt(xn), sy=Math.sqrt(yn), sz=Math.sqrt(zn);
    const lam = sx*sy + sy*sz + sz*sx;
    S += pow4/(sz*(zn+lam));  pow4 /= 4;
    xn = (xn+lam)/4; yn=(yn+lam)/4; zn=(zn+lam)/4;
  }
  const A = (xn+yn+3*zn)/5;
  const X = 1-xn/A, Y = 1-yn/A, Z = 1-zn/A;
  const E2 = X*Y - 6*Z*Z;
  const E3 = 3*X*Y*Z - 8*Z*Z*Z;
  const E4 = 3*(X*Y - Z*Z)*Z*Z;
  const E5 = X*Y*Z*Z*Z;
  const poly = 1 - 3*E2/14 + E3/6 + 9*E2*E2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26;
  return 3*S + pow4*poly/(A*Math.sqrt(A));
}
function ellipticE_complete(m) {
  // Direct Carlson form for E(π/2 | m), independent of ellipticE to avoid recursion.
  // E(m) = RF(0, 1-m, 1) - (m/3) RD(0, 1-m, 1)
  if (m >= 1) return 1;
  return carlsonRF(0, 1 - m, 1) - (m/3)*carlsonRD(0, 1 - m, 1);
}
function ellipticE(phi, m) {
  // Incomplete E(φ | m) extended to all real φ via the identity
  //   E(φ + kπ | m) = 2k E(m) + E(φ | m),   |φ| ≤ π/2.
  // The Carlson-form below is only valid on [-π/2, π/2]; reducing first fixes
  // the spurious negative / oscillating values we'd otherwise see for E > π/2.
  if (phi === 0) return 0;
  const sign = phi < 0 ? -1 : 1;
  const absPhi = Math.abs(phi);
  const k = Math.round(absPhi / Math.PI);
  const phi0 = absPhi - k*Math.PI;     // phi0 ∈ (-π/2, π/2]
  const s = Math.sin(phi0), c = Math.cos(phi0), d2 = 1 - m*s*s;
  const E0 = s*carlsonRF(c*c, d2, 1) - (m*s*s*s/3)*carlsonRD(c*c, d2, 1);
  return sign * (E0 + 2*k*ellipticE_complete(m));
}

// ───────────── data loading ─────────────────────────────────────────────
async function loadJSON(path) {
  try {
    const r = await fetch(path);
    if (!r.ok) throw new Error(r.statusText);
    return await r.json();
  } catch (err) {
    console.warn(`Could not load ${path}: ${err}`);
    return null;
  }
}

// synthetic fallback if the subset file can't be loaded
function syntheticSubset(n = 3000) {
  const bodies = [];
  const groups = [
    ['MBA',0.45, 2.15,2.5, 0.3, 7], ['MBA',0.30, 2.5,2.82, 0.25, 8],
    ['MBA',0.25, 2.82,3.27, 0.3, 10], ['HIL',0.05, 3.9,4.0, 0.2, 8],
    ['TJN',0.10, 5.05,5.35, 0.12, 15], ['HUN',0.04, 1.78,2.0, 0.2, 20],
    ['MCA',0.05, 1.7,2.6, 0.5, 15], ['NEA',0.06, 0.6,1.6, 0.8, 18],
  ];
  for (const [cls, frac, a0, a1, emax, isig] of groups) {
    const m = Math.floor(n*frac);
    for (let k = 0; k < m; k++) {
      const a = a0 + Math.random()*(a1-a0);
      const e = Math.min(emax, Math.pow(Math.random(), 4) * emax * 1.4);
      const i = Math.abs(randn())*isig*DEG;
      bodies.push({
        a, e, i,
        Om: Math.random()*2*Math.PI,
        w: Math.random()*2*Math.PI,
        M0: Math.random()*2*Math.PI,
        c: cls,
      });
    }
  }
  return {count: bodies.length, epoch_jd: 2460000.5, bodies};
}
function randn() {
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);
}

// ───────────── responsive SVG width helper ──────────────────────────────
function figW(el) { return el.parentElement.clientWidth - 48; }

// clipPath helper — returns a <g> that clips to the plot rectangle
let _clipSeq = 0;
function clipGroup(g, m, W, H) {
  const id = `clip-${++_clipSeq}`;
  g.append('defs').append('clipPath').attr('id', id)
    .append('rect')
      .attr('x', m.l).attr('y', m.t)
      .attr('width', Math.max(0, W - m.l - m.r))
      .attr('height', Math.max(0, H - m.t - m.b));
  return g.append('g').attr('clip-path', `url(#${id})`);
}

// ═══════════════════════════════════════════════════════════════════════
// §1 Dataset — population scatter + histograms
// ═══════════════════════════════════════════════════════════════════════
function drawPopulation(bodies) {
  const svg = document.getElementById('fig-pop');
  const W = figW(svg), H = 380;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  // extra top margin to carry the planet axis above the plot
  const m = {t:48, r:16, b:44, l:54};

  const g = d3.select(svg); g.selectAll('*').remove();

  const colorMode = document.getElementById('pop-color').value;
  const x = d3.scaleLinear().domain([0, 6]).range([m.l, W-m.r]);
  const y = d3.scaleLinear().domain([0, 1]).range([H-m.b, m.t]);

  g.append('g').attr('transform', `translate(0,${H-m.b})`)
   .call(d3.axisBottom(x).ticks(6).tickFormat(d => d+'AU'))
   .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  g.append('g').attr('transform', `translate(${m.l},0)`)
   .call(d3.axisLeft(y).ticks(5))
   .call(s => s.selectAll('text').attr('font-family', 'Source Sans 3').attr('font-size', 11));
  g.append('text').attr('x', (W+m.l)/2).attr('y', H-8)
   .attr('font-family','Source Sans 3').attr('font-size', 12)
   .attr('text-anchor','middle').attr('fill','#555')
   .text('semi-major axis a [AU]');
  g.append('text').attr('transform', `translate(14,${(H+m.t)/2})rotate(-90)`)
   .attr('font-family','Source Sans 3').attr('font-size', 12)
   .attr('text-anchor','middle').attr('fill','#555')
   .text('eccentricity e');

  // ── planet markers on a second axis above the plot ──
  const PLANETS_AU = [
    {n:'Mercury', a:0.387, c:'#9e9e9e'},
    {n:'Venus',   a:0.723, c:'#d4a35a'},
    {n:'Earth',   a:1.000, c:'#3d78b8'},
    {n:'Mars',    a:1.524, c:'#c0602a'},
    {n:'Jupiter', a:5.204, c:'#b98b5a'},
  ];
  const pg = g.append('g').attr('class','planet-axis');
  // horizontal baseline above the plot
  pg.append('line')
    .attr('x1', m.l).attr('x2', W-m.r)
    .attr('y1', m.t-6).attr('y2', m.t-6)
    .attr('stroke', '#d0d0d0').attr('stroke-width', 0.8);
  for (const p of PLANETS_AU) {
    if (p.a < x.domain()[0] || p.a > x.domain()[1]) continue;
    const xp = x(p.a);
    // tick into the plot area
    pg.append('line')
      .attr('x1', xp).attr('x2', xp)
      .attr('y1', m.t-6).attr('y2', H-m.b)
      .attr('stroke', p.c).attr('stroke-width', 0.8)
      .attr('stroke-dasharray', '2,3').attr('opacity', 0.55);
    // dot on the planet axis
    pg.append('circle')
      .attr('cx', xp).attr('cy', m.t-6).attr('r', 3.2)
      .attr('fill', p.c).attr('stroke', '#fff').attr('stroke-width', 1);
    // label above, with AU value
    pg.append('text')
      .attr('x', xp).attr('y', m.t-14)
      .attr('text-anchor','middle')
      .attr('font-family','Source Sans 3').attr('font-size', 10.5)
      .attr('font-weight', 600).attr('fill', '#444')
      .text(p.n);
    pg.append('text')
      .attr('x', xp).attr('y', m.t-25)
      .attr('text-anchor','middle')
      .attr('font-family','JetBrains Mono').attr('font-size', 9.5)
      .attr('fill', '#888')
      .text(p.a.toFixed(2) + ' AU');
  }

  const colorFn = colorMode === 'class'
    ? (b) => CLASS_COLOR[b.c] || '#666'
    : colorMode === 'e'
      ? d3.scaleSequential(d3.interpolateTurbo).domain([0, 0.7])
      : d3.scaleSequential(d3.interpolateTurbo).domain([0, 35*DEG]);

  // all data drawn inside clip rectangle so nothing leaks past the axes
  const plot = clipGroup(g, m, W, H);

  // draw all bodies — subsample if huge for responsiveness
  const stride = Math.max(1, Math.floor(bodies.length / 4000));
  plot.append('g').selectAll('circle')
    .data(bodies.filter((_, i) => i % stride === 0))
    .join('circle')
      .attr('cx', d => x(d.a))
      .attr('cy', d => y(d.e))
      .attr('r', 1.3)
      .attr('fill', d => colorMode === 'class' ? colorFn(d) : colorFn(d[colorMode]))
      .attr('fill-opacity', 0.65);

  // Kirkwood gap markers
  const gaps = [
    [2.50, '3:1'], [2.82, '5:2'], [2.95, '7:3'], [3.27, '2:1'],
  ];
  plot.append('g').selectAll('line').data(gaps).join('line')
    .attr('x1', d => x(d[0])).attr('x2', d => x(d[0]))
    .attr('y1', m.t).attr('y2', H-m.b)
    .attr('stroke', '#888').attr('stroke-dasharray', '2,3').attr('stroke-width', 0.7);
  plot.append('g').selectAll('text').data(gaps).join('text')
    .attr('x', d => x(d[0])+3).attr('y', m.t+11)
    .attr('font-family','Source Sans 3').attr('font-size', 10)
    .attr('fill', '#888').text(d => d[1]);

  // legend
  const leg = document.getElementById('pop-legend');
  leg.innerHTML = '';
  if (colorMode === 'class') {
    const classes = Array.from(new Set(bodies.map(b => b.c))).sort();
    for (const c of classes) {
      const span = document.createElement('span');
      span.className = 'leg';
      span.innerHTML = `<span class="dot" style="background:${CLASS_COLOR[c]||'#666'}"></span>${CLASS_NAME[c]||c}`;
      leg.appendChild(span);
    }
  } else {
    const dom = colorMode === 'e' ? [0, 0.7] : [0, 35];
    const unit = colorMode === 'e' ? '' : '°';
    leg.innerHTML = `<span class="leg">0${unit}</span><span class="leg" style="width:140px;height:10px;
      background:linear-gradient(90deg,${Array.from({length:8}, (_,k)=>d3.interpolateTurbo(k/7)).join(',')})"></span><span class="leg">${dom[1]}${unit}</span>`;
  }
}

function drawHistograms(bodies) {
  const svg = document.getElementById('fig-hist');
  const W = figW(svg), H = 200;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const g = d3.select(svg); g.selectAll('*').remove();
  const panelW = (W - 40) / 2;

  function hist(data, panelX, label, dom, bins, color) {
    const m = {t:10, b:32, l:40, r:10};
    const H2 = H;
    const x = d3.scaleLinear().domain(dom).range([panelX+m.l, panelX+panelW-m.r]);
    const h = d3.bin().domain(dom).thresholds(bins)(data);
    const y = d3.scaleLinear().domain([0, d3.max(h, d=>d.length)]).range([H2-m.b, m.t]);
    g.append('g').attr('transform', `translate(0,${H2-m.b})`)
      .call(d3.axisBottom(x).ticks(5))
      .call(s => s.selectAll('text').attr('font-family','Source Sans 3').attr('font-size',10));
    g.append('g').attr('transform', `translate(${panelX+m.l},0)`)
      .call(d3.axisLeft(y).ticks(4))
      .call(s => s.selectAll('text').attr('font-family','Source Sans 3').attr('font-size',10));
    g.append('g').selectAll('rect').data(h).join('rect')
      .attr('x', d => x(d.x0)+0.5).attr('y', d => y(d.length))
      .attr('width', d => Math.max(0.5, x(d.x1)-x(d.x0)-1))
      .attr('height', d => (H2-m.b) - y(d.length))
      .attr('fill', color).attr('fill-opacity', 0.75);
    g.append('text').attr('x', panelX+panelW/2).attr('y', H2-8)
      .attr('font-family','Source Sans 3').attr('font-size', 11)
      .attr('text-anchor','middle').attr('fill','#555').text(label);
  }

  hist(bodies.map(b=>b.e),         0,           'eccentricity',    [0, 0.95], 40, '#1565c0');
  hist(bodies.map(b=>b.i*180/Math.PI), panelW+40, 'inclination [°]', [0, 45],   40, '#e65100');
}

// ═══════════════════════════════════════════════════════════════════════
// §2 Benchmark charts
// ═══════════════════════════════════════════════════════════════════════
const BACKEND_COLOR = { numpy_cpu:'#9ca3af', torch_cpu:'#1565c0', torch_cuda:'#e65100', jax_cpu:'#2e7d32', jax_cuda:'#00695c' };

function drawKPIs(bm) {
  const row = document.getElementById('kpi-row');
  if (!bm) return;
  const cuda = bm.results.filter(r => r.backend === 'torch_cuda');
  const numpy = bm.results.filter(r => r.backend === 'numpy_cpu');
  const topT = cuda.reduce((a,b) => a.throughput_meps > b.throughput_meps ? a : b, {throughput_meps:0});
  const topS = cuda.reduce((a,b) => (a.speedup_vs_numpy||0) > (b.speedup_vs_numpy||0) ? a : b, {speedup_vs_numpy:0});
  const maxBE = Math.max(...bm.results.map(r => r.body_epochs || r.bodies*r.epochs));
  const dev = bm.device_info?.cuda_device || '(no CUDA device reported)';

  const fmt = (n) => n >= 1e9 ? (n/1e9).toFixed(1)+'B' : n >= 1e6 ? (n/1e6).toFixed(1)+'M' : (n/1e3).toFixed(0)+'k';
  row.innerHTML = `
    <div class="kpi"><div class="label">Peak throughput</div>
      <div class="value">${topT.throughput_meps.toFixed(0)}<span class="unit">M body-epochs/s</span></div>
      <div class="sub">${topT.bodies?.toLocaleString()||''} bodies × ${topT.epochs||''} epochs</div></div>
    <div class="kpi"><div class="label">Max speedup vs NumPy</div>
      <div class="value">${(topS.speedup_vs_numpy||0).toFixed(1)}<span class="unit">×</span></div>
      <div class="sub">${topS.bodies?.toLocaleString()||''} × ${topS.epochs||''}</div></div>
    <div class="kpi"><div class="label">Largest problem tested</div>
      <div class="value">${fmt(maxBE)}<span class="unit">body-epochs</span></div>
      <div class="sub">${bm.results.length} configurations</div></div>
    <div class="kpi"><div class="label">Device</div>
      <div class="value" style="font-size:1rem;line-height:1.4">${dev}</div>
      <div class="sub">${bm.device_info?.platform || ''}</div></div>
  `;
}

function drawThroughput(bm) {
  const svg = document.getElementById('fig-throughput');
  const W = figW(svg), H = 300;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const g = d3.select(svg); g.selectAll('*').remove();
  if (!bm) return;
  const m = {t:14, r:16, b:40, l:64};

  const select = document.getElementById('bm-epochs');
  const epochs = Array.from(new Set(bm.results.map(r=>r.epochs))).sort((a,b)=>a-b);
  if (select.options.length === 0) {
    epochs.forEach(e => { const o = document.createElement('option'); o.value = o.textContent = e; select.appendChild(o); });
    select.value = epochs[epochs.length-1];
    select.onchange = () => drawThroughput(bm);
  }
  const T = +select.value;

  const rows = bm.results.filter(r => r.epochs === T);
  const x = d3.scaleLog().domain([d3.min(rows,r=>r.bodies)/1.3, d3.max(rows,r=>r.bodies)*1.3]).range([m.l, W-m.r]);
  const y = d3.scaleLog().domain([0.5, d3.max(rows,r=>r.throughput_meps)*1.3]).range([H-m.b, m.t]);

  g.append('g').attr('transform', `translate(0,${H-m.b})`)
    .call(d3.axisBottom(x).ticks(5, d=>d>=1e6?d/1e6+'M':d>=1e3?d/1e3+'k':d))
    .call(s=>s.selectAll('text').attr('font-family','Source Sans 3').attr('font-size',11));
  g.append('g').attr('transform', `translate(${m.l},0)`)
    .call(d3.axisLeft(y).ticks(5))
    .call(s=>s.selectAll('text').attr('font-family','Source Sans 3').attr('font-size',11));
  g.append('text').attr('x',(W+m.l)/2).attr('y',H-8)
    .attr('font-family','Source Sans 3').attr('font-size',12)
    .attr('text-anchor','middle').attr('fill','#555').text('N bodies');
  g.append('text').attr('transform', `translate(16,${(H+m.t)/2})rotate(-90)`)
    .attr('font-family','Source Sans 3').attr('font-size',12)
    .attr('text-anchor','middle').attr('fill','#555').text('throughput [M body-epochs / s]');

  const plot = clipGroup(g, m, W, H);
  const byB = d3.group(rows, r => r.backend);
  for (const [backend, arr] of byB) {
    const sorted = arr.slice().sort((a,b)=>a.bodies-b.bodies);
    const color = BACKEND_COLOR[backend] || '#666';
    plot.append('path').attr('fill','none').attr('stroke',color).attr('stroke-width',2)
      .attr('d', d3.line().x(r=>x(r.bodies)).y(r=>y(r.throughput_meps))(sorted));
    plot.append('g').selectAll('circle').data(sorted).join('circle')
      .attr('cx', r=>x(r.bodies)).attr('cy', r=>y(r.throughput_meps))
      .attr('r', 3.5).attr('fill', color);
  }
}

function drawSpeedup(bm) {
  const svg = document.getElementById('fig-speedup');
  const W = figW(svg), H = 260;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const g = d3.select(svg); g.selectAll('*').remove();
  if (!bm) return;
  const m = {t:14, r:16, b:40, l:64};

  const rows = bm.results.filter(r => r.backend === 'torch_cuda' && r.speedup_vs_numpy);
  if (!rows.length) return;
  const x = d3.scaleLog().domain([
    d3.min(rows, r=>r.bodies*r.epochs)/1.3,
    d3.max(rows, r=>r.bodies*r.epochs)*1.3]).range([m.l, W-m.r]);
  const y = d3.scaleLinear().domain([0, d3.max(rows,r=>r.speedup_vs_numpy)*1.15]).range([H-m.b, m.t]);

  for (const k of [10, 30]) {
    g.append('line').attr('x1',m.l).attr('x2',W-m.r).attr('y1',y(k)).attr('y2',y(k))
      .attr('stroke','#bbb').attr('stroke-dasharray','4,4');
    g.append('text').attr('x',W-m.r-4).attr('y',y(k)-4).attr('text-anchor','end')
      .attr('font-family','Source Sans 3').attr('font-size',10).attr('fill','#999').text(k+'×');
  }
  g.append('g').attr('transform', `translate(0,${H-m.b})`)
    .call(d3.axisBottom(x).ticks(5, d=>d>=1e9?d/1e9+'B':d>=1e6?d/1e6+'M':d>=1e3?d/1e3+'k':d))
    .call(s=>s.selectAll('text').attr('font-family','Source Sans 3').attr('font-size',11));
  g.append('g').attr('transform', `translate(${m.l},0)`)
    .call(d3.axisLeft(y).ticks(6))
    .call(s=>s.selectAll('text').attr('font-family','Source Sans 3').attr('font-size',11));
  g.append('text').attr('x',(W+m.l)/2).attr('y',H-8)
    .attr('font-family','Source Sans 3').attr('font-size',12)
    .attr('text-anchor','middle').attr('fill','#555').text('body-epochs = N × T');
  g.append('text').attr('transform', `translate(16,${(H+m.t)/2})rotate(-90)`)
    .attr('font-family','Source Sans 3').attr('font-size',12)
    .attr('text-anchor','middle').attr('fill','#555').text('GPU speedup vs NumPy');

  const sorted = rows.slice().sort((a,b) => a.bodies*a.epochs - b.bodies*b.epochs);
  const plot = clipGroup(g, m, W, H);
  plot.append('path').attr('fill','none').attr('stroke','#e65100').attr('stroke-width',2)
    .attr('d', d3.line().x(r=>x(r.bodies*r.epochs)).y(r=>y(r.speedup_vs_numpy))(sorted));
  plot.append('g').selectAll('circle').data(sorted).join('circle')
    .attr('cx', r=>x(r.bodies*r.epochs)).attr('cy', r=>y(r.speedup_vs_numpy))
    .attr('r', 4).attr('fill', '#e65100');
}

function drawBreakdown(bm) {
  const svg = document.getElementById('fig-breakdown');
  const W = figW(svg), H = 220;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const g = d3.select(svg); g.selectAll('*').remove();
  if (!bm) return;
  // Use the largest N × smallest T (common across all backends)
  const nMax = d3.max(bm.results, r=>r.bodies);
  const rows = bm.results.filter(r => r.bodies === nMax).sort((a,b)=>a.epochs-b.epochs);
  if (!rows.length) return;
  const pick = rows.filter(r => r.epochs === rows[0].epochs);
  const m = {t:20, r:16, b:60, l:80};

  const y = d3.scaleBand().domain(pick.map(r=>r.backend)).range([m.t, H-m.b]).padding(0.3);
  const maxT = d3.max(pick, r=>r.t_total_ms);
  const x = d3.scaleLog().domain([Math.max(1, maxT/1000), maxT*1.2]).range([m.l, W-m.r]);

  g.append('g').attr('transform', `translate(0,${H-m.b})`)
    .call(d3.axisBottom(x).ticks(5, d=>d>=1000?d/1000+'s':d+'ms'))
    .call(s=>s.selectAll('text').attr('font-family','Source Sans 3').attr('font-size',11));
  g.append('g').attr('transform', `translate(${m.l},0)`)
    .call(d3.axisLeft(y))
    .call(s=>s.selectAll('text').attr('font-family','Source Sans 3').attr('font-size',11));
  g.append('text').attr('x',(W+m.l)/2).attr('y',H-22)
    .attr('font-family','Source Sans 3').attr('font-size',12)
    .attr('text-anchor','middle').attr('fill','#555')
    .text(`runtime breakdown at ${nMax.toLocaleString()} bodies × ${pick[0].epochs} epochs`);

  const plotB = clipGroup(g, m, W, H);
  const xMin = m.l, xMax = W - m.r;
  for (const r of pick) {
    const cx0 = Math.max(xMin, Math.min(xMax, x(Math.max(1, r.t_propagate_ms))));
    const cx1 = Math.max(xMin, Math.min(xMax, x(Math.max(1, r.t_total_ms))));
    plotB.append('rect').attr('x', xMin).attr('y', y(r.backend))
      .attr('width', Math.max(0, cx0 - xMin)).attr('height', y.bandwidth())
      .attr('fill','#1565c0').attr('fill-opacity',0.85);
    plotB.append('rect').attr('x', cx0).attr('y', y(r.backend))
      .attr('width', Math.max(0, cx1 - cx0)).attr('height', y.bandwidth())
      .attr('fill','#e65100').attr('fill-opacity',0.85);
    // label outside clipping so it's always visible, and kept inside the SVG
    const tx = Math.min(cx1 + 6, W - m.r - 50);
    g.append('text').attr('x', tx).attr('y', y(r.backend)+y.bandwidth()/2+3)
      .attr('font-family','JetBrains Mono').attr('font-size',11).attr('fill','#333')
      .text(`${r.t_total_ms.toFixed(0)} ms`);
  }
  // legend
  g.append('rect').attr('x',m.l).attr('y',4).attr('width',10).attr('height',10).attr('fill','#1565c0');
  g.append('text').attr('x',m.l+14).attr('y',13).attr('font-family','Source Sans 3').attr('font-size',11).text('Kepler propagation');
  g.append('rect').attr('x',m.l+140).attr('y',4).attr('width',10).attr('height',10).attr('fill','#e65100');
  g.append('text').attr('x',m.l+154).attr('y',13).attr('font-family','Source Sans 3').attr('font-size',11).text('elliptic-integral arc length');
}

function drawBenchmarkTable(bm) {
  const tb = document.getElementById('bm-tbody');
  tb.innerHTML = '';
  if (!bm) return;
  for (const r of bm.results) {
    const cls = (r.speedup_vs_numpy||1) >= 10 ? 'hi' : (r.speedup_vs_numpy||1) >= 2 ? 'mid' : '';
    const tr = document.createElement('tr');
    if (cls) tr.className = cls;
    tr.innerHTML = `<td>${r.backend}</td>
      <td>${r.bodies.toLocaleString()}</td><td>${r.epochs}</td>
      <td>${r.t_propagate_ms.toFixed(1)}</td>
      <td>${r.t_arc_length_ms.toFixed(1)}</td>
      <td>${r.t_total_ms.toFixed(1)}</td>
      <td>${r.throughput_meps.toFixed(1)}</td>
      <td>${r.speedup_vs_numpy ? r.speedup_vs_numpy.toFixed(2)+'×' : '—'}</td>`;
    tb.appendChild(tr);
  }
}

// ═══════════════════════════════════════════════════════════════════════
// §3 Three.js swarm
// ═══════════════════════════════════════════════════════════════════════
function initSwarm(subset) {
  const canvas = document.getElementById('swarm-canvas');
  const hud = document.getElementById('swarm-hud');
  document.getElementById('sw-count').textContent = subset.bodies.length.toLocaleString();

  const renderer = new THREE.WebGLRenderer({canvas, antialias:true, alpha:false});
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  const scene = new THREE.Scene();
  scene.background = new THREE.Color('#04040e');

  function resize() {
    const w = canvas.clientWidth, h = canvas.clientHeight;
    renderer.setSize(w, h, false);
    camera.aspect = w/h; camera.updateProjectionMatrix();
  }
  const camera = new THREE.PerspectiveCamera(35, 1, 0.1, 2000);
  camera.position.set(7, 6, 8);

  const controls = new OrbitControls(camera, canvas);
  controls.enableDamping = true; controls.dampingFactor = 0.08;

  // Sun
  const sunGeom = new THREE.SphereGeometry(0.08, 24, 16);
  const sunMat  = new THREE.MeshBasicMaterial({color:'#ffe082'});
  scene.add(new THREE.Mesh(sunGeom, sunMat));
  const sunHalo = new THREE.PointLight('#ffca28', 1.2, 0, 0.2);
  scene.add(sunHalo);
  scene.add(new THREE.AmbientLight(0xffffff, 0.35));

  // Orbital plane reference grid
  const grid = new THREE.GridHelper(14, 14, 0x223355, 0x101c33);
  grid.rotation.x = Math.PI/2; grid.material.opacity = 0.4; grid.material.transparent = true;
  scene.add(grid);

  // Planets — J2000 mean elements (JPL, ecliptic / equinox of J2000).
  // Om = Ω (longitude of ascending node), ob = ϖ (longitude of perihelion),
  // L = mean longitude at J2000 epoch.
  const PLANETS = [
    {name:'Mercury', a:0.38710, e:0.20564, i:7.005, Om:48.331, ob:77.457, L:252.251, color:'#9e9e9e'},
    {name:'Venus',   a:0.72333, e:0.00678, i:3.395, Om:76.680, ob:131.602, L:181.979, color:'#e8c97a'},
    {name:'Earth',   a:1.00000, e:0.01671, i:0.000, Om:0.0,    ob:102.938, L:100.465, color:'#4fa0e8'},
    {name:'Mars',    a:1.52371, e:0.09339, i:1.850, Om:49.560, ob:336.056, L:355.453, color:'#c0602a'},
    {name:'Jupiter', a:5.20289, e:0.04839, i:1.304, Om:100.474,ob:14.728,  L:34.396,  color:'#c4956a'},
  ];
  // Two groups so the controls can toggle planets (spheres+labels) and
  // reference orbit rings independently.
  const planetsGroup = new THREE.Group(); scene.add(planetsGroup);
  const orbitsGroup  = new THREE.Group(); scene.add(orbitsGroup);

  function makeLabelSprite(text, color) {
    const canvas = document.createElement('canvas');
    canvas.width = 256; canvas.height = 64;
    const ctx = canvas.getContext('2d');
    ctx.font = '600 40px "Source Sans 3", system-ui, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    // soft dark halo so white text stays readable against bright planets
    ctx.lineWidth = 6; ctx.strokeStyle = 'rgba(0,0,10,0.85)';
    ctx.strokeText(text, 128, 32);
    ctx.fillStyle = color;
    ctx.fillText(text, 128, 32);
    const tex = new THREE.CanvasTexture(canvas);
    tex.anisotropy = 4;
    const mat = new THREE.SpriteMaterial({map: tex, transparent: true, depthTest: false});
    const sprite = new THREE.Sprite(mat);
    sprite.scale.set(0.55, 0.14, 1);   // world units
    sprite.position.set(0, 0, 0.12);   // slight z-offset above the planet
    return sprite;
  }

  const planetMeshes = PLANETS.map(p => {
    const g = new THREE.Group();
    const mesh = new THREE.Mesh(
      new THREE.SphereGeometry(0.045, 18, 12),
      new THREE.MeshBasicMaterial({color:p.color}));
    g.add(mesh);
    const label = makeLabelSprite(p.name, '#ffffff');
    label.visible = false;          // fades in when the camera is close
    g.add(label);
    planetsGroup.add(g);
    // reference orbit ring
    const N = 180, pts = [];
    for (let k = 0; k <= N; k++) {
      const nu = (k/N)*2*Math.PI;
      const r = p.a*(1-p.e*p.e)/(1+p.e*Math.cos(nu));
      const u = nu + (p.ob - p.Om)*DEG;
      const cO = Math.cos(p.Om*DEG), sO = Math.sin(p.Om*DEG);
      const ci = Math.cos(p.i*DEG),  si = Math.sin(p.i*DEG);
      pts.push(new THREE.Vector3(
        r*(cO*Math.cos(u) - sO*Math.sin(u)*ci),
        r*(sO*Math.cos(u) + cO*Math.sin(u)*ci),
        r*si*Math.sin(u)));
    }
    const ring = new THREE.Line(
      new THREE.BufferGeometry().setFromPoints(pts),
      new THREE.LineBasicMaterial({color: p.color, transparent:true, opacity:0.45}));
    orbitsGroup.add(ring);
    return {p, g, label};
  });

  // Instanced asteroid points
  const N = subset.bodies.length;
  const posArr = new Float32Array(N*3);
  const colArr = new Float32Array(N*3);
  const geom = new THREE.BufferGeometry();
  geom.setAttribute('position', new THREE.BufferAttribute(posArr, 3));
  geom.setAttribute('color',    new THREE.BufferAttribute(colArr, 3));
  const mat = new THREE.PointsMaterial({
    size: 1.5, sizeAttenuation:false, vertexColors:true,
    transparent:true, opacity:0.85,
  });
  const points = new THREE.Points(geom, mat); scene.add(points);

  function refreshColors(mode) {
    const arr = geom.attributes.color.array;
    for (let k = 0; k < N; k++) {
      const b = subset.bodies[k];
      let col;
      if (mode === 'class') col = new THREE.Color(CLASS_COLOR[b.c] || '#999');
      else if (mode === 'e') col = new THREE.Color(d3.interpolateTurbo(Math.min(1, b.e/0.7)));
      else if (mode === 'i') col = new THREE.Color(d3.interpolateTurbo(Math.min(1, b.i/(35*DEG))));
      else col = new THREE.Color(d3.interpolateTurbo(Math.min(1, (b.a-0.5)/5)));
      arr[3*k]=col.r; arr[3*k+1]=col.g; arr[3*k+2]=col.b;
    }
    geom.attributes.color.needsUpdate = true;
  }
  refreshColors('class');

  // Controls
  const speedInput = document.getElementById('sw-speed');
  const speedOut   = document.getElementById('sw-speed-v');
  const colorInput = document.getElementById('sw-color');
  const sizeInput  = document.getElementById('sw-size');
  const sizeOut    = document.getElementById('sw-size-v');
  const planetsChk = document.getElementById('sw-planets');
  const orbitsChk  = document.getElementById('sw-orbits');

  // initial values reflect the current control states so the UI matches what is drawn
  speedOut.textContent = `${speedInput.value} d/s`;
  orbitsGroup.visible  = orbitsChk.checked;
  planetsGroup.visible = planetsChk.checked;

  speedInput.addEventListener('input', ()=>{ speedOut.textContent = `${speedInput.value} d/s`; });
  colorInput.addEventListener('change', ()=> refreshColors(colorInput.value));
  sizeInput.addEventListener('input', ()=>{
    sizeOut.textContent = (+sizeInput.value).toFixed(1);
    mat.size = +sizeInput.value;
  });
  planetsChk.addEventListener('change', ()=>{ planetsGroup.visible = planetsChk.checked; });
  orbitsChk .addEventListener('change', ()=>{ orbitsGroup.visible  = orbitsChk.checked;  });

  // Time state
  let simTime = 0;   // days
  let lastReal = performance.now();

  function step() {
    const now = performance.now();
    const dt = (now - lastReal)/1000; lastReal = now;
    simTime += dt * (+speedInput.value);

    // propagate all bodies
    for (let k = 0; k < N; k++) {
      const b = subset.bodies[k];
      const s = bodyXYZ(b, simTime);
      posArr[3*k] = s.x; posArr[3*k+1] = s.y; posArr[3*k+2] = s.z;
    }
    geom.attributes.position.needsUpdate = true;

    // propagate planets (Kepler: M = L - ϖ + n·t)
    const camDist = camera.position.length();
    // show planet name labels only when the camera is close enough that they
    // don't clutter the wide-field view
    const showLabels = camDist < 3.2;
    for (const {p, g, label} of planetMeshes) {
      const n = Math.sqrt(GM_SUN / (p.a*p.a*p.a));
      const M = ((p.L*DEG - p.ob*DEG) + n*simTime);
      const E = solveKepler(M, p.e);
      const r = p.a*(1 - p.e*Math.cos(E));
      const nu = 2*Math.atan2(Math.sqrt(1+p.e)*Math.sin(E/2), Math.sqrt(1-p.e)*Math.cos(E/2));
      const u = nu + (p.ob - p.Om)*DEG;
      const cO = Math.cos(p.Om*DEG), sO = Math.sin(p.Om*DEG);
      const ci = Math.cos(p.i*DEG),  si = Math.sin(p.i*DEG);
      g.position.set(
        r*(cO*Math.cos(u) - sO*Math.sin(u)*ci),
        r*(sO*Math.cos(u) + cO*Math.sin(u)*ci),
        r*si*Math.sin(u));
      label.visible = showLabels && planetsGroup.visible;
    }

    controls.update();
    renderer.render(scene, camera);
    hud.textContent = `${N.toLocaleString()} asteroids   t = ${(simTime/365.25).toFixed(2)} yr`;
    requestAnimationFrame(step);
  }

  resize(); window.addEventListener('resize', resize);
  step();
}

// ═══════════════════════════════════════════════════════════════════════
// §4 Arc length panel
// ═══════════════════════════════════════════════════════════════════════
function initArcPanel(subset) {
  const sel = document.getElementById('arc-asteroid');
  // Pick a handful of interesting bodies: one low-e MBA, one high-e NEA, one Hilda, etc.
  const PICKS = {};
  const categories = [
    {label:'low-e main-belt',   filter:b => MAIN_BELT.has(b.c) && b.e < 0.05},
    {label:'typical main-belt', filter:b => MAIN_BELT.has(b.c) && b.e > 0.15 && b.e < 0.25},
    {label:'Hilda (3:2)',       filter:b => b.c==='HIL'},
    {label:'Jupiter Trojan',    filter:b => b.c==='TJN'},
    {label:'Mars crosser',      filter:b => b.c==='MCA'},
    {label:'high-e NEA',        filter:b => NEA_CLASSES.has(b.c) && b.e > 0.4},
  ];
  for (const c of categories) {
    const hit = subset.bodies.find(c.filter);
    if (hit) PICKS[c.label] = hit;
  }
  sel.innerHTML = '';
  Object.keys(PICKS).forEach(lab => {
    const o = document.createElement('option'); o.value = lab; o.textContent =
      `${lab} — a=${PICKS[lab].a.toFixed(2)} AU, e=${PICKS[lab].e.toFixed(3)}`;
    sel.appendChild(o);
  });

  const slider = document.getElementById('arc-E');
  const out    = document.getElementById('arc-E-v');
  const svg    = document.getElementById('fig-arc');

  function draw() {
    const body = PICKS[sel.value];
    const E = +slider.value;
    out.textContent = E.toFixed(2);
    const W = figW(svg), H = 340;
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    const g = d3.select(svg); g.selectAll('*').remove();
    const panelW = (W - 30)/2;

    // -------- left: orbit in 2D --------
    const m1 = {t:16, b:32, l:16, r:16};
    const a = body.a, e = body.e;
    const b = a*Math.sqrt(1-e*e);
    const sx = d3.scaleLinear()
      .domain([-a*(1+e)*1.1, a*(1+e)*1.1])
      .range([m1.l, panelW - m1.r]);
    const sy = d3.scaleLinear()
      .domain([-b*1.2, b*1.2])
      .range([H - m1.b, m1.t]);
    const k = Math.min((sx.range()[1]-sx.range()[0])/(sx.domain()[1]-sx.domain()[0]),
                       (sy.range()[0]-sy.range()[1])/(sy.domain()[1]-sy.domain()[0]));
    const cxPix = (sx(0));
    const cyPix = (sy(0));

    // clip the left panel drawing area so arc/planet can't leak past the box
    const leftM = {t:m1.t, b:m1.b, l:m1.l, r:W - panelW + m1.r};
    const plotL = clipGroup(g, leftM, W, H);

    // draw full ellipse (centred at origin in orbital plane; x_op from -a(1+e) to a(1-e))
    const ORBIT = [];
    for (let k2 = 0; k2 <= 240; k2++) {
      const th = (k2/240)*2*Math.PI;
      ORBIT.push([a*(Math.cos(th)-e), b*Math.sin(th)]);
    }
    plotL.append('path')
      .attr('d', d3.line().x(p=>sx(p[0])).y(p=>sy(p[1]))(ORBIT))
      .attr('fill','none').attr('stroke','#1565c0').attr('stroke-width',1.5);

    // focus / Sun
    const focusX = -a*e;
    plotL.append('circle').attr('cx',sx(focusX)).attr('cy',sy(0)).attr('r',4).attr('fill','#ffa000');
    plotL.append('text').attr('x',sx(focusX)+7).attr('y',sy(0)-7)
      .attr('font-family','Source Sans 3').attr('font-size',11).attr('fill','#555').text('Sun');

    // Arc up to E (red)
    const ARC = [];
    const steps = 120;
    for (let k2 = 0; k2 <= steps; k2++) {
      const th = (k2/steps)*E;
      ARC.push([a*(Math.cos(th)-e), b*Math.sin(th)]);
    }
    plotL.append('path').attr('d', d3.line().x(p=>sx(p[0])).y(p=>sy(p[1]))(ARC))
      .attr('fill','none').attr('stroke','#c62828').attr('stroke-width',3);

    // planet position at E
    const px = a*(Math.cos(E)-e), py = b*Math.sin(E);
    plotL.append('circle').attr('cx', sx(px)).attr('cy', sy(py)).attr('r',5).attr('fill','#c62828');

    g.append('text').attr('x', m1.l).attr('y', H-10)
      .attr('font-family','Source Sans 3').attr('font-size',11).attr('fill','#555')
      .text(`a=${a.toFixed(3)} AU, e=${e.toFixed(3)}  ·  E=${E.toFixed(3)} rad`);

    // -------- right: L(E) curve --------
    const m2 = {t:16, b:40, l:52, r:16};
    const x0 = panelW + 30;
    const xE = d3.scaleLinear().domain([0, 2*Math.PI]).range([x0+m2.l, W-m2.r]);
    const perim = 4*a*ellipticE_complete(e*e);
    // domain must cover BOTH the exact curve (→ perim) and the circular proxy (→ 2πa),
    // otherwise the proxy line slips above the top of the plot for e > 0.
    const yMax = Math.max(perim, 2*Math.PI*a) * 1.05;
    const yL = d3.scaleLinear().domain([0, yMax]).range([H-m2.b, m2.t]);

    // exact L(E) curve — ellipticE is now extended via E(φ+kπ|m)=2k·E(m)+E(φ|m)
    // so the curve is monotone on [0, 2π] and terminates at the orbit perimeter.
    const Epts = d3.range(0, 2*Math.PI*1.0001, 2*Math.PI/120);
    const Lexact = Epts.map(x => [x, a * ellipticE(x, e*e)]);
    const Lproxy = Epts.map(x => [x, a * x]);   // circular proxy 2πa at E=2π

    g.append('g').attr('transform',`translate(0,${H-m2.b})`)
      .call(d3.axisBottom(xE).ticks(5).tickFormat(d => (d/Math.PI).toFixed(1)+'π'))
      .call(s=>s.selectAll('text').attr('font-family','Source Sans 3').attr('font-size',11));
    g.append('g').attr('transform',`translate(${x0+m2.l},0)`)
      .call(d3.axisLeft(yL).ticks(5))
      .call(s=>s.selectAll('text').attr('font-family','Source Sans 3').attr('font-size',11));
    g.append('text').attr('x', (W+x0+m2.l)/2).attr('y', H-12)
      .attr('font-family','Source Sans 3').attr('font-size',12).attr('text-anchor','middle').attr('fill','#555')
      .text('eccentric anomaly E');
    g.append('text').attr('transform',`translate(${x0+16},${(H+m2.t)/2})rotate(-90)`)
      .attr('font-family','Source Sans 3').attr('font-size',12).attr('text-anchor','middle').attr('fill','#555')
      .text('arc length L(E) [AU]');

    // clip data to the right-panel rectangle so neither curve can escape
    const rightM = {t:m2.t, b:m2.b, l:x0+m2.l, r:m2.r};
    const plotR = clipGroup(g, rightM, W, H);

    plotR.append('path').attr('d', d3.line().x(p=>xE(p[0])).y(p=>yL(p[1]))(Lproxy))
      .attr('fill','none').attr('stroke','#9e9e9e').attr('stroke-dasharray','4,3');
    plotR.append('path').attr('d', d3.line().x(p=>xE(p[0])).y(p=>yL(p[1]))(Lexact))
      .attr('fill','none').attr('stroke','#c62828').attr('stroke-width',2);

    // perimeter reference line
    plotR.append('line')
      .attr('x1', x0+m2.l).attr('x2', W-m2.r)
      .attr('y1', yL(perim)).attr('y2', yL(perim))
      .attr('stroke','#c62828').attr('stroke-width',0.8)
      .attr('stroke-dasharray','2,3').attr('opacity',0.55);

    // current E marker
    const Lcur = a*ellipticE(E, e*e);
    const Lpr  = a*E;
    plotR.append('circle').attr('cx', xE(E)).attr('cy', yL(Lcur)).attr('r',4).attr('fill','#c62828');
    plotR.append('line').attr('x1', xE(E)).attr('x2', xE(E)).attr('y1', yL(Lcur)).attr('y2', yL(Lpr))
      .attr('stroke','#9e9e9e').attr('stroke-width',1);
    const relErr = Math.abs(Lcur - Lpr) / (Lcur || 1e-12);
    g.append('text').attr('x', x0+m2.l+10).attr('y', m2.t+14)
      .attr('font-family','JetBrains Mono').attr('font-size',11).attr('fill','#333')
      .text(`exact L = ${Lcur.toFixed(4)} AU`);
    g.append('text').attr('x', x0+m2.l+10).attr('y', m2.t+28)
      .attr('font-family','JetBrains Mono').attr('font-size',11).attr('fill','#666')
      .text(`circular ≈ ${Lpr.toFixed(4)} AU`);
    g.append('text').attr('x', x0+m2.l+10).attr('y', m2.t+42)
      .attr('font-family','JetBrains Mono').attr('font-size',11).attr('fill', relErr>0.01?'#c62828':'#555')
      .text(`rel err = ${(100*relErr).toFixed(2)}%`);
    g.append('text').attr('x', x0+m2.l+10).attr('y', m2.t+58)
      .attr('font-family','JetBrains Mono').attr('font-size',11).attr('fill','#333')
      .text(`perim P = ${perim.toFixed(4)} AU`);
  }

  sel.addEventListener('change', draw);
  slider.addEventListener('input', draw);
  window.addEventListener('resize', draw);
  draw();
}

// ═══════════════════════════════════════════════════════════════════════
// §5 Validation table
// ═══════════════════════════════════════════════════════════════════════
function drawValidation(val) {
  const tb = document.getElementById('val-tbody');
  tb.innerHTML = '';
  if (!val) return;
  for (const s of val.samples) {
    let cls = 'hi';
    if (s.err_au >= 1e-4) cls = 'mid';
    if (s.err_au >= 1e-2) cls = 'low';
    const pred = s.pred_xyz_au ? Math.sqrt(s.pred_xyz_au.reduce((a,b)=>a+b*b,0)).toFixed(5)+' AU' : '—';
    const ref  = s.horizons_xyz_au ? Math.sqrt(s.horizons_xyz_au.reduce((a,b)=>a+b*b,0)).toFixed(5)+' AU' : '—';
    const tr = document.createElement('tr');
    tr.className = cls;
    tr.innerHTML = `<td>${s.name}</td>
      <td>${s.dt_years.toFixed ? s.dt_years.toFixed(2) : s.dt_years}</td>
      <td>${pred}</td><td>${ref}</td>
      <td>${s.err_au.toExponential(2)}</td>`;
    tb.appendChild(tr);
  }
}

// ═══════════════════════════════════════════════════════════════════════
// Figure fullscreen toggle — adds a small ⛶ button to every .fig-box
// ═══════════════════════════════════════════════════════════════════════
function installFigureFullscreen() {
  document.querySelectorAll('.fig-box').forEach(box => {
    if (box.querySelector('.fig-expand-btn')) return;
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'fig-expand-btn';
    btn.title = 'Expand to fullscreen';
    btn.setAttribute('aria-label', 'Expand to fullscreen');
    btn.textContent = '⛶';
    btn.addEventListener('click', () => {
      const inFs = document.fullscreenElement === box;
      if (!document.fullscreenElement) {
        const req = box.requestFullscreen?.bind(box);
        if (req) {
          req().catch(() => box.classList.add('fig-expanded'));
        } else {
          box.classList.toggle('fig-expanded');
          window.dispatchEvent(new Event('resize'));
        }
      } else if (inFs) {
        document.exitFullscreen?.();
      } else {
        box.classList.toggle('fig-expanded');
        window.dispatchEvent(new Event('resize'));
      }
    });
    box.appendChild(btn);
  });
  // on fullscreen enter/exit we fire a resize so every draw() re-fits
  document.addEventListener('fullscreenchange', () => {
    window.dispatchEvent(new Event('resize'));
  });
  // Esc key exits the CSS fallback
  window.addEventListener('keydown', (ev) => {
    if (ev.key === 'Escape') {
      document.querySelectorAll('.fig-box.fig-expanded').forEach(el => el.classList.remove('fig-expanded'));
      window.dispatchEvent(new Event('resize'));
    }
  });
}

// Reproduce-section language tabs (Python ↔ MATLAB/Octave)
function installCodeTabs() {
  document.querySelectorAll('.code-tabs').forEach(wrap => {
    const btns = wrap.querySelectorAll('.tab-btn');
    const pans = wrap.querySelectorAll('.tab-panel');
    btns.forEach(btn => btn.addEventListener('click', () => {
      const key = btn.dataset.tab;
      btns.forEach(b => b.classList.toggle('active', b.dataset.tab === key));
      pans.forEach(p => p.classList.toggle('active', p.dataset.tab === key));
    }));
  });
}

// ═══════════════════════════════════════════════════════════════════════
// Main
// ═══════════════════════════════════════════════════════════════════════
(async () => {
  const [subsetRaw, bm, val] = await Promise.all([
    loadJSON('data/asteroids_subset.json'),
    loadJSON('data/benchmark_results.json'),
    loadJSON('data/validation_sample.json'),
  ]);
  const subset = subsetRaw || syntheticSubset(2500);

  // Population + hist
  drawPopulation(subset.bodies);
  drawHistograms(subset.bodies);
  document.getElementById('pop-color').addEventListener('change',
    () => drawPopulation(subset.bodies));
  window.addEventListener('resize', () => { drawPopulation(subset.bodies); drawHistograms(subset.bodies); });

  // Benchmark
  drawKPIs(bm);
  drawThroughput(bm);
  drawSpeedup(bm);
  drawBreakdown(bm);
  drawBenchmarkTable(bm);
  window.addEventListener('resize', () => { drawThroughput(bm); drawSpeedup(bm); drawBreakdown(bm); });

  // 3D swarm
  initSwarm(subset);

  // Arc length
  initArcPanel(subset);

  // Validation
  drawValidation(val);

  // UI chrome: fullscreen buttons, tabbed code panels
  installFigureFullscreen();
  installCodeTabs();
})();
