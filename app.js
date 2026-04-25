/* elliptic — site behavior (landing page)
 *
 * Renders the function finder grid + API reference from window.FUNCTIONS
 * (see data/functions.js). Handles:
 *   - orbit animation in the hero
 *   - global MATLAB/Python language toggle (persists to localStorage)
 *   - copy-to-clipboard for code-style buttons
 *   - finder search + category tabs
 *   - modal (<dialog>) with per-function detail, language tabs, numerical
 *     example, deep-dive + API reference links
 *   - collapsible API reference accordion
 */

/* ─── Language preference ─────────────────────────────────────────── */
const LANG_KEY = 'elliptic.lang';
function getLang() {
  const v = localStorage.getItem(LANG_KEY);
  return (v === 'matlab' || v === 'python') ? v : 'python';
}
function setLang(lang) {
  if (lang !== 'matlab' && lang !== 'python') return;
  localStorage.setItem(LANG_KEY, lang);
  applyLangEverywhere(lang);
}

/* ─── Orbit animation (hero SVG) ──────────────────────────────────── */
(function orbit() {
  const dot = document.getElementById('orbit-dot');
  const arc = document.getElementById('arc-path');
  if (!dot || !arc) return;
  const cx = 120, cy = 90, rx = 88, ry = 62;
  const pt = t => ({ x: cx + rx * Math.cos(t), y: cy - ry * Math.sin(t) });
  function buildArc(endAngle) {
    const n = 64; const pts = [];
    for (let i = 0; i <= n; i++) pts.push(pt((endAngle / n) * i));
    return 'M ' + pts.map(p => `${p.x.toFixed(2)} ${p.y.toFixed(2)}`).join(' L ');
  }
  let t = 0;
  (function tick() {
    t += 0.008;
    const a = t % (2 * Math.PI);
    const p = pt(a);
    dot.setAttribute('cx', p.x.toFixed(2));
    dot.setAttribute('cy', p.y.toFixed(2));
    arc.setAttribute('d', buildArc(a));
    requestAnimationFrame(tick);
  })();
})();

/* ─── Copy-to-clipboard for .btn-code[data-copy] ──────────────────── */
function wireCopyButtons(root = document) {
  root.querySelectorAll('.btn-code[data-copy]').forEach(btn => {
    if (btn.dataset.wired) return;
    btn.dataset.wired = '1';
    btn.addEventListener('click', async () => {
      const text = btn.getAttribute('data-copy');
      try { await navigator.clipboard.writeText(text); }
      catch {
        const ta = document.createElement('textarea');
        ta.value = text; document.body.appendChild(ta); ta.select();
        try { document.execCommand('copy'); } catch {}
        ta.remove();
      }
      const hint = btn.querySelector('.copy-hint');
      const original = hint ? hint.textContent : '';
      if (hint) hint.textContent = 'copied ✓';
      btn.classList.add('copied');
      setTimeout(() => {
        if (hint) hint.textContent = original || 'copy';
        btn.classList.remove('copied');
      }, 1400);
    });
  });
}

/* ─── HTML escape ─────────────────────────────────────────────────── */
function esc(s) {
  return String(s)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}

/* ─── Function finder ─────────────────────────────────────────────── */
function renderFinder() {
  const grid = document.getElementById('fnGrid');
  if (!grid || !window.FUNCTIONS) return;
  const catLabel = id => (window.CATEGORIES.find(c => c.id === id) || {}).title || id;
  grid.innerHTML = window.FUNCTIONS.map(fn => `
    <button type="button" class="fn-card" data-fn="${esc(fn.id)}" data-cat="${esc(fn.category)}"
            data-name="${esc((fn.searchTerms || '') + ' ' + fn.name + ' ' + fn.desc).toLowerCase()}">
      <span class="fn-card-name">${esc(fn.name)}</span>
      <span class="fn-card-desc">${esc(fn.desc)}</span>
      <span class="fn-card-cat">${esc(catLabel(fn.category))}</span>
    </button>
  `).join('');

  // Wire click → open modal
  grid.querySelectorAll('.fn-card').forEach(btn => {
    btn.addEventListener('click', () => openModal(btn.dataset.fn));
  });
}

function filterFinder() {
  const search  = document.getElementById('fnSearch');
  const grid    = document.getElementById('fnGrid');
  const tabs    = document.getElementById('catTabs');
  if (!search || !grid || !tabs) return;
  const q = search.value.toLowerCase().trim();
  const activeCat = tabs.querySelector('.cat-tab.active')?.dataset.cat || 'all';
  let visible = 0;
  grid.querySelectorAll('.fn-card').forEach(card => {
    const okCat = activeCat === 'all' || card.dataset.cat === activeCat;
    const okText = !q || card.dataset.name.includes(q);
    const show = okCat && okText;
    card.hidden = !show;
    if (show) visible++;
  });
  let empty = grid.querySelector('.finder-empty');
  if (visible === 0) {
    if (!empty) {
      empty = document.createElement('div');
      empty.className = 'finder-empty';
      grid.appendChild(empty);
    }
    empty.textContent = `No functions match "${search.value}"`;
  } else if (empty) {
    empty.remove();
  }
}

function wireFinder() {
  const search = document.getElementById('fnSearch');
  const tabs   = document.getElementById('catTabs');
  if (search) search.addEventListener('input', filterFinder);
  if (tabs) tabs.querySelectorAll('.cat-tab').forEach(btn => {
    btn.addEventListener('click', () => {
      tabs.querySelectorAll('.cat-tab').forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      filterFinder();
    });
  });
}

/* ─── Modal ──────────────────────────────────────────────────────── */
function renderMath(el, latex) {
  const doRender = () => {
    try { el.innerHTML = katex.renderToString(latex, { displayMode: true, throwOnError: false }); }
    catch(e) { el.textContent = latex; }
  };
  if (typeof katex !== 'undefined') { doRender(); return; }
  // KaTeX script not yet executed — show fallback, re-render on load
  el.textContent = latex;
  window.addEventListener('katex-ready', doRender, { once: true });
}

function openModal(fnId) {
  const fn = window.FUNCTIONS.find(f => f.id === fnId);
  if (!fn) return;
  const dlg = document.getElementById('fnModal');
  if (!dlg) return;

  const lang = getLang();
  const langs = ['python', 'matlab'];
  const langLabels = { python: 'Python', matlab: 'MATLAB / Octave' };
  const activeLang = (fn[lang] && !fn[lang].disabled) ? lang : (lang === 'python' ? 'matlab' : 'python');

  const deepDive = fn.deepDive
    ? `<a class="modal-link" href="${esc(fn.deepDive.href)}">Deep dive → ${esc(fn.deepDive.label)}</a>`
    : '';
  const dlmf = fn.dlmf
    ? `<a class="modal-link" href="${esc(fn.dlmf)}" target="_blank" rel="noopener">DLMF reference ↗</a>`
    : '';

  // Definition block — prefer KaTeX-rendered LaTeX, fall back to ASCII pre
  const formulaHtml = fn.formulaTeX
    ? `<div class="modal-sig-label">DEFINITION</div>
       <div class="modal-math" data-latex="${esc(fn.formulaTeX)}"></div>`
    : fn.formula
      ? `<div class="modal-sig-label">DEFINITION</div>
         <pre class="modal-formula">${esc(fn.formula)}</pre>`
      : '';

  // Relations section
  const relationsHtml = fn.relations && fn.relations.length
    ? `<div class="modal-sig-label">RELATIONS</div>
       <ul class="modal-relations">
         ${fn.relations.map(r => `
           <li class="modal-rel-item">
             <button type="button" class="modal-rel-btn" data-rel-fn="${esc(r.fn)}">${esc(r.fn)}</button>
             <span class="modal-rel-note">${esc(r.note)}</span>
           </li>`).join('')}
       </ul>`
    : '';

  const tabs = langs.map(L => {
    const available = fn[L] && !fn[L].disabled;
    const label = langLabels[L] + (available ? '' : ' · n/a');
    return `<button type="button" class="modal-lang-tab ${L === activeLang ? 'active' : ''}"
               data-lang="${L}" ${available ? '' : 'data-disabled="1"'}>${esc(label)}</button>`;
  }).join('');

  dlg.innerHTML = `
    <button type="button" class="modal-close" aria-label="Close" data-close>×</button>
    <div class="modal-header">
      <div class="modal-cat">${esc(cat(fn.category))}</div>
      <h3 class="modal-name">${esc(fn.name)}</h3>
      <p class="modal-desc">${esc(fn.longDesc || fn.desc)}</p>
    </div>
    <div class="modal-lang-tabs">${tabs}</div>
    <div class="modal-body">
      <div class="modal-sig-label">SIGNATURE</div>
      <pre class="modal-sig" data-sig></pre>
      <div class="modal-sig-label">EXAMPLE</div>
      <div class="modal-example-wrap">
        <pre class="modal-example" data-example></pre>
        <button type="button" class="btn btn-code btn-code-sm" data-copy-example aria-label="Copy example">
          <span class="copy-hint">copy</span>
        </button>
      </div>
      <div class="modal-sig-label">OUTPUT</div>
      <pre class="modal-output">${esc(fn.output)}</pre>
      ${formulaHtml}
      ${relationsHtml}
      <div class="modal-note" data-note ${fn[activeLang]?.note ? '' : 'hidden'}>${esc(fn[activeLang]?.note || '')}</div>
    </div>
    <div class="modal-footer">
      ${deepDive}
      ${dlmf}
    </div>
  `;

  // Render KaTeX for all [data-latex] elements
  dlg.querySelectorAll('[data-latex]').forEach(el => renderMath(el, el.dataset.latex));

  // Wire relation buttons → open sibling modal
  dlg.querySelectorAll('[data-rel-fn]').forEach(btn => {
    btn.addEventListener('click', () => {
      dlg.close();
      // Small delay so the closing animation finishes
      setTimeout(() => openModal(btn.dataset.relFn), 80);
    });
  });

  // helper: populate sig/example for the active language
  function paint(L) {
    const sig = dlg.querySelector('[data-sig]');
    const ex  = dlg.querySelector('[data-example]');
    const note = dlg.querySelector('[data-note]');
    const copyBtn = dlg.querySelector('[data-copy-example]');
    const src = fn[L];
    if (!src || src.disabled) {
      sig.textContent = src?.sig || '(not available)';
      ex.textContent  = src?.example || '';
    } else {
      sig.textContent = src.sig;
      ex.textContent  = src.example;
    }
    if (copyBtn) copyBtn.setAttribute('data-copy', src?.example || '');
    if (note) {
      const n = src?.note;
      if (n) { note.textContent = n; note.hidden = false; }
      else   { note.textContent = ''; note.hidden = true; }
    }
    // highlight tab
    dlg.querySelectorAll('.modal-lang-tab').forEach(t => {
      t.classList.toggle('active', t.dataset.lang === L);
    });
  }

  paint(activeLang);
  wireCopyButtons(dlg);

  // tab clicks
  dlg.querySelectorAll('.modal-lang-tab').forEach(t => {
    t.addEventListener('click', () => {
      if (t.dataset.disabled) return;
      paint(t.dataset.lang);
    });
  });
  // close button
  dlg.querySelectorAll('[data-close]').forEach(b =>
    b.addEventListener('click', () => dlg.close()));
  // backdrop click closes
  dlg.addEventListener('click', e => {
    const r = dlg.getBoundingClientRect();
    const inDialog = e.clientX >= r.left && e.clientX <= r.right &&
                     e.clientY >= r.top && e.clientY <= r.bottom;
    if (!inDialog) dlg.close();
  }, { once: true });

  dlg.showModal();
}
function cat(id) { return (window.CATEGORIES.find(c => c.id === id) || {}).title || id; }

/* ─── API reference accordion (rendered from FUNCTIONS + CATEGORIES) ─ */
function renderApiRef() {
  const root = document.getElementById('apiReference');
  if (!root || !window.FUNCTIONS || !window.CATEGORIES) return;

  const groups = window.CATEGORIES.map(c => {
    const fns = window.FUNCTIONS.filter(f => f.category === c.id);
    if (fns.length === 0) return '';
    const sections = fns.map(f => {
      const pySig = f.python && !f.python.disabled ? f.python.sig : f.python?.sig || '—';
      const mlSig = f.matlab && !f.matlab.disabled ? f.matlab.sig : f.matlab?.sig || '—';
      const formula = f.formula ? `<pre>${esc(f.formula)}</pre>` : '';
      return `
        <div class="fn-section" id="ref-${esc(f.id)}">
          <h3>${esc(f.name)}</h3>
          <div class="fn-sigs">
            <div>
              <div class="fn-sig-label">MATLAB / Octave</div>
              <div class="fn-sig">${esc(mlSig)}</div>
            </div>
            <div>
              <div class="fn-sig-label">Python</div>
              <div class="fn-sig">${esc(pySig)}</div>
            </div>
          </div>
          <p class="fn-summary">${esc(f.longDesc || f.desc)}</p>
          ${formula}
          <div class="fn-row">
            <button type="button" class="btn btn-ghost btn-sm" data-open-modal="${esc(f.id)}">Open example →</button>
            ${f.deepDive ? `<a class="btn btn-ghost btn-sm" href="${esc(f.deepDive.href)}">Deep dive: ${esc(f.deepDive.label)}</a>` : ''}
          </div>
        </div>
      `;
    }).join('');
    const note = c.note ? `<p class="api-group-note">${esc(c.note)}</p>` : '';
    return `
      <div class="api-group ${c.open ? 'open' : ''}" id="g-${esc(c.id)}">
        <div class="api-group-header" data-toggle>
          <span class="api-group-title">${esc(c.title)} <span class="api-group-count">${fns.length} function${fns.length>1?'s':''}</span></span>
          <svg class="api-group-chevron" viewBox="0 0 16 16" fill="currentColor"><path d="M8 10.94L2.06 5 3 4.06 8 9.06l5-5 .94.94z"/></svg>
        </div>
        <div class="api-group-body">
          ${note}
          ${sections}
        </div>
      </div>
    `;
  }).join('');

  root.innerHTML = groups;

  // Wire accordion toggle
  root.querySelectorAll('[data-toggle]').forEach(h => {
    h.addEventListener('click', () => h.parentElement.classList.toggle('open'));
  });
  // Wire "Open example →" → modal
  root.querySelectorAll('[data-open-modal]').forEach(b => {
    b.addEventListener('click', () => openModal(b.dataset.openModal));
  });
}

/* ─── Language toggle — update hero install chip + future bindings ── */
function applyLangEverywhere(lang) {
  // Nav toggle
  document.querySelectorAll('[data-lang-toggle]').forEach(btn => {
    btn.classList.toggle('active', btn.dataset.langToggle === lang);
  });
  // Hero install chip swaps between MATLAB git-clone and pip install
  const chip = document.getElementById('heroInstallChip');
  if (chip) {
    const span = chip.querySelector('.btn-code-text');
    if (lang === 'matlab') {
      chip.setAttribute('data-copy', 'git clone https://github.com/moiseevigor/elliptic.git && cd elliptic/matlab && matlab -r setup');
      if (span) span.textContent = 'git clone ... && cd elliptic/matlab && setup';
    } else {
      chip.setAttribute('data-copy', 'pip install elliptic-functions');
      if (span) span.textContent = 'pip install elliptic-functions';
    }
  }
  // Modal (if open): just re-paint by clicking active tab
  const dlg = document.getElementById('fnModal');
  if (dlg?.open) {
    const tab = dlg.querySelector(`.modal-lang-tab[data-lang="${lang}"]`);
    if (tab && !tab.dataset.disabled) tab.click();
  }
}

function wireLangToggle() {
  document.querySelectorAll('[data-lang-toggle]').forEach(btn => {
    btn.addEventListener('click', () => setLang(btn.dataset.langToggle));
  });
  applyLangEverywhere(getLang());
}

/* ─── Init ────────────────────────────────────────────────────────── */
document.addEventListener('DOMContentLoaded', () => {
  wireCopyButtons(document);
  renderFinder();
  wireFinder();
  renderApiRef();
  wireLangToggle();

  // Keyboard shortcut: press "/" to focus finder
  document.addEventListener('keydown', e => {
    if (e.key === '/' && !(e.target.matches('input, textarea, button'))) {
      const s = document.getElementById('fnSearch');
      if (s) { e.preventDefault(); s.focus(); }
    }
    if (e.key === 'Escape') {
      const dlg = document.getElementById('fnModal');
      if (dlg?.open) dlg.close();
    }
  });
});
