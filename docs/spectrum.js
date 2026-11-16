// Evren spectrum panel. Code written with help from an AI assistant.
(() => {
  const box = document.getElementById('spectrum');
  if (!box) return;
  const src = box.dataset.src;
  const NS = 'http://www.w3.org/2000/svg';
  const el = (tag, attrs, parent) => {
    const n = document.createElementNS(NS, tag);
    for (const k in attrs) n.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(n);
    return n;
  };

  fetch(src + (src.includes('?') ? '&' : '?') + 'v=1').then(r => r.json()).then(draw).catch(() => {
    box.querySelector('.spec-chart').textContent = 'Spectrum not available.';
  });

  function draw(S) {
    const pts = S.points.filter(p => isFinite(p.r));
    const bins = S.binned || [];
    if (!pts.length) return;
    const W = 700, H = 330, m = { l: 52, r: 14, t: 18, b: 44 };
    const wMax = Math.max(...pts.map(p => p.w));
    const x0 = 0.7, x1 = Math.max(1.2, Math.ceil((wMax + 0.12) * 10) / 10);
    let yTop = 1.25;
    for (const b of bins) yTop = Math.max(yTop, b.r + b.e);
    yTop = Math.min(2.2, Math.ceil(yTop * 11) / 10);
    const X = w => m.l + (w - x0) / (x1 - x0) * (W - m.l - m.r);
    const Y = r => H - m.b - (r / yTop) * (H - m.t - m.b);

    const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': `${S.target} infrared spectrum` });
    const covered = w => pts.some(p => Math.abs(p.w - w) < 0.06);

    // methane ice absorption bands (laboratory band centers, microns)
    const bands = (box.dataset.bands || '').split(',').map(Number).filter(w => w && covered(w));
    for (const w of bands) {
      const half = Math.max(0.035, w / 60);
      el('rect', { x: X(w - half), y: m.t, width: X(w + half) - X(w - half), height: H - m.t - m.b, class: 'spec-band' }, svg);
    }
    const label = box.dataset.bandLabel;
    if (label && bands.length) {
      const w = bands[bands.length - 1];
      const t = el('text', { x: X(w), y: m.t + 12, class: 'spec-band-label', 'text-anchor': 'middle' }, svg);
      t.textContent = label;
    }

    // axes and grid
    for (let r = 0; r <= yTop + 1e-9; r += 0.5) {
      el('line', { x1: m.l, x2: W - m.r, y1: Y(r), y2: Y(r), class: 'spec-grid' }, svg);
      const t = el('text', { x: m.l - 8, y: Y(r) + 4, class: 'spec-tick', 'text-anchor': 'end' }, svg);
      t.textContent = r.toFixed(1);
    }
    const step = x1 - x0 > 2.5 ? 0.5 : 0.1;
    for (let w = Math.ceil(x0 / step) * step; w <= x1 + 1e-9; w += step) {
      el('line', { x1: X(w), x2: X(w), y1: H - m.b, y2: H - m.b + 5, class: 'spec-axis' }, svg);
      const t = el('text', { x: X(w), y: H - m.b + 19, class: 'spec-tick', 'text-anchor': 'middle' }, svg);
      t.textContent = (Math.round(w * 10) / 10).toString();
    }
    el('line', { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: 'spec-axis' }, svg);
    const xl = el('text', { x: (m.l + W - m.r) / 2, y: H - 6, class: 'spec-title', 'text-anchor': 'middle' }, svg);
    xl.textContent = 'Wavelength, microns';
    const yl = el('text', { x: 14, y: (m.t + H - m.b) / 2, class: 'spec-title', 'text-anchor': 'middle',
      transform: `rotate(-90 14 ${(m.t + H - m.b) / 2})` }, svg);
    yl.textContent = 'Reflected light';
    el('line', { x1: m.l, x2: W - m.r, y1: Y(1), y2: Y(1), class: 'spec-flat' }, svg);

    // single exposures
    for (const p of pts) {
      if (p.r < -0.2 || p.r > yTop) continue;
      el('circle', { cx: X(p.w), cy: Y(Math.max(0, p.r)), r: 1.7, class: 'spec-pt' }, svg);
    }
    // binned spectrum, broken across gaps in coverage
    let seg = [];
    const flush = () => {
      if (seg.length > 1) el('polyline', { points: seg.map(b => `${X(b.w)},${Y(Math.max(0, b.r))}`).join(' '), class: 'spec-line' }, svg);
      seg = [];
    };
    bins.forEach((b, i) => {
      if (i && b.w - bins[i - 1].w > 0.15) flush();
      seg.push(b);
    });
    flush();
    for (const b of bins) {
      if (b.n < 2) continue;
      const lo = Math.max(0, b.r - b.e), hi = Math.min(yTop, b.r + b.e);
      el('line', { x1: X(b.w), x2: X(b.w), y1: Y(lo), y2: Y(hi), class: 'spec-err' }, svg);
      el('circle', { cx: X(b.w), cy: Y(Math.max(0, b.r)), r: 3, class: 'spec-bin' }, svg);
    }
    const chart = box.querySelector('.spec-chart');
    chart.textContent = '';
    chart.appendChild(svg);

    // stamps
    const strip = box.querySelector('.spec-stamps');
    if (strip && S.stamps && S.stamps.length) {
      const base = src.replace(/[^/]*$/, '');
      for (const s of S.stamps) {
        const f = document.createElement('figure');
        const img = document.createElement('img');
        img.src = base + s.file;
        img.alt = `${S.target} at ${s.w.toFixed(2)} microns`;
        img.width = 168; img.height = 168;
        const c = document.createElement('figcaption');
        c.textContent = `${s.w.toFixed(2)} µm`;
        f.appendChild(img); f.appendChild(c);
        strip.appendChild(f);
      }
    }
  }
})();
