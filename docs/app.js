// Evren viewer. Code written with help from an AI assistant.
(() => {
  const $ = (id) => document.getElementById(id);
  const cv = $('sky');
  const ctx = cv.getContext('2d');
  const S = {
    targets: [], meta: null, imgs: [], resid: [], median: null,
    k: 0, mode: 'frames', playTimer: null, blinkTimer: null, blinkPhase: false, blinkB: null,
    showKnown: true, showCand: true, smooth: true, focus: null, cursor: null, first: true,
  };
  const HELP = {
    frames: 'Each frame is one SPHEREx exposure, aligned so north is up and east is left. Stars stay put. Press Play to flip through the year, or Blink to flip between two dates.',
    changes: 'Each frame minus the static sky. White means something appeared, black means it was gone. A real mover shows up in every exposure of one visit.',
    static: 'The median of every frame: the part of the sky that never changes.',
  };
  const FONT = '-apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif';

  const fmtDate = (iso) => iso.slice(0, 16).replace('T', ' ') + ' UTC';
  const shortDate = (iso) => new Date(iso).toLocaleDateString('en-GB', { day: 'numeric', month: 'short', year: 'numeric', timeZone: 'UTC' });
  const stamp = (iso) => shortDate(iso) + ', ' + iso.slice(11, 16);
  const clsName = (c) => {
    c = c || '';
    if (c.startsWith('MB')) return 'main belt asteroid';
    if (c.startsWith('NEA')) return 'near Earth asteroid';
    if (c.startsWith('Mars')) return 'Mars crosser';
    if (c.startsWith('Hilda')) return 'Hilda asteroid';
    if (c.startsWith('Trojan')) return 'Jupiter trojan';
    if (c.startsWith('KBO') || c.startsWith('TNO') || c.startsWith('Centaur')) return 'outer solar system object';
    if (c.toLowerCase().includes('comet')) return 'comet';
    return 'small body';
  };

  // ------------------------------------------------------------- loading
  async function init() {
    try {
      S.targets = await (await fetch('data/targets.json')).json();
    } catch (e) {
      $('status').textContent = 'Could not load the data.';
      return;
    }
    const sel = $('field');
    const groups = new Map();
    for (const t of S.targets) {
      const g = t.group || 'Fields';
      if (!groups.has(g)) {
        const og = document.createElement('optgroup');
        og.label = g;
        groups.set(g, og);
        sel.appendChild(og);
      }
      const o = document.createElement('option');
      o.value = t.id;
      o.textContent = t.name;
      groups.get(g).appendChild(o);
    }
    sel.addEventListener('change', () => loadField(sel.value));
    const want = location.hash.slice(1);
    summarize();
    await loadField(S.targets.some((t) => t.id === want) ? want : S.targets[0].id);
  }

  const metaCache = {};
  const getMeta = (id) => (metaCache[id] ||= fetch(`data/${id}/meta.json`).then((r) => r.json()));

  // One honest line of totals across every field on the site.
  async function summarize() {
    const metas = await Promise.all(S.targets.map((t) => getMeta(t.id)));
    const frames = metas.reduce((a, m) => a + m.frames.length, 0);
    const plutoVisits = metas.reduce((a, m) => a + m.tracklets.filter((t) => /pluto/i.test(t.known || '')).length, 0);
    const known = new Set();
    metas.forEach((m) => m.known.forEach((e) => { if (e.recovered && !/pluto/i.test(e.name)) known.add(e.name); }));
    $('stats').innerHTML = `Searched <b>${frames}</b> SPHEREx frames in ${metas.length} fields. ` +
      `Pluto detected in <b>${plutoVisits}</b> separate visits. <b>${known.size}</b> known asteroids recovered where their orbits put them.`;
  }

  async function loadField(id) {
    stopAll();
    const meta = await getMeta(id);
    const base = `data/${id}/`;
    const total = meta.frames.length * 2 + 1;
    let done = 0;
    const tick = () => { done += 1; $('status').textContent = `Loading ${done} of ${total} images`; };
    const load = (src) => new Promise((ok, bad) => {
      const im = new Image();
      im.onload = () => { tick(); ok(im); };
      im.onerror = () => bad(new Error(src));
      im.src = base + src;
    });
    const [imgs, resid, median] = await Promise.all([
      Promise.all(meta.frames.map((f) => load(f.file))),
      Promise.all(meta.frames.map((f) => load(f.resid))),
      load('median.png'),
    ]);
    Object.assign(S, { meta, imgs, resid, median, k: 0, focus: null, blinkB: null });
    $('field').value = id;
    history.replaceState(null, '', '#' + id);
    $('fieldnote').textContent = meta.note || '';
    const sc = Math.max(4, Math.round(700 / meta.size_px));   // whole number scale keeps raw pixels square
    cv.width = cv.height = meta.size_px * sc;
    $('slider').max = meta.frames.length - 1;
    $('status').textContent = '';
    buildTicks();
    buildFindings();
    pickBlinkPair();
    if (S.first) {               // greet the visitor with the blink itself
      S.first = false;
      startBlink();
    }
  }

  // The classic view: the same stars on two dates. Prefer two visits of the
  // brightest known mover so the jump is easy to see.
  function pickBlinkPair() {
    const m = S.meta;
    const visitsOf = (frames) => {
      const out = [];
      for (const k of frames) {
        const last = out[out.length - 1];
        if (last && Math.abs(m.frames[k].mjd - m.frames[last[0]].mjd) < 0.5) last.push(k); else out.push([k]);
      }
      return out;
    };
    const movers = m.known.filter((e) => e.recovered)
      .map((e) => ({ e, ks: Object.keys(e.snr).filter((k) => e.snr[k] >= 8).map(Number).sort((a, b) => a - b) }))
      .filter((o) => visitsOf(o.ks).length >= 2)
      .sort((a, b) => Math.max(...Object.values(b.e.snr)) - Math.max(...Object.values(a.e.snr)));
    if (movers.length) {
      const { e, ks } = movers[0];
      const v = visitsOf(ks);
      const a = v[0][0], b = v[v.length - 1][0];
      const all = Object.keys(e.pos).map(Number).sort((x, y) => x - y);
      S.focus = { frames: all, xs: all.map((k) => e.pos[k][0]), ys: all.map((k) => e.pos[k][1]), name: e.name };
      S.blinkB = b;
      markFinding(e.name);
      setFrame(a);
      strip(all, S.focus.xs, S.focus.ys);
      return;
    }
    let b = 0, best = -1;
    m.frames.forEach((f, i) => { const d = Math.abs(f.mjd - m.frames[0].mjd); if (d > best) { best = d; b = i; } });
    S.blinkB = b;
    setFrame(0);
  }

  function setFrame(k) {
    const m = S.meta;
    if (!m) return;
    S.k = (k + m.frames.length) % m.frames.length;
    $('slider').value = S.k;
    renderInfo();
    markTick();
    blinkCaption();
    draw();
  }

  // ------------------------------------------------------------- drawing
  function shownFrame() {
    return S.blinkTimer && S.blinkPhase && S.blinkB != null ? S.blinkB : S.k;
  }

  function draw() {
    const m = S.meta;
    if (!m) return;
    const n = m.size_px, W = cv.width, s = W / n;
    ctx.imageSmoothingEnabled = S.smooth;
    ctx.imageSmoothingQuality = 'high';
    ctx.fillStyle = '#000';
    ctx.fillRect(0, 0, W, W);
    const k = shownFrame();
    const staticView = S.mode === 'static';
    const src = staticView ? S.median : (S.mode === 'changes' ? S.resid[k] : S.imgs[k]);
    ctx.drawImage(src, 0, 0, W, W);
    if (!staticView) {
      overlays(k, s);
      if (S.focus) focusRing(k, s);
    }
    if (S.blinkTimer) blinkLabel(k);
    compass(W, s);
  }

  function overlays(k, s) {
    const m = S.meta;
    ctx.lineWidth = 2;
    ctx.font = `600 14px ${FONT}`;
    if (S.showKnown) {
      for (const e of m.known) {
        if (!e.recovered) continue;
        const ks = Object.keys(e.pos).map(Number).sort((a, b) => a - b);
        if (ks.length > 1) {
          ctx.strokeStyle = 'rgba(124,199,255,0.35)';
          ctx.setLineDash([2, 5]);
          ctx.beginPath();
          ks.forEach((kk, i) => {
            const q = e.pos[String(kk)];
            ctx[i ? 'lineTo' : 'moveTo']((q[0] + 0.5) * s, (q[1] + 0.5) * s);
          });
          ctx.stroke();
          ctx.setLineDash([]);
        }
        const p = e.pos[String(k)];
        if (!p) continue;
        const x = (p[0] + 0.5) * s, y = (p[1] + 0.5) * s, r = Math.max(3.2 * s, 14);
        ctx.strokeStyle = ctx.fillStyle = '#7cc7ff';
        ctx.beginPath();
        ctx.arc(x, y, r, 0, Math.PI * 2);
        ctx.stroke();
        label(e.name, x, y - r * 0.6, r);
      }
    }
    if (S.showCand) {
      for (const t of m.tracklets) {
        if (t.known) continue;
        const j = t.frames.indexOf(k);
        if (j < 0) continue;
        const x = (t.x[j] + 0.5) * s, y = (t.y[j] + 0.5) * s, r = Math.max(3.2 * s, 14);
        const variable = t.kind === 'variable';
        ctx.strokeStyle = ctx.fillStyle = variable ? '#d9c26a' : '#ff5a1f';
        ctx.strokeRect(x - r, y - r, 2 * r, 2 * r);
        label((variable ? 'V' : 'C') + t.id, x, y + r, r);
      }
    }
  }

  // Text next to a marker, flipped to the left side near the right edge.
  function label(txt, x, y, r) {
    const w = ctx.measureText(txt).width;
    const lx = x + r + 6 + w < cv.width - 6 ? x + r + 6 : x - r - 6 - w;
    ctx.fillText(txt, lx, Math.max(y, 16));
  }

  function focusRing(k, s) {
    const f = S.focus;
    let j = f.frames.indexOf(k);
    const inside = j >= 0;
    if (!inside) j = 0;
    const x = (f.xs[j] + 0.5) * s, y = (f.ys[j] + 0.5) * s;
    ctx.strokeStyle = 'rgba(255,255,255,0.8)';
    ctx.lineWidth = 1.5;
    ctx.setLineDash(inside ? [] : [3, 5]);
    ctx.beginPath();
    ctx.arc(x, y, Math.max(7 * s, 30), 0, Math.PI * 2);
    ctx.stroke();
    ctx.setLineDash([]);
  }

  function blinkLabel(k) {
    const m = S.meta;
    const isB = S.blinkPhase && S.blinkB != null;
    const txt = `${isB ? 'B' : 'A'}  ${stamp(m.frames[k].t)}`;
    ctx.font = `600 16px ${FONT}`;
    const w = ctx.measureText(txt).width;
    ctx.fillStyle = 'rgba(0,0,0,0.6)';
    ctx.fillRect(12, 12, w + 20, 32);
    ctx.fillStyle = isB ? '#ff5a1f' : '#ffffff';
    ctx.fillText(txt, 22, 34);
  }

  function compass(W, s) {
    const m = S.meta;
    ctx.save();
    ctx.strokeStyle = ctx.fillStyle = 'rgba(255,255,255,0.75)';
    ctx.lineWidth = 1.5;
    ctx.font = `12px ${FONT}`;
    const ox = W - 26, oy = 50, L = 30;
    ctx.beginPath();
    ctx.moveTo(ox, oy); ctx.lineTo(ox, oy - L);
    ctx.moveTo(ox, oy); ctx.lineTo(ox - L, oy);
    ctx.stroke();
    ctx.fillText('N', ox - 4, oy - L - 6);
    ctx.fillText('E', ox - L - 13, oy + 4);
    const bar = (60 / m.pix_arcsec) * s;
    ctx.beginPath();
    ctx.moveTo(16, W - 18); ctx.lineTo(16 + bar, W - 18);
    ctx.stroke();
    ctx.fillText('1 arcmin', 16, W - 26);
    ctx.restore();
  }

  // --------------------------------------------------------------- info
  function renderInfo() {
    const m = S.meta, f = m.frames[S.k];
    const rows = [
      ['Frame', `${S.k + 1} of ${m.frames.length}`],
      ['Taken', fmtDate(f.t)],
      ['Band', `SPHEREx ${m.band}, ${f.um[0].toFixed(2)} to ${f.um[1].toFixed(2)} µm`],
      ['Field', `RA ${m.ra.toFixed(3)}°, Dec ${m.dec.toFixed(3)}°, ${Math.round(m.size_px * m.pix_arcsec / 60)} arcmin wide`],
      ['Exposure', `${f.obs_id}, ${f.collection === 'spherex_dr1' ? 'Data Release 1' : 'Quick Release 3'}`],
    ];
    if (f.noisy) rows.push(['Note', 'Stray light in this frame, so it was not searched']);
    if (S.cursor) rows.push(['Cursor', S.cursor]);
    $('info').innerHTML = rows.map(([a, b]) => `<dt>${a}</dt><dd>${b}</dd>`).join('');
  }

  function blinkCaption() {
    const m = S.meta;
    if (!m || S.blinkB == null) { $('blinkinfo').textContent = ''; return; }
    if (S.blinkB === S.k) {
      $('blinkinfo').textContent = `B is set to ${stamp(m.frames[S.k].t)}. Step to another date and press Blink.`;
      return;
    }
    $('blinkinfo').textContent = `Blink compares ${stamp(m.frames[S.k].t)} (A) with ${stamp(m.frames[S.blinkB].t)} (B).`;
  }

  // ------------------------------------------------------------ timeline
  function buildTicks() {
    const m = S.meta, svg = $('ticks');
    const t0 = m.frames[0].mjd, t1 = m.frames[m.frames.length - 1].mjd;
    const X = (mjd) => 10 + ((mjd - t0) / Math.max(t1 - t0, 1e-6)) * 580;
    let html = '<line x1="10" y1="26" x2="590" y2="26" stroke="#23272e" stroke-width="1"/>';
    const start = new Date((t0 - 40587) * 864e5), end = new Date((t1 - 40587) * 864e5);
    const d = new Date(Date.UTC(start.getUTCFullYear(), start.getUTCMonth() + 1, 1));
    while (d <= end) {
      const x = X(d.getTime() / 864e5 + 40587);
      const label = d.getUTCMonth() === 0 ? String(d.getUTCFullYear()) : d.toLocaleString('en-GB', { month: 'short', timeZone: 'UTC' });
      html += `<line x1="${x}" y1="30" x2="${x}" y2="34" stroke="#3a404a"/>`;
      html += `<text x="${x}" y="44" fill="#9aa0a8" font-size="9" text-anchor="middle">${label}</text>`;
      d.setUTCMonth(d.getUTCMonth() + 1);
    }
    m.frames.forEach((f, i) => {
      html += `<line class="tk" data-i="${i}" x1="${X(f.mjd)}" y1="14" x2="${X(f.mjd)}" y2="26" stroke="#6b7280" stroke-width="1.5"/>`;
    });
    html += '<line id="curb" x1="-10" y1="10" x2="-10" y2="28" stroke="#ff5a1f" stroke-width="1.5" stroke-dasharray="2 2"/>';
    html += '<line id="cur" x1="0" y1="6" x2="0" y2="28" stroke="#ffffff" stroke-width="2.5"/>';
    svg.innerHTML = html;
    svg.onclick = (ev) => {
      const r = svg.getBoundingClientRect();
      const x = ((ev.clientX - r.left) / r.width) * 600;
      let best = 0, bd = 1e9;
      m.frames.forEach((f, i) => { const dd = Math.abs(X(f.mjd) - x); if (dd < bd) { bd = dd; best = i; } });
      stopAll();
      setFrame(best);
    };
    S.tickX = X;
  }

  function markTick() {
    if (!S.tickX) return;
    const set = (id, k) => {
      const el = document.getElementById(id);
      if (!el) return;
      const x = k == null ? -10 : S.tickX(S.meta.frames[k].mjd);
      el.setAttribute('x1', x);
      el.setAttribute('x2', x);
    };
    set('cur', S.k);
    set('curb', S.blinkB);
    const fr = S.focus ? new Set(S.focus.frames) : null;
    for (const el of document.querySelectorAll('#ticks .tk')) {
      el.setAttribute('stroke', fr && fr.has(+el.dataset.i) ? '#7cc7ff' : '#6b7280');
    }
  }

  // ------------------------------------------------------------ findings
  function item(kind, title, detail, key) {
    const li = document.createElement('li');
    li.className = kind;
    if (key) li.dataset.key = key;
    li.innerHTML = `<div class="t">${title}</div>` + (detail ? `<div class="d">${detail}</div>` : '');
    return li;
  }

  function buildFindings() {
    const m = S.meta, ol = $('findings');
    ol.innerHTML = '';
    $('strip').innerHTML = '';
    const best = (e) => Math.round(Math.max(...Object.values(e.snr)));
    for (const e of m.known.filter((e) => e.recovered && /pluto/i.test(e.name))) {
      const li = item('known', 'Pluto, the original Planet X',
        `Seen in ${Object.values(e.snr).filter((v) => v >= 4).length} exposures, drifting across the stars exactly where its orbit says. Up to ${best(e)}σ`, e.name);
      li.onclick = () => focusKnown(e, li);
      ol.appendChild(li);
    }
    for (const t of m.tracklets.filter((t) => !t.known && t.kind === 'mover')) {
      const li = item('cand', `Candidate C${t.id}`,
        `${t.frames.length} exposures on ${shortDate(m.frames[t.frames[0]].t)}, up to ${Math.round(Math.max(...t.snr))}σ, no known object at this spot`);
      li.onclick = () => focusTracklet(t, li);
      ol.appendChild(li);
    }
    const groups = [];
    for (const t of m.tracklets.filter((t) => !t.known && t.kind === 'variable')) {
      const g = groups.find((g) => Math.hypot(g[0].x[0] - t.x[0], g[0].y[0] - t.y[0]) < 2);
      if (g) g.push(t); else groups.push([t]);
    }
    for (const g of groups) {
      const dates = [...new Set(g.map((t) => shortDate(m.frames[t.frames[0]].t)))];
      const li = item('var', `Changing spot V${g[0].id}`,
        `Brightened in place on ${dates.join(', ')}. It does not move, so it is not an asteroid. Possibly a variable star`);
      li.onclick = () => focusTracklet(g[0], li);
      ol.appendChild(li);
    }
    const rec = m.known.filter((e) => e.recovered && !/pluto/i.test(e.name)).sort((a, b) => (a.vmag ?? 99) - (b.vmag ?? 99));
    for (const e of rec) {
      const li = item('known', e.name,
        `Known ${clsName(e.cls)}, V ${e.vmag}. Found where its orbit puts it, ${best(e)}σ`, e.name);
      li.onclick = () => focusKnown(e, li);
      ol.appendChild(li);
    }
    if (!ol.children.length) ol.appendChild(item('empty', 'Nothing moved above the threshold here', 'Most of the sky looks like this: every frame was searched and nothing new turned up.'));
  }

  function select(li) {
    for (const el of document.querySelectorAll('#findings li')) el.classList.remove('sel');
    if (li) li.classList.add('sel');
  }

  function markFinding(key) {
    select([...document.querySelectorAll('#findings li')].find((el) => el.dataset.key === key));
  }

  function focusTracklet(t, li) {
    stopAll();
    select(li);
    S.focus = { frames: t.frames, xs: t.x, ys: t.y };
    const m = S.meta;
    let b = null, bd = -1;     // compare against the same spot on the farthest other date
    m.frames.forEach((f, i) => { const d = Math.abs(f.mjd - m.frames[t.frames[0]].mjd); if (d > 0.5 && d > bd) { bd = d; b = i; } });
    S.blinkB = b;
    setFrame(t.frames[0]);
    strip(t.frames, t.x, t.y);
  }

  function focusKnown(e, li) {
    stopAll();
    select(li);
    const ks = Object.keys(e.pos).map(Number).sort((a, b) => a - b);
    S.focus = { frames: ks, xs: ks.map((k) => e.pos[k][0]), ys: ks.map((k) => e.pos[k][1]), name: e.name };
    const good = Object.keys(e.snr).filter((k) => e.snr[k] >= 4).map(Number).sort((a, b) => a - b);
    const a = good.length ? good[0] : ks[0];
    const last = good.length ? good[good.length - 1] : ks[ks.length - 1];
    if (Math.abs(S.meta.frames[last].mjd - S.meta.frames[a].mjd) > 0.5) {
      S.blinkB = last;
    } else {
      let b = null, bd = -1;
      S.meta.frames.forEach((f, i) => { const d = Math.abs(f.mjd - S.meta.frames[a].mjd); if (d > 0.5 && d > bd) { bd = d; b = i; } });
      S.blinkB = b;
    }
    setFrame(a);
    strip(ks, S.focus.xs, S.focus.ys);
  }

  // Stamps of the object in its own exposures, plus the same spot on other dates.
  function strip(frames, xs, ys) {
    const m = S.meta, box = $('strip');
    box.innerHTML = '';
    const others = m.frames.map((f, i) => i).filter((i) => !frames.includes(i));
    const before = others.filter((i) => i < frames[0]).slice(-1);
    const after = others.filter((i) => i > frames[frames.length - 1]).slice(0, 1);
    const own = frames.length > 6 ? frames.filter((_, j) => j % Math.ceil(frames.length / 6) === 0) : frames;
    const list = [...before.map((i) => [i, false]), ...own.map((i) => [i, true]), ...after.map((i) => [i, false])];
    const half = 7;
    for (const [i, isOwn] of list) {
      const j = isOwn ? frames.indexOf(i) : 0;
      const cx = Math.round(xs[j]), cy = Math.round(ys[j]);
      const c = document.createElement('canvas');
      c.width = c.height = 2 * half + 1;
      const g = c.getContext('2d');
      g.imageSmoothingEnabled = false;
      g.fillStyle = '#000';
      g.fillRect(0, 0, c.width, c.height);
      g.drawImage(S.imgs[i], cx - half, cy - half, c.width, c.height, 0, 0, c.width, c.height);
      const fig = document.createElement('figure');
      const cap = document.createElement('figcaption');
      cap.textContent = (isOwn ? '' : 'Other date: ') + shortDate(m.frames[i].t);
      c.title = fmtDate(m.frames[i].t);
      if (isOwn) c.classList.add('own');
      c.onclick = () => { stopAll(); setFrame(i); };
      fig.append(c, cap);
      box.appendChild(fig);
    }
  }

  // ------------------------------------------------------------ playback
  function stopAll() {
    if (S.playTimer) { clearInterval(S.playTimer); S.playTimer = null; $('play').textContent = 'Play'; }
    if (S.blinkTimer) { clearInterval(S.blinkTimer); S.blinkTimer = null; S.blinkPhase = false; $('blink').classList.remove('on'); }
  }

  function togglePlay() {
    if (S.playTimer) { stopAll(); draw(); return; }
    stopAll();
    $('play').textContent = 'Pause';
    S.playTimer = setInterval(() => setFrame(S.k + 1), 1000 / +$('speed').value);
  }

  function startBlink() {
    stopAll();
    if (S.blinkB == null || S.blinkB === S.k) {
      const m = S.meta;
      let b = null, bd = -1;
      m.frames.forEach((f, i) => { const d = Math.abs(f.mjd - m.frames[S.k].mjd); if (d > bd) { bd = d; b = i; } });
      S.blinkB = b;
      markTick();
      blinkCaption();
    }
    $('blink').classList.add('on');
    S.blinkTimer = setInterval(() => { S.blinkPhase = !S.blinkPhase; draw(); }, 550);
    draw();
  }

  function toggleBlink() {
    if (S.blinkTimer) { stopAll(); draw(); return; }
    startBlink();
  }

  // -------------------------------------------------------------- events
  $('prev').onclick = () => { stopAll(); setFrame(S.k - 1); };
  $('next').onclick = () => { stopAll(); setFrame(S.k + 1); };
  $('play').onclick = togglePlay;
  $('blink').onclick = toggleBlink;
  $('setb').onclick = () => { stopAll(); S.blinkB = S.k; markTick(); blinkCaption(); draw(); };
  $('speed').oninput = () => {
    if (S.playTimer) { clearInterval(S.playTimer); S.playTimer = setInterval(() => setFrame(S.k + 1), 1000 / +$('speed').value); }
  };
  $('slider').oninput = (e) => { stopAll(); setFrame(+e.target.value); };
  $('showKnown').onchange = (e) => { S.showKnown = e.target.checked; draw(); };
  $('showCand').onchange = (e) => { S.showCand = e.target.checked; draw(); };
  $('rawPix').onchange = (e) => { S.smooth = !e.target.checked; cv.classList.toggle('raw', e.target.checked); draw(); };
  for (const b of document.querySelectorAll('.views button')) {
    b.onclick = () => {
      for (const o of document.querySelectorAll('.views button')) o.classList.remove('on');
      b.classList.add('on');
      S.mode = b.dataset.mode;
      $('modehelp').textContent = HELP[S.mode];
      draw();
    };
  }
  $('modehelp').textContent = HELP.frames;
  window.addEventListener('hashchange', () => {
    const id = location.hash.slice(1);
    if (S.meta && S.meta.id !== id && S.targets.some((t) => t.id === id)) loadField(id);
  });
  cv.addEventListener('mousemove', (ev) => {
    const m = S.meta;
    if (!m) return;
    const r = cv.getBoundingClientRect();
    const n = m.size_px;
    const px = ((ev.clientX - r.left) / r.width) * n - 0.5;
    const py = ((ev.clientY - r.top) / r.height) * n - 0.5;
    const c = (n - 1) / 2, d2r = Math.PI / 180;
    const dra = ((c - px) * m.pix_arcsec / 3600) / Math.cos(m.dec * d2r);
    const ddec = (c - py) * m.pix_arcsec / 3600;
    S.cursor = `RA ${(m.ra + dra).toFixed(4)}°, Dec ${(m.dec + ddec).toFixed(4)}°`;
    renderInfo();
  });
  cv.addEventListener('mouseleave', () => { S.cursor = null; if (S.meta) renderInfo(); });
  document.addEventListener('keydown', (ev) => {
    if (ev.target.tagName === 'SELECT' || ev.target.tagName === 'INPUT') return;
    if (ev.key === 'ArrowRight') { stopAll(); setFrame(S.k + 1); }
    else if (ev.key === 'ArrowLeft') { stopAll(); setFrame(S.k - 1); }
    else if (ev.key === ' ') { ev.preventDefault(); togglePlay(); }
    else if (ev.key === 'b' || ev.key === 'B') toggleBlink();
  });

  init();
})();
